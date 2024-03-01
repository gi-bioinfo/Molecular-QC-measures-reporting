#!/usr/bin/env python3

import json
import os
import csv
import glob
from argparse import ArgumentParser
import sys
import subprocess
from collections import OrderedDict
import functools
import pandas as pd
import numpy as np
from datetime import date
import tarfile

pd.options.mode.chained_assignment = None  # default='warn'

total_size = {
    'wgs': 3088269832,
    'wxs': 148544048
}

variant_calling_stats_fields = {
    'study_id': 'study_id',
    'donor_id': 'donor_id',
    'submitter_donor_id': 'submitter_donor_id',
    'gender': 'gender',
    'geno_infer_gender_match': 'tumour.sanger.genotype_inference.frac_match_gender',
    'experimental_strategy': 'experimental_strategy',
    'normal_aligned': 'flags.normal_aligned',
    'tumour_aligned': 'flags.tumour_aligned',
    'sanger_called': 'flags.sanger_called',
    'mutect2_called': 'flags.mutect2_called',
    'normal_sample_id': 'normal.sample_id',
    'normal_error_rate': 'normal.alignment.error_rate',
    'normal_duplicate_rate': 'normal.alignment.duplicate_rate',
    'normal_pairs_on_different_chromosomes_rate': 'normal.alignment.pairs_on_different_chromosomes_rate',
    'normal_oxoQ_score': 'normal.alignment.oxoQ_score',
    'normal_insert_size_mean': 'normal.alignment.average_insert_size',
    'normal_avg_depth': 'normal.sanger.contamination.avg_depth',
    'normal_estimated_coverage': 'normal.alignment.estimated_coverage',
    'normal_sanger_contamination': 'normal.sanger.contamination.contamination',
    'normal_mutect2_contamination': 'normal.mutect2.contamination.contamination',
    'normal_properly_paired_reads': 'normal.alignment.properly_paired_reads',
    'tumour_sample_id': 'tumour.sample_id',
    'tumour_error_rate': 'tumour.alignment.error_rate',
    'tumour_duplicate_rate': 'tumour.alignment.duplicate_rate',
    'tumour_pairs_on_different_chromosomes_rate': 'tumour.alignment.pairs_on_different_chromosomes_rate',
    'tumour_oxoQ_score': 'tumour.alignment.oxoQ_score',
    'tumour_insert_size_mean': 'tumour.alignment.average_insert_size',
    'tumour_avg_depth': 'tumour.sanger.contamination.avg_depth',
    'tumour_estimated_coverage': 'tumour.alignment.estimated_coverage',
    'tumour_sanger_contamination': 'tumour.sanger.contamination.contamination',
    'tumour_mutect2_contamination': 'tumour.mutect2.contamination.contamination',
    'tumour_properly_paired_reads': 'tumour.alignment.properly_paired_reads',
    'ascat_normal_contamination': 'tumour.sanger.ascat_metrics.NormalContamination',
    'ascat_ploidy': 'tumour.sanger.ascat_metrics.Ploidy',
    'ascat_purity': 'tumour.sanger.ascat_metrics.rho',
    'mutect2_callable': 'tumour.mutect2.callable'
}

def get_extra_metrics(fname, extra_metrics, metrics):
    if not os.path.isfile(fname): 
        return metrics
    collected_sum_fields = {
        'insert size standard deviation': 'insert_size_sd'
    }
    
    unzip_dir = 'data/qc_metrics/unzip'
    if os.path.isdir(unzip_dir): 
        cmd = 'rm -rf %s && mkdir %s && tar -C %s -xzf %s' % (unzip_dir, unzip_dir, unzip_dir, fname)
    else:
        cmd = 'mkdir %s && tar -C %s -xzf %s' % (unzip_dir, unzip_dir, fname)
    run_cmd(cmd)

    for fn in glob.glob(os.path.join(unzip_dir, '*.aln.cram.bamstat')):    
        with open(fn, 'r') as f:
            for row in f:
                if not row.startswith('SN\t'): continue
                cols = row.replace(':', '').strip().split('\t')

                if not cols[1] in collected_sum_fields: continue
                metrics.update({
                    collected_sum_fields[cols[1]]: float(cols[2]) if ('.' in cols[2] or 'e' in cols[2]) else int(cols[2])
                    })
        
    return metrics

def get_extra_calling_metrics(fname):
    metrics = {}
    if not os.path.isfile(fname): 
        return metrics

    tar = tarfile.open(fname)
    for member in tar.getmembers():
        if member.name.endswith('.extra_info.json'):
            f = tar.extractfile(member)
            extra_info = json.load(f)
            metrics = extra_info.get('metrics')
            break
                
    return metrics

def report(donor, report_name):
    report_dir = os.path.dirname(report_name)
    if not os.path.exists(report_dir):
        os.makedirs(report_dir)
        
    keys = donor[0].keys()
    with open(report_name, 'w') as output_file:
        dict_writer = csv.DictWriter(output_file, keys, delimiter="\t")
        dict_writer.writeheader()
        dict_writer.writerows(donor)

def run_cmd(cmd):
    try:
        p = subprocess.run([cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           shell=True, check=True)

    except subprocess.CalledProcessError as e:   # this is triggered when cmd returned non-zero code
        print(e.stdout.decode("utf-8"))
        print('Execution returned non-zero code: %s. Additional error message: %s' %
              (e.returncode, e.stderr.decode("utf-8")), file=sys.stderr)
        sys.exit(e.returncode)

    except Exception as e:  # all other errors go here, like unable to find the command
        sys.exit('Execution failed: %s' % e)

    return p  # in case the caller of this funtion needs p.stdout, p.stderr etc


def udf(x, y):
    if not x:
        return
    if isinstance(x, dict):
        value = x.get(y)
    elif isinstance(x, list):
        value = [udf(a,y) for a in x]

    if value is not None: 
        return value
    else:
        return

def get_dict_value(fields, json_obj, field_map):
    tsv_obj = OrderedDict()
    if fields is None: fields = field_map.keys()
    for f in fields:
        fm = field_map.get(f, None)
        if fm is None: continue
        value = functools.reduce(udf, fm.split('.'), json_obj)
        tsv_obj[f] = value
    return tsv_obj 

def download(song_dump, file_type, ACCESSTOKEN, METADATA_URL, STORAGE_URL, include=None, subfolder=None):

    file_type_map = { # [analysisType, dataType, data_category]
        "qc_metrics": ['qc_metrics', ['Analysis QC', 'Sample QC'], 'Quality Control Metrics'],
        "timing_metrics": ['variant_calling_supplement', 'Variant Calling Supplement', None],
        "snv": ['variant_calling', 'Raw SNV Calls', 'Simple Nucleotide Variation'],
        "indel": ['variant_calling', 'Raw InDel Calls', 'Simple Nucleotide Variation'],
        "sv": ['variant_calling', 'Raw SV Calls', 'Structural Variation'],
        "cnv": ['variant_calling', 'Raw CNV Calls', 'Copy Number Variation']
    }

    if subfolder:
        data_dir = os.path.join("data", subfolder)
    else:
        data_dir = os.path.join("data", file_type_map[file_type][0])
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    downloaded = []
    for fn in glob.glob(os.path.join(data_dir, "*-*", "*.*"), recursive=True):
        downloaded.append(os.path.basename(fn))

    download_flist = set()
    with open(song_dump, 'r') as fp:
        for fline in fp:
            analysis = json.loads(fline)
            # print(analysis['analysisId'])
            if include and not analysis['analysisId'] in include: continue
            if not analysis.get('analysisState') == 'PUBLISHED': continue
            if not analysis['analysisType']['name'] == file_type_map[file_type][0]: continue
            output_dir = os.path.join(data_dir, analysis['studyId'])
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            for fl in analysis['files']:
                if not fl['dataType'] in file_type_map[file_type][1]: continue
                if file_type_map[file_type][2] is None and 'data_category' in fl['info']: continue
                if file_type_map[file_type][2] and not fl['info']['data_category'] == file_type_map[file_type][2]: continue
                download_flist.add(fl['fileName'])
                if fl['fileName'] in downloaded: continue

                cmd = 'export ACCESSTOKEN=%s && export METADATA_URL=%s \
                    && export STORAGE_URL=%s && export TRANSPORT_PARALLEL=3 \
                    && export TRANSPORT_MEMORY=8 \
                    && docker run -it --rm  -u $(id -u):$(id -g) \
                    -e ACCESSTOKEN -e METADATA_URL -e STORAGE_URL \
                    -e TRANSPORT_PARALLEL -e TRANSPORT_MEMORY \
                    -v "$PWD":"$PWD" -w "$PWD" overture/score:latest /score-client/bin/score-client download \
                    --study-id %s --object-id %s --output-dir %s' \
                    % (ACCESSTOKEN, METADATA_URL, STORAGE_URL, fl['studyId'], fl['objectId'], output_dir)

                run_cmd(cmd)
    return download_flist


def process_qc_metrics(song_dump, variant_calling_stats):
    sample_map = {}
    with open(song_dump, 'r') as fp:
        for fline in fp:
            analysis = json.loads(fline)
            if not analysis.get('analysisState') == 'PUBLISHED': continue
            if not analysis['samples'][0]['specimen']['tumourNormalDesignation'] == 'Tumour': continue
            studyId = analysis['studyId']
            sampleId = analysis['samples'][0]['sampleId']
            submitterSampleId = analysis['samples'][0]['submitterSampleId']
            matchedNormal = analysis['samples'][0]['matchedNormalSubmitterSampleId']
            experimental_strategy = analysis['experiment']['experimental_strategy'] if analysis['experiment'].get('experimental_strategy') else analysis['experiment']['library_strategy']
            normal_sample_id = '_'.join([studyId, experimental_strategy, matchedNormal])
            if not sample_map.get(normal_sample_id): 
                sample_map[normal_sample_id] = []
            sample_map[normal_sample_id].append(experimental_strategy+"_"+sampleId)

            donorId = analysis['samples'][0]['donor']['donorId']
            gender = analysis['samples'][0]['donor']['gender']
            
            unique_sampleId = experimental_strategy+"_"+sampleId
            if not variant_calling_stats.get(unique_sampleId): variant_calling_stats[unique_sampleId] = {
                'study_id': studyId,
                'donor_id': donorId,
                'submitter_donor_id': analysis['samples'][0]['donor']['submitterDonorId'],
                'gender': gender,
                'experimental_strategy': experimental_strategy,
                'flags': {
                    'normal_aligned': False,
                    'tumour_aligned': False,
                    'sanger_called': False,
                    'mutect2_called': False,
                    'open_filter': False
                },
                'normal': {
                    'alignment': {},
                    'sanger': {
                        'contamination': {}
                    },
                    'mutect2': {
                        'contamination': {}
                    }
                },
                'tumour': {
                    'sample_id': sampleId,
                    'submitterSampleId': submitterSampleId,
                    'alignment': {},
                    'sanger': {
                        'contamination': {},
                        'ascat_metrics': {},
                        'genotype_inference': {}
                    },
                    'mutect2': {
                        'contamination': {}
                    },
                    'open_filter_count': 0
                }
            }

            if analysis['analysisType']['name'] == 'variant_calling': 
                if analysis['workflow']['workflow_short_name'] in ['sanger-wgs', 'sanger-wxs']:
                    variant_calling_stats[unique_sampleId]['flags']['sanger_called'] = True
                if analysis['workflow']['workflow_short_name'] == 'gatk-mutect2':
                    variant_calling_stats[unique_sampleId]['flags']['mutect2_called'] = True
            elif analysis['analysisType']['name'] == 'sequencing_alignment':
                variant_calling_stats[unique_sampleId]['flags']['tumour_aligned'] = True
            elif analysis['analysisType']['name'] == 'variant_processing':
                open_filter_count = variant_calling_stats[unique_sampleId]['tumour']['open_filter_count'] + 1
                if open_filter_count == 4: 
                  variant_calling_stats[unique_sampleId]['flags']['open_filter'] = True
                variant_calling_stats[unique_sampleId]['tumour']['open_filter_count'] = open_filter_count
            elif not analysis['analysisType']['name'] == 'qc_metrics': 
                continue
            
            for fl in analysis['files']:
                if fl.get('info') and fl['info'].get('data_subtypes') and 'Cross Sample Contamination' in fl['info']['data_subtypes']:
                    fname = os.path.join("data", 'qc_metrics', analysis['studyId'], fl['fileName'])
                    metrics = get_extra_calling_metrics(fname)
                    if metrics['sample_id'] == sampleId:
                        if 'sanger' in fl['fileName']:
                            variant_calling_stats[unique_sampleId]['tumour']['sanger']['contamination'].update(metrics)
                        elif 'gatk-mutect2' in fl['fileName']:
                            variant_calling_stats[unique_sampleId]['tumour']['mutect2']['contamination'].update(metrics)
                        else:
                            pass
                    else:
                        if 'sanger' in fl['fileName']:
                            variant_calling_stats[unique_sampleId]['normal']['sanger']['contamination'].update(metrics)
                        elif 'gatk-mutect2' in fl['fileName']:
                            variant_calling_stats[unique_sampleId]['normal']['mutect2']['contamination'].update(metrics)
                        else:
                            pass
                elif fl.get('info') and fl['info'].get('data_subtypes') and 'Ploidy' in fl['info']['data_subtypes'] and 'Tumour Purity' in fl['info']['data_subtypes']:
                    fname = os.path.join("data", 'qc_metrics', analysis['studyId'], fl['fileName'])
                    metrics = get_extra_calling_metrics(fname)
                    variant_calling_stats[unique_sampleId]['tumour']['sanger']['ascat_metrics'].update(metrics)
                elif fl.get('info') and fl['info'].get('data_subtypes') and 'Genotyping Stats' in fl['info']['data_subtypes']:
                    fname = os.path.join("data", 'qc_metrics', analysis['studyId'], fl['fileName'])
                    metrics = get_extra_calling_metrics(fname)
                    variant_calling_stats[unique_sampleId]['tumour']['sanger']['genotype_inference'].update(metrics['tumours'][0]['gender'])
                elif fl.get('info') and fl['info'].get('data_subtypes') and 'Alignment Metrics' in fl['info']['data_subtypes'] and 'qc_metrics' in fl['fileName']:
                    metrics = {}
                    for fn in ['error_rate', 'properly_paired_reads', 'total_reads', 'average_insert_size', 'average_length', 'pairs_on_different_chromosomes']:
                        metrics.update({fn: fl['info']['metrics'][fn]})
                    if fl['info']['metrics']['total_reads']==0: continue
                    metrics.update({
                        'duplicate_rate': round(fl['info']['metrics']['duplicated_bases']/(fl['info']['metrics']['total_reads']*fl['info']['metrics']['average_length']), 3)
                        })
                    metrics.update({
                        'pairs_on_different_chromosomes_rate': round(fl['info']['metrics']['pairs_on_different_chromosomes']*2/(fl['info']['metrics']['paired_reads']), 3)
                    })
                    metrics.update({
                        'estimated_coverage': round(fl['info']['metrics']['total_bases']/total_size.get(experimental_strategy.lower()), 3)
                    })
                    fname = os.path.join("data", 'qc_metrics', analysis['studyId'], fl['fileName'])
                    extra_metrics = ['insert_size_sd']
                    metrics = get_extra_metrics(fname, extra_metrics, metrics)
                    variant_calling_stats[unique_sampleId]['tumour']['alignment'].update(metrics)
                elif fl.get('info') and fl['info'].get('data_subtypes') and 'OxoG Metrics' in fl['info']['data_subtypes']:
                    variant_calling_stats[unique_sampleId]['tumour']['alignment'].update({'oxoQ_score': fl['info']['metrics']['oxoQ_score'] if fl['info']['metrics'].get('oxoQ_score') else None})
                    
                elif fl.get('info') and fl['info'].get('data_subtypes') and 'Variant Callable Stats' in fl['info']['data_subtypes']:
                    fname = os.path.join("data", 'qc_metrics', analysis['studyId'], fl['fileName'])
                    metrics = get_extra_calling_metrics(fname)
                    variant_calling_stats[unique_sampleId]['tumour']['mutect2'].update(metrics)

                elif fl['dataType'] == 'Aligned Reads':
                    variant_calling_stats[unique_sampleId]['tumour']['alignment'].update({"file_size": round(fl['fileSize']/(1024*1024*1024), 3)})

                else:
                    continue
                


    with open(song_dump, 'r') as fp:
        for fline in fp:
            analysis = json.loads(fline)
            if not analysis.get('analysisState') == 'PUBLISHED': continue
            if not analysis['analysisType']['name'] in ['qc_metrics', 'sequencing_alignment']: continue
            if not analysis['samples'][0]['specimen']['tumourNormalDesignation'] == 'Normal': continue

            experimental_strategy = analysis['experiment']['experimental_strategy'] if analysis['experiment'].get('experimental_strategy') else analysis['experiment']['library_strategy']           
            studyId = analysis['studyId']
            submitSampleId = analysis['samples'][0]['submitterSampleId']
            normal_sample_id = '_'.join([studyId, experimental_strategy, submitSampleId])
            
            if not normal_sample_id in sample_map: continue
            
            for fl in analysis['files']:
                if fl.get('info') and fl['info'].get('data_subtypes') and 'Alignment Metrics' in fl['info']['data_subtypes'] and 'qc_metrics' in fl['fileName']:
                    metrics = {}
                    for fn in ['error_rate', 'properly_paired_reads', 'total_reads', 'average_insert_size', 'average_length', 'pairs_on_different_chromosomes']:
                        metrics.update({fn: fl['info']['metrics'][fn]})
                    if fl['info']['metrics']['total_reads'] == 0: continue    
                    metrics.update({
                        'duplicate_rate': round(fl['info']['metrics']['duplicated_bases']/(fl['info']['metrics']['total_reads']*fl['info']['metrics']['average_length']), 3)
                        })
                    if fl['info']['metrics']['paired_reads']>0:
                        metrics.update({
                            'pairs_on_different_chromosomes_rate': round(fl['info']['metrics']['pairs_on_different_chromosomes']*2/(fl['info']['metrics']['paired_reads']), 3)
                        })
                    metrics.update({
                        'estimated_coverage': round(fl['info']['metrics']['mapped_bases_cigar']/total_size.get(experimental_strategy.lower()), 3)
                    })
                    fname = os.path.join("data", 'qc_metrics', analysis['studyId'], fl['fileName'])
                    extra_metrics = ['insert_size_sd']
                    metrics = get_extra_metrics(fname, extra_metrics, metrics)

                    for sa in sample_map[normal_sample_id]:
                        variant_calling_stats[sa]['normal']['sample_id'] = analysis['samples'][0]['sampleId']
                        variant_calling_stats[sa]['normal']['submitterSampleId'] = analysis['samples'][0]['submitterSampleId']  
                        variant_calling_stats[sa]['normal']['alignment'].update(metrics)
                        variant_calling_stats[sa]['flags']['normal_aligned'] = True 
                elif fl.get('info') and fl['info'].get('data_subtypes') and 'OxoG Metrics' in fl['info']['data_subtypes']:
                    for sa in sample_map[normal_sample_id]:  
                        variant_calling_stats[sa]['normal']['alignment'].update({'oxoQ_score': fl['info']['metrics'].get('oxoQ_score', None)})                 
                elif fl['dataType'] == 'Aligned Reads':
                    for sa in sample_map[normal_sample_id]:  
                        variant_calling_stats[sa]['normal']['alignment'].update({"file_size": round(fl['fileSize']/(1024*1024*1024), 3)})                    
                else:
                    continue                   

    return variant_calling_stats


def main():
    parser = ArgumentParser()
    parser.add_argument("-d", "--dump_path", dest="dump_path", type=str, default="data/rdpc-song.jsonl", help="path to song dump jsonl file")
    parser.add_argument("-m", "--metadata_url", dest="metadata_url", type=str, default="https://song.rdpc-prod.cumulus.genomeinformatics.org")
    parser.add_argument("-s", "--storage_url", dest="storage_url", type=str, default="https://score.rdpc-prod.cumulus.genomeinformatics.org")
    parser.add_argument("-t", "--token", dest="token", type=str, required=True)
    args = parser.parse_args()

    song_dump = args.dump_path
    variant_calling_stats = {}

    #download qc_metrics
    download(song_dump, 'qc_metrics', args.token, args.metadata_url, args.storage_url)

    variant_calling_stats = process_qc_metrics(song_dump, variant_calling_stats)

    report_dir = 'report'
    if not os.path.exists(report_dir):
        os.makedirs(report_dir)
    study_id = args.dump_path.split('.')[-3]
    with open(os.path.join(report_dir, study_id+'.variant_calling_stats.json'), 'w') as f:
        f.write(json.dumps(variant_calling_stats, indent=2))

    # generate tsv file
    date_str = date.today().strftime("%Y-%m-%d")
    variant_calling_stats_tsv = []

    for d, v in variant_calling_stats.items():
        variant_calling_stats_tsv.append(get_dict_value(None, v, variant_calling_stats_fields))
    report(variant_calling_stats_tsv, os.path.join(report_dir, '.'.join([study_id, date_str, 'qc.tsv'])))


if __name__ == "__main__":
    main()
