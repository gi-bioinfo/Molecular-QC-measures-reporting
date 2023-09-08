#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
  Copyright (C) 2022,  icgc-argo

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Authors:
    Edmund Su
"""

import pandas as pd
import requests
import numpy as np
import os
import argparse
import plotly
from plotly.subplots import make_subplots
import plotly.graph_objs as go
import pickle
import sys

def main():
    """

    """
    parser = argparse.ArgumentParser(description='Retrieve stats from SONG API and generate plots')
    parser.add_argument('-p', '--project', dest="project", help="projects to query", required=True,type=str,nargs="+")
    parser.add_argument('-u', '--url', dest="rdpc_url", help="SONG RDPC URL", required=True,type=str)
    parser.add_argument('-o', '--output_directory', dest="out_dir", help="SONG RDPC URL", default=os.getcwd(),type=str)
    parser.add_argument('-x', '--exclude_analyses', dest="excluded_analyses", help="analyses to exclude",nargs="+",type=str)
    parser.add_argument('-e', '--experiment', dest="experiment", help="experiment type", required=True,nargs="+",type=str,choices=['RNA-Seq', 'WGS', 'WXS'])
    parser.add_argument('-z', '--plot', dest="plot", help="make pretty plots", default=True,type=bool)
    parser.add_argument('-d', '--debug', dest="debug", help="debug", default=False,type=bool)
    parser.add_argument('-s', '--state', dest="state", help="analysis state to query : PUBLISHED,SUPPRESSED,UNPUBLISHED",
                        nargs="+",
                        default=['PUBLISHED'],
                        choices=["PUBLISHED","SUPPRESSED","UNPUBLISHED"],
                        type=str)

    cli_input= parser.parse_args()

    metadata={}
    for project in cli_input.project:
        metadata[project]={}
        for state in cli_input.state:
            response=song_phone_home(project,cli_input.rdpc_url,state)
            for experiment in cli_input.experiment:
                write_dir="%s/%s_%s_%s" % (cli_input.out_dir,state,project,experiment)
                if not os.path.exists(write_dir):
                    os.mkdir(write_dir)

                tsv_dir="%s/%s" % (write_dir,"tsv")
                if not os.path.exists(tsv_dir):
                    os.mkdir(tsv_dir)

                metadata[project][experiment]={}
                
                metadata[project][experiment][state]=\
                generate_rdpc_metadata(response,experiment)
                
                metadata[project][experiment][state].iloc[:,:-1].to_csv(
                    "%s/%s_%s_%s_fileIDs.tsv" % (tsv_dir,state,project,experiment),
                    sep="\t"
                )
                plots={}
                metrics={}
                ###RNA-Seq
                if experiment=='RNA-Seq':
                    metrics['Picard:CollectRnaSeqMetrics']=aggreate_picard_collect_rnaseq_metrics(
                        response,
                        metadata[project][experiment][state],
                        cli_input.excluded_analyses,
                        cli_input.debug
                    )
                    metrics['biobambam2:bammarkduplicates2']=aggregate_picard_mark_duplicates_metrics(
                        response,
                        metadata[project][experiment][state],
                        cli_input.excluded_analyses,
                        cli_input.debug
                    )                    
                
                    
                    for ind,item in enumerate(['TOTAL_READS','DUPLICATION_PCT','MAPPING_PCT']):
                        title="%s %s %s" % (project,experiment,item)
                        plots["fig.%s.%s.%s" % (1,ind+1,title.replace(" ","_"))]=generate_plot(
                            metrics['biobambam2:bammarkduplicates2'],
                            1000,
                            600,
                            ["STAR","HISAT2"],
                            [item],
                            title
                        )
                    for ind,item in enumerate([
                        "median_3prime_bias",
                        "median_5prime_bias",
                        "median_5prime_to_3prime_bias",
                        "median_cv_coverage",
                        "pct_coding_bases",
                        "pct_correct_strand_reads",
                        "pct_intergenic_bases",
                        "pct_intronic_bases",
                        "pct_mrna_bases",
                        "pct_r1_transcript_strand_reads",
                        "pct_r2_transcript_strand_reads",
                        "pct_ribosomal_bases",
                        "pct_usable_bases",
                        "pct_utr_bases"]):
                        title="%s %s %s" % (project,experiment,item)
                        plots["fig.%s.%s.%s" % (2,ind+1,title.replace(" ","_"))]=generate_plot(
                            metrics['Picard:CollectRnaSeqMetrics'],
                            1000,
                            600,
                            ["STAR","HISAT2"],
                            [item],
                            title
                        )

                    metrics['Picard:CollectRnaSeqMetrics'].to_csv("%s/%s_%s_%s_rnaMetrics.tsv" % (tsv_dir,state,project,experiment),sep="\t")
                    metrics['biobambam2:bammarkduplicates2'].to_csv("%s/%s_%s_%s_libraryMetrics.tsv" % (tsv_dir,state,project,experiment),sep="\t")

                    save_pkl_plots(write_dir,plots,cli_input.plot)
                if experiment=='WGS' or experiment=='WXS':
                    metrics['biobambam2:bammarkduplicates2']=aggregate_picard_mark_duplicates_metrics(
                        response,
                        cli_input.excluded_analyses,
                        cli_input.debug
                    )
                    metrics['GATK:CollectOxoGMetrics']=aggregate_gatk_oxo_metrics(
                        response,
                        cli_input.excluded_analyses,
                        cli_input.debug
                    )
                    metrics['Samtools:stats']=aggregate_samtools_stats_metrics(
                        response,
                        cli_input.excluded_analyses,
                        cli_input.debug
                    )
                    metrics['Picard:CollectQualityYieldMetrics']=aggregate_gatk_quality_yield_metrics(
                        response,
                        cli_input.excluded_analyses,
                        cli_input.debug
                    )
                    metrics['Sanger:verifyBamHomChk']=aggregate_sanger_verifyBamHomChk_metrics(
                        response,
                        cli_input.excluded_analyses,
                        cli_input.debug
                    )
                    metrics['Sanger:compareBamGenotypes']=aggregate_sanger_compareBamGenotypes_metrics(
                        response,
                        cli_input.excluded_analyses,
                        cli_input.debug
                    )

                    for ind,item in enumerate(['TOTAL_READS','DUPLICATION_PCT','MAPPING_PCT']):
                        title="%s %s %s" % (project,experiment,item)
                        plots["fig.%s.%s.%s" % (1,ind+1,title.replace(" ","_"))]=generate_plot(
                            metrics['biobambam2:bammarkduplicates2'],
                            1000,
                            600,
                            ["BWA-MEM"],
                            [item],
                            title
                        )
                    for ind,item in enumerate(['oxoQ_score']):
                        title="%s %s %s" % (project,experiment,item)
                        plots["fig.%s.%s.%s" % (1,ind+1,title.replace(" ","_"))]=generate_plot(
                            metrics['GATK:CollectOxoGMetrics'],
                            1000,
                            600,
                            ["BWA-MEM"],
                            [item],
                            title
                        )
                    for ind,item in enumerate(["average_insert_size",
                        "average_length",
                        "duplicated_bases",
                        "error_rate",
                        "mapped_bases_cigar",
                        "mapped_reads",
                        "mismatch_bases",
                        "paired_reads",
                        "pairs_on_different_chromosomes",
                        "properly_paired_reads",
                        "total_bases",
                        "total_reads"]):
                        title="%s %s %s" % (project,experiment,item)
                        plots["fig.%s.%s.%s" % (1,ind+1,title.replace(" ","_"))]=generate_plot(
                            metrics['Samtools:stats'],
                            1000,
                            600,
                            ["BWA-MEM"],
                            [item],
                            title
                        )
                    for ind,item in enumerate(['total_reads','read_length','pf_reads']):
                        title="%s %s %s" % (project,experiment,item)
                        plots["fig.%s.%s.%s" % (1,ind+1,title.replace(" ","_"))]=generate_plot(
                            metrics['Picard:CollectQualityYieldMetrics'],
                            1000,
                            600,
                            ["BWA-MEM"],
                            [item],
                            title
                        )

                    for ind,item in enumerate(['avg_depth','contamination','reads_used','snps_used']):
                        title="%s %s %s" % (project,experiment,item)
                        plots["fig.%s.%s.%s" % (1,ind+1,title.replace(" ","_"))]=generate_plot(
                            metrics['Sanger:verifyBamHomChk'],
                            1000,
                            600,
                            ["BWA-MEM"],
                            [item],
                            title
                        )

                    for ind,item in enumerate(['total_loci_genotype','frac_match_gender','frac_informative_genotype','frac_matched_genotype']):
                        title="%s %s %s" % (project,experiment,item)
                        plots["fig.%s.%s.%s" % (1,ind+1,title.replace(" ","_"))]=generate_plot(
                            metrics['Sanger:compareBamGenotypes'],
                            1000,
                            600,
                            ["BWA-MEM"],
                            [item],
                            title
                        )
                    metrics['biobambam2:bammarkduplicates2'].to_csv("%s/%s_%s_%s_markDupMetrics.tsv" % (tsv_dir,state,project,experiment),sep="\t")
                    metrics['GATK:CollectOxoGMetrics'].to_csv("%s/%s_%s_%s_oxoMetrics.tsv" % (tsv_dir,state,project,experiment),sep="\t")
                    metrics['Samtools:stats'].to_csv("%s/%s_%s_%s_samtoolsMetrics.tsv" % (tsv_dir,state,project,experiment),sep="\t")
                    metrics['Picard:CollectQualityYieldMetrics'].to_csv("%s/%s_%s_%s_readGroupMetrics.tsv" % (tsv_dir,state,project,experiment),sep="\t")
                    metrics['Sanger:verifyBamHomChk'].to_csv("%s/%s_%s_%s_verifyBamMetrics.tsv" % (tsv_dir,state,project,experiment),sep="\t")
                    metrics['Sanger:compareBamGenotypes'].to_csv("%s/%s_%s_%s_bamGenotypesMetrics.tsv" % (tsv_dir,state,project,experiment),sep="\t")
                    save_pkl_plots(write_dir,plots,cli_input.plot)   

def save_pkl_plots(out_dir,gen_plots,plot):
    print("Saving plots...")
    svg_dir="%s/%s" % (out_dir,"svg")
    pkl_dir="%s/%s" % (out_dir,"pkl")

    if not os.path.exists(svg_dir):
        os.mkdir(svg_dir)
    if not os.path.exists(pkl_dir):
        os.mkdir(pkl_dir)
    for gen_plot in gen_plots.keys():
        file = open("%s/%s.pkl" % (pkl_dir,gen_plot),"wb")
        pickle.dump(gen_plots[gen_plot],file)
        file.close()
    print("Saving plots...Complete")

    if plot:
        print("Saving plots SVGs...")
        for gen_plot in gen_plots.keys():
            gen_plots[gen_plot].write_image("%s/%s.svg" % (svg_dir,gen_plot))
        print("Saving plots SVGs...Complete")

def generate_plot(metrics,x_dim,y_dim,cols,rows,title):
    print("Generating plot for %s" % (title))
    fig=plotly.subplots.make_subplots(
        cols=len(cols),
        rows=len(rows),
        subplot_titles=cols
    )
    
    for row_ind,row in enumerate(rows):
        for col_ind,col in enumerate(cols):
            fig.append_trace(
                go.Scatter(
                    x=metrics.query("PIPELINE==@col").sort_values(row)['sampleId'].values.tolist(),
                    y=metrics.query("PIPELINE==@col").sort_values(row)[row].values.tolist(),
                    mode='markers+lines',
                    showlegend=False),
                row_ind+1,
                col_ind+1
            )
            
            for val in [25,5,75]:
                fig.append_trace(
                    go.Scatter(
                        x=[
                            metrics.query("PIPELINE==@col").sort_values(row)['sampleId'].values.tolist()[0],
                            metrics.query("PIPELINE==@col").sort_values(row)['sampleId'].values.tolist()[-1]
                        ],
                        y=[np.percentile(metrics.query("PIPELINE==@col").sort_values(row)[row].values.tolist(),val,axis=0)]*2,
                        mode='lines',
                        line=dict(dash="dash",color="black"),
                        opacity=0.3,
                        showlegend=False),
                    row_ind+1,
                    col_ind+1
                )
                
    fig['layout'].update(
        width=x_dim,
        height=y_dim,
        title=title,
        xaxis=dict(title="analysisId"),
        yaxis=dict(title=""),
        showlegend=True,
        titlefont=dict(size=20)
    )
    return(fig)
def aggregate_gatk_quality_yield_metrics(response,analysis_exclude_list,debug):
    print("Aggregating metrics from : %s" % ('Picard:CollectQualityYieldMetrics'))
    metrics=pd.DataFrame()
    debug=False
    analysis_exclude_list=None
    for count,analysis in enumerate(response.json()):
        if analysis_exclude_list!=None:
            if analysis['analysisId'] in analysis_exclude_list:
                continue
        if count%50==0 and debug:
            print(count)
        for file in analysis['files']:
            if file['info'].get('analysis_tools'):
                if file['info']['analysis_tools'][0]=='Picard:CollectQualityYieldMetrics':
                    
                    sampleId=analysis['samples'][0]['sampleId']
                    readGroupId=file['info']['metrics']['read_group_id']
                    uniqueId=sampleId+"."+readGroupId
                    
                    if "star" in file['fileName']:
                        metrics.loc[uniqueId,"PIPELINE"]="STAR"
                    elif "hisat2" in file['fileName']:
                        metrics.loc[uniqueId,"PIPELINE"]="HISAT2"
                    else :
                        metrics.loc[uniqueId,"PIPELINE"]="BWA-MEM"
                        

                    metrics.loc[uniqueId,"sampleId"]=analysis['samples'][0]['sampleId']
                    metrics.loc[uniqueId,"readGroupId"]=file['info']['metrics']['read_group_id']
                    metrics.loc[uniqueId,"analysisId"]=analysis['analysisId']
                    metrics.loc[uniqueId,"objectId"]=file['objectId']
                    metrics.loc[uniqueId,'total_reads']=file['info']['metrics']['total_reads']
                    metrics.loc[uniqueId,'read_length']=file['info']['metrics']['read_length']
                    metrics.loc[uniqueId,'pf_reads']=file['info']['metrics']['pf_reads']
            else:
                continue

    return(metrics)

def aggregate_samtools_stats_metrics(response,analysis_exclude_list,debug):
    print("Aggregating metrics from : %s" % ('Samtools:stats'))
    metrics=pd.DataFrame()
    #total=len([response.json()[int(ind)] for ind in metadata_df.query("analysisId!=@analysis_exclude_list")["ind"].values.tolist()])
    for count,analysis in enumerate(response.json()):
        if analysis_exclude_list!=None:
            if analysis['analysisId'] in analysis_exclude_list:
                continue
        if count%50==0 and debug:
            print(count)
        for file in analysis['files']:
            if file['info'].get('analysis_tools'):
                if file['info']['analysis_tools'][0]=='Samtools:stats':
                    metrics.loc[analysis['analysisId'],'sampleId']=analysis['samples'][0]['sampleId']
                    if "star" in file['fileName']:
                        metrics.loc[analysis['analysisId'],"PIPELINE"]="STAR"
                    elif "hisat2" in file['fileName']:
                        metrics.loc[analysis['analysisId'],"PIPELINE"]="HISAT2"
                    else :
                        metrics.loc[analysis['analysisId'],"PIPELINE"]="BWA-MEM"
                    for query in [
                        "average_insert_size",
                        "average_length",
                        "duplicated_bases",
                        "error_rate",
                        "mapped_bases_cigar",
                        "mapped_reads",
                        "mismatch_bases",
                        "paired_reads",
                        "pairs_on_different_chromosomes",
                        "properly_paired_reads",
                        "total_bases",
                        "total_reads",
                    ]:
                        metrics.loc[analysis['analysisId'],query]=file['info']['metrics'][query]
            else:
                continue

    return(metrics)

def aggregate_sanger_compareBamGenotypes_metrics(response,analysis_exclude_list,debug):
    metrics=pd.DataFrame()
    debug=False
    analysis_exclude_list=None
    for count,analysis in enumerate(response.json()):
            if analysis_exclude_list!=None:
                if analysis['analysisId'] in analysis_exclude_list:
                    continue
            if count%50==0 and debug:
                print(count)
            for file in analysis['files']:
                if file['info'].get('analysis_tools'):
                    if file['info']['analysis_tools'][0]=='Sanger:compareBamGenotypes':
                        if file['info'].get('metrics'):

                            if "star" in file['fileName']:
                                metrics.loc[analysis['analysisId'],"PIPELINE"]="STAR"
                            elif "hisat2" in file['fileName']:
                                metrics.loc[analysis['analysisId'],"PIPELINE"]="HISAT2"
                            else :
                                metrics.loc[analysis['analysisId'],"PIPELINE"]="BWA-MEM"

                            metrics.loc[analysis['analysisId'],"sampleId"]=analysis['samples'][0]['sampleId']
                            metrics.loc[analysis['analysisId'],"compared_against"]=file['info']['metrics']['compared_against']
                            metrics.loc[analysis['analysisId'],"compared_against"]=file['info']['metrics']['total_loci_gender']
                            metrics.loc[analysis['analysisId'],"total_loci_genotype"]=file['info']['metrics']['total_loci_genotype']
                            metrics.loc[analysis['analysisId'],'frac_match_gender']=file['info']['metrics']['tumours'][0]['gender']['frac_match_gender']
                            metrics.loc[analysis['analysisId'],'gender']=file['info']['metrics']['tumours'][0]['gender']['gender']
                            metrics.loc[analysis['analysisId'],'frac_informative_genotype']=file['info']['metrics']['tumours'][0]['genotype']['frac_informative_genotype']
                            metrics.loc[analysis['analysisId'],'frac_matched_genotype']=file['info']['metrics']['tumours'][0]['genotype']['frac_matched_genotype']
            else:
                continue

    return(metrics)

def aggregate_sanger_verifyBamHomChk_metrics(response,analysis_exclude_list,debug):
    metrics=pd.DataFrame()
    debug=False
    analysis_exclude_list=None
    for count,analysis in enumerate(response.json()):
            if analysis_exclude_list!=None:
                if analysis['analysisId'] in analysis_exclude_list:
                    continue
            if count%50==0 and debug:
                print(count)
            for file in analysis['files']:
                if file['info'].get('analysis_tools'):
                    if file['info']['analysis_tools'][0]=='Sanger:verifyBamHomChk':
                        if file['info'].get('metrics'):
                            metrics.loc[analysis['analysisId'],"sampleId"]=analysis['samples'][0]['sampleId']
                            if "star" in file['fileName']:
                                metrics.loc[analysis['analysisId'],"PIPELINE"]="STAR"
                            elif "hisat2" in file['fileName']:
                                metrics.loc[analysis['analysisId'],"PIPELINE"]="HISAT2"
                            else :
                                metrics.loc[analysis['analysisId'],"PIPELINE"]="BWA-MEM"
                                
                            for key in file['info']['metrics'].keys():
                                if key=='sample_id':
                                    continue
                                else:
                                    metrics.loc[analysis['analysisId'],key]=file['info']['metrics'][key]
            else:
                continue

    return(metrics)

def aggregate_gatk_oxo_metrics(response,analysis_exclude_list,debug):
    print("Aggregating metrics from : %s" % ('GATK:CollectOxoGMetrics'))
    metrics=pd.DataFrame()
    #total=len([response.json()[int(ind)] for ind in metadata_df.query("analysisId!=@analysis_exclude_list")["ind"].values.tolist()])
    for count,analysis in enumerate(response.json()):
        if analysis_exclude_list!=None:
            if analysis['analysisId'] in analysis_exclude_list:
                continue
        if count%50==0 and debug:
            print(count)
        for file in analysis['files']:
            if file['info'].get('analysis_tools'):
                if file['info']['analysis_tools'][0]=='GATK:CollectOxoGMetrics':
                    metrics.loc[analysis['analysisId'],'sampleId']=analysis['samples'][0]['sampleId']
                    if "star" in file['fileName']:
                        metrics.loc[analysis['analysisId'],"PIPELINE"]="STAR"
                    elif "hisat2" in file['fileName']:
                        metrics.loc[analysis['analysisId'],"PIPELINE"]="HISAT2"
                    else :
                        metrics.loc[analysis['analysisId'],"PIPELINE"]="BWA-MEM"
                    for query in [
                                'oxoQ_score'
                    ]:
                        metrics.loc[analysis['analysisId'],query]=file['info']['metrics'][query]
            else:
                continue

    return(metrics)

def aggregate_picard_mark_duplicates_metrics(response,analysis_exclude_list,debug):
    print("Aggregating metrics from : %s" % ('biobambam2:bammarkduplicates2'))
    metrics=pd.DataFrame()
    #total=len([response.json()[int(ind)] for ind in metadata_df.query("analysisId!=@analysis_exclude_list")["ind"].values.tolist()])
    for count,analysis in enumerate(response.json()):
        if analysis_exclude_list!=None:
            if analysis['analysisId'] in analysis_exclude_list:
                continue
        if count%50==0 and debug:
            print(count)
        for file in analysis['files']:
            if file['info'].get('analysis_tools'):
                if file['info']['analysis_tools'][0]=='biobambam2:bammarkduplicates2':
                    metrics.loc[analysis['analysisId'],'sampleId']=analysis['samples'][0]['sampleId']
                    if "star" in file['fileName']:
                        metrics.loc[analysis['analysisId'],"PIPELINE"]="STAR"
                    elif "hisat2" in file['fileName']:
                        metrics.loc[analysis['analysisId'],"PIPELINE"]="HISAT2"
                    else :
                        metrics.loc[analysis['analysisId'],"PIPELINE"]="BWA-MEM"
                    for query in [
                                'READ_PAIRS_EXAMINED',
                                'READ_PAIR_DUPLICATES',
                                'READ_PAIR_OPTICAL_DUPLICATES',
                                'UNMAPPED_READS',
                                'UNPAIRED_READS_EXAMINED',
                                'UNPAIRED_READ_DUPLICATES'
                    ]:
                        metrics.loc[analysis['analysisId'],query]=sum(z[query] for z in file['info']['metrics']['libraries'])
            else:
                continue

    metrics['TOTAL_READS']=(metrics['READ_PAIRS_EXAMINED']*2)+metrics['UNPAIRED_READS_EXAMINED']
    metrics['DUPLICATION_PCT']=((metrics['READ_PAIR_DUPLICATES']*2)+metrics['UNPAIRED_READ_DUPLICATES'])/metrics['TOTAL_READS']*100
    metrics['MAPPING_PCT']=(metrics['TOTAL_READS']-metrics['UNMAPPED_READS'])/metrics['TOTAL_READS']*100
    return(metrics)
    
def aggreate_picard_collect_rnaseq_metrics(response,metadata_df,analysis_exclude_list,debug):
    print("Aggregating metrics from : %s " % ('Picard:CollectRnaSeqMetrics'))
    metrics=pd.DataFrame()
    for count,analysis in enumerate(response.json()):
        if analysis_exclude_list!=None:
            if analysis['analysisId'] in analysis_exclude_list:
                continue
        if count%50==0 and debug:
            print(count)
        for file in analysis['files']:
            if file['info'].get('analysis_tools'):
                if file['info']['analysis_tools'][0]=='Picard:CollectRnaSeqMetrics':
                    if "star" in file['fileName']:
                        metrics.loc[analysis['analysisId'],"PIPELINE"]="STAR"
                    else:
                        metrics.loc[analysis['analysisId'],"PIPELINE"]="HISAT2"
                    for query in [key for key in file['info']['metrics'].keys() if "pct" in key or "median" in key]:
                        metrics.loc[analysis['analysisId'],query]=file['info']['metrics'][query]
                    metrics.loc[analysis['analysisId'],"sampleId"]=analysis['samples'][0]['sampleId']
            else:
                continue
    return(metrics)
            
            
def song_phone_home(project,rdpc_url,state):
    print("Calling Song API...")
    combined_url="%s/studies/%s/analysis?analysisState=%s" % (rdpc_url,project,state)

    response=requests.get(combined_url)

    if response.status_code!=200:
        sys.exit("Query response failed, return status_code :%s" % response.status_code)

    print("Calling Song API...Complete")  
    return(response)

def generate_rdpc_metadata(response,experiment):
    print("Aggregating Files and IDs...")  
    metadata=pd.DataFrame()
    count=0
    for ind,analysis in enumerate(response.json()):
        if analysis['experiment']['experimental_strategy']==experiment:

            for file in analysis['files']:
                metadata.loc[count,'fileDataType']=file['dataType']
                metadata.loc[count,'objectId']=file['objectId']
                metadata.loc[count,'submitterSampleId']=analysis['samples'][0]['submitterSampleId']
                metadata.loc[count,'submitterSpecimenId']=analysis['samples'][0]['specimen']['submitterSpecimenId']
                metadata.loc[count,'submitterDonorId']=analysis['samples'][0]['donor']['submitterDonorId']

                metadata.loc[count,'donorId']=analysis['samples'][0]['donor']['donorId']
                metadata.loc[count,'specimenId']=analysis['samples'][0]['specimen']['specimenId']
                metadata.loc[count,'sampleId']=analysis['samples'][0]['sampleId']

                metadata.loc[count,'analysisId']=analysis['analysisId']
                metadata.loc[count,'runId']=analysis['workflow']['run_id']
                metadata.loc[count,'ind']=ind
                count+=1
    print("Aggregating Files and IDs...Complete")        
    return(metadata)


if __name__ == "__main__":
    main()
