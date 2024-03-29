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
import warnings
warnings.filterwarnings('ignore')

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
    parser.add_argument('-y', '--plot_level', dest="plot_level",default=['sample','donor'],help="sample or donor level plots",nargs="+",type=str,choices=['sample','donor'])
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
                    os.makedirs(write_dir)

                tsv_dir="%s/%s" % (write_dir,"tsv")
                if not os.path.exists(tsv_dir):
                    os.makedirs(tsv_dir)

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
                        cli_input.excluded_analyses,
                        cli_input.debug
                    )
                    metrics['biobambam2:bammarkduplicates2']=aggregate_picard_mark_duplicates_metrics(
                        response,
                        cli_input.excluded_analyses,
                        cli_input.debug
                    )                    
                    for plot_level in cli_input.plot_level:
                        if len(metrics['biobambam2:bammarkduplicates2'])>0:
                            for ind,item in enumerate(['TOTAL_READS','DUPLICATION_PCT','MAPPING_PCT']):
                                title="%s %s %s %s" % (project,experiment,plot_level+"Lvl",item)
                                plots["fig.%s.%s.%s" % (1,ind+1,title.replace(" ","_"))]=generate_plot(
                                    metadata[project][experiment][state],
                                    metrics['biobambam2:bammarkduplicates2'],
                                    1000,
                                    600,
                                    ["STAR","HISAT2"],
                                    [item],
                                    title,
                                    plot_level
                                )
                        if len(metrics['Picard:CollectRnaSeqMetrics'])>0:
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
                                title="%s %s %s %s" % (project,experiment,plot_level+"Lvl",item)
                                plots["fig.%s.%s.%s" % (2,ind+1,title.replace(" ","_"))]=generate_plot(
                                    metadata[project][experiment][state],
                                    metrics['Picard:CollectRnaSeqMetrics'],
                                    1000,
                                    600,
                                    ["STAR","HISAT2"],
                                    [item],
                                    title,
                                    plot_level
                                )
                    for key,name in zip(
                        ['Picard:CollectRnaSeqMetrics','biobambam2:bammarkduplicates2'],
                        ["rnaMetrics","libraryMetrics"]
                    ):
                        if len(metrics[key])>0:
                            metrics[key].to_csv("%s/%s_%s_%s_%s.tsv" % (tsv_dir,state,project,experiment,name),sep="\t",index=True)

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

                    for plot_level in cli_input.plot_level:
                        if len(metrics['biobambam2:bammarkduplicates2'])>0:
                            for ind,item in enumerate(['TOTAL_READS','DUPLICATION_PCT','MAPPING_PCT']):
                                title="%s %s %s %s" % (project,experiment,plot_level+"Lvl",item)
                                plots["fig.%s.%s.%s" % (1,ind+1,title.replace(" ","_"))]=generate_plot(
                                    metadata[project][experiment][state],
                                    metrics['biobambam2:bammarkduplicates2'],
                                    1000,
                                    600,
                                    ["BWA-MEM"],
                                    [item],
                                    title,
                                    plot_level
                                )
                        if len(metrics['GATK:CollectOxoGMetrics'])>0:
                            for ind,item in enumerate(['oxoQ_score']):
                                title="%s %s %s %s" % (project,experiment,plot_level+"Lvl",item)
                                plots["fig.%s.%s.%s" % (1,ind+1,title.replace(" ","_"))]=generate_plot(
                                    metadata[project][experiment][state],
                                    metrics['GATK:CollectOxoGMetrics'],
                                    1000,
                                    600,
                                    ["BWA-MEM"],
                                    [item],
                                    title,
                                    plot_level
                                )

                        if len(metrics['Samtools:stats'])>0:
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
                                title="%s %s %s %s" % (project,experiment,plot_level+"Lvl",item)
                                plots["fig.%s.%s.%s" % (1,ind+1,title.replace(" ","_"))]=generate_plot(
                                    metadata[project][experiment][state],
                                    metrics['Samtools:stats'],
                                    1000,
                                    600,
                                    ["BWA-MEM"],
                                    [item],
                                    title,
                                    plot_level
                                )

                        if len(metrics['Picard:CollectQualityYieldMetrics'])>0 and plot_level=='sample':
                            for ind,item in enumerate(['total_reads','read_length','pf_reads']):
                                title="%s %s %s %s" % (project,experiment,plot_level+"Lvl",item)
                                plots["fig.%s.%s.%s" % (1,ind+1,title.replace(" ","_"))]=generate_plot(
                                    metadata[project][experiment][state],
                                    metrics['Picard:CollectQualityYieldMetrics'],
                                    1000,
                                    600,
                                    ["BWA-MEM"],
                                    [item],
                                    title,
                                    plot_level
                                )

                        if len(metrics['Sanger:verifyBamHomChk'])>0 and plot_level=='sample':
                            for ind,item in enumerate(['avg_depth','contamination','reads_used','snps_used']):
                                title="%s %s %s %s" % (project,experiment,plot_level+"Lvl",item)
                                plots["fig.%s.%s.%s" % (1,ind+1,title.replace(" ","_"))]=generate_plot(
                                    metadata[project][experiment][state],
                                    metrics['Sanger:verifyBamHomChk'],
                                    1000,
                                    600,
                                    ["BWA-MEM"],
                                    [item],
                                    title,
                                    plot_level
                                )

                        if len(metrics['Sanger:compareBamGenotypes'])>0 and plot_level=='sample':
                            for ind,item in enumerate(['total_loci_genotype','frac_match_gender','frac_informative_genotype','frac_matched_genotype']):
                                title="%s %s %s %s" % (project,experiment,plot_level+"Lvl",item)
                                plots["fig.%s.%s.%s" % (1,ind+1,title.replace(" ","_"))]=generate_plot(
                                    metadata[project][experiment][state],
                                    metrics['Sanger:compareBamGenotypes'],
                                    1000,
                                    600,
                                    ["BWA-MEM"],
                                    [item],
                                    title,
                                    plot_level
                                )
                    for key,name in zip(
                        [
                            'biobambam2:bammarkduplicates2',
                            'GATK:CollectOxoGMetrics',
                            'Samtools:stats',
                            'Picard:CollectQualityYieldMetrics',
                            'Sanger:verifyBamHomChk',
                            'Sanger:compareBamGenotypes'
                        ],
                        ["markDupMetrics","oxoMetrics","samtoolsMetrics","readGroupMetrics","verifyBamMetrics","bamGenotypesMetrics"]
                    ):
                        if len(metrics[key])>0:
                            metrics[key].to_csv("%s/%s_%s_%s_%s.tsv" % (tsv_dir,state,project,experiment,name),sep="\t",index=True)

                    save_pkl_plots(write_dir,plots,cli_input.plot)   

def save_pkl_plots(out_dir,gen_plots,plot):
    print("Saving plots...")
    svg_dir="%s/%s" % (out_dir,"svg")
    pkl_dir="%s/%s" % (out_dir,"pkl")

    if not os.path.exists(svg_dir):
        os.makedirs(svg_dir)
    if not os.path.exists(pkl_dir):
        os.makedirs(pkl_dir)
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

def generate_plot(metadata,metrics,x_dim,y_dim,cols,rows,title,plot_level):
    if plot_level=='sample':
        return(generate_sample_plot(metrics,x_dim,y_dim,cols,rows,title))
    else:
        return(generate_donor_plot(metadata,metrics,x_dim,y_dim,cols,rows,title))

def generate_sample_plot(metrics,x_dim,y_dim,cols,rows,title):
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
        xaxis=dict(title="sampleId"),
        yaxis=dict(title=""),
        showlegend=True,
        titlefont=dict(size=20)
    )
    return(fig)

def generate_donor_plot(metadata,metrics,x_dim,y_dim,cols,rows,title):
    print("Generating plot for %s" % (title))
    fig=plotly.subplots.make_subplots(
        cols=len(cols),
        rows=len(rows),
        subplot_titles=cols
    )

    for row_ind,row in enumerate(rows):
        for col_ind,col in enumerate(cols):
            samples=metrics.query("PIPELINE==@col")['sampleId'].values.tolist()
            ###Return tumour samples
            subset_metadata=metadata.query("sampleId==@samples and tumourNormalDesignation=='Tumour'").loc[:,["sampleId","matchedNormalSampleId","donorId"]].drop_duplicates().set_index("sampleId")
            ###Reorder tumour samples according to metric value
            subset_metadata=subset_metadata.loc[metrics.set_index("sampleId").loc[subset_metadata.index.values.tolist(),row].sort_values().index.values.tolist(),:]
            ###Return metric value for tumour and normal samples
            tumour_values=metrics.set_index("sampleId").loc[subset_metadata.index.values.tolist(),row].values.tolist()
            normal_values=[metrics.set_index("sampleId").loc[ind,row] if ind in metrics.set_index("sampleId").index else None for ind in subset_metadata['matchedNormalSampleId'].values.tolist()]
            
            ###Plot donor vs metrics from normal
            fig.append_trace(
                go.Scatter(
                    x=subset_metadata['donorId'].values.tolist(),
                    y=normal_values,
                    mode='markers+lines',
                    name='Normal',
                    line=dict(dash="dash",color="green"),
                    showlegend=True),
                row_ind+1,
                col_ind+1
            )
            ###Plot donor vs metrics from tumour
            fig.append_trace(
                go.Scatter(
                    x=subset_metadata['donorId'].values.tolist(),
                    y=tumour_values,
                    mode='markers+lines',
                    name='Tumour',
                    line=dict(dash="dash",color="red"),
                    showlegend=True),
                row_ind+1,
                col_ind+1
            )
            for val in [25,50,75]:
                ###Print percentiles for tumours
                fig.append_trace(
                    go.Scatter(
                        x=[
                            subset_metadata['donorId'].values.tolist()[0],
                            subset_metadata['donorId'].values.tolist()[-1]
                        ],
                        y=[np.percentile([z for z in tumour_values if z!=None],val,axis=0)]*2,
                        mode='lines',
                        line=dict(dash="dash",color="black"),
                        opacity=0.6,
                        name="Tumour Percentiles : 25,50,75",
                        showlegend=True if val==25 else False),
                    row_ind+1,
                    col_ind+1
                )
                ###Print percentiles for normals
                fig.append_trace(
                    go.Scatter(
                        x=[
                            subset_metadata['donorId'].values.tolist()[0],
                            subset_metadata['donorId'].values.tolist()[-1]
                        ],
                        y=[np.percentile([z for z in normal_values if z!=None],val,axis=0)]*2,
                        mode='lines',
                        line=dict(dash="dot",color="black"),
                        opacity=0.4,
                        name="Normal Percentiles : 25,50,75",
                        showlegend=True if val==25 else False),
                    row_ind+1,
                    col_ind+1
                )

    fig['layout'].update(
        width=x_dim,
        height=y_dim,
        title=title,
        xaxis=dict(title="donorId"),
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

                    for query in [
                        "total_reads",
                        "read_length",
                        'pf_reads'
                    ]:
                        if query in file['info']['metrics']:
                            metrics.loc[uniqueId,query]=file['info']['metrics'][query]
                        else:
                            print("Analysis %s is missing field %s" % (analysis['analysisId'],query))
                            metrics.loc[uniqueId,query]=None
            else:
                continue

    if debug:
        print(metrics)

    if len(metrics)==0:
        print("No Data Found!")

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
                        if query in file['info']['metrics']:
                            metrics.loc[analysis['analysisId'],query]=file['info']['metrics'][query]
                        else:
                            print("Analysis %s is missing field %s" % (analysis['analysisId'],query))
                            metrics.loc[analysis['analysisId'],query]=None
            else:
                continue

    if debug:
        print(metrics)
        
    if len(metrics)==0:
        print("No Data Found!")

    return(metrics)

def aggregate_sanger_compareBamGenotypes_metrics(response,analysis_exclude_list,debug):
    print("Aggregating metrics from : %s" % ('Sanger:compareBamGenotypes'))
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

                            for query in [
                                'frac_match_gender',
                                'gender'
                            ]:
                                if query in file['info']['metrics']['tumours'][0]['gender']:
                                    metrics.loc[analysis['analysisId'],query]=file['info']['metrics']['tumours'][0]['gender'][query]
                                else:
                                    print("Analysis %s is missing field %s" % (analysis['analysisId'],query))
                                    metrics.loc[analysis['analysisId'],query]=None

                            for query in [
                                'compared_against',
                                'total_loci_gender',
                                'total_loci_genotype'
                            ]:
                                if query in file['info']['metrics']:
                                    metrics.loc[analysis['analysisId'],query]=file['info']['metrics'][query]
                                else:
                                    print("Analysis %s is missing field %s" % (analysis['analysisId'],query))
                                    metrics.loc[analysis['analysisId'],query]=None

                            for query in [
                                'frac_informative_genotype',
                                'frac_matched_genotype'

                            ]:
                                if query in file['info']['metrics']['tumours'][0]['genotype']:
                                    metrics.loc[analysis['analysisId'],query]=file['info']['metrics']['tumours'][0]['genotype'][query]
                                else:
                                    print("Analysis %s is missing field %s" % (analysis['analysisId'],query))
                                    metrics.loc[analysis['analysisId'],query]=None
            else:
                continue

    if debug:
        print(metrics)
        
    if len(metrics)==0:
        print("No Data Found!")

    return(metrics)

def aggregate_sanger_verifyBamHomChk_metrics(response,analysis_exclude_list,debug):
    print("Aggregating metrics from : %s" % ('Sanger:verifyBamHomChk'))
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
    if debug:
        print(metrics)

    if len(metrics)==0:
        print("No Data Found!")

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
                        if query in file['info']['metrics']:
                            metrics.loc[analysis['analysisId'],query]=file['info']['metrics'][query]
                        else:
                            print("Analysis %s is missing field %s" % (analysis['analysisId'],query))
                            metrics.loc[analysis['analysisId'],query]=None

            else:
                continue

    if debug:
        print(metrics)

    if len(metrics)==0:
        print("No Data Found!")

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
    

    if debug:
        print(metrics)
        
    if len(metrics)==0:
        print("No Data Found!")

    return(metrics)
    
def aggreate_picard_collect_rnaseq_metrics(response,analysis_exclude_list,debug):
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

    if debug:
        print(metrics)
        
    if len(metrics)==0:
        print("No Data Found!")

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

                metadata.loc[count,'tumourNormalDesignation']=analysis['samples'][0]['specimen']['tumourNormalDesignation']
                metadata.loc[count,'matchedNormalSubmitterSampleId']=analysis['samples'][0]['matchedNormalSubmitterSampleId']

                metadata.loc[count,'donorId']=analysis['samples'][0]['donor']['donorId']
                metadata.loc[count,'specimenId']=analysis['samples'][0]['specimen']['specimenId']
                metadata.loc[count,'sampleId']=analysis['samples'][0]['sampleId']

                metadata.loc[count,'analysisId']=analysis['analysisId']
                metadata.loc[count,'runId']=analysis['workflow']['run_id']
                metadata.loc[count,'ind']=ind
                count+=1

    warnings=[]
    normalSubmitterToArgo={z[0]:z[1] for z in metadata.loc[:,["submitterSampleId","sampleId"]].drop_duplicates().values.tolist()}
    for ind in metadata.query("tumourNormalDesignation=='Tumour'").index.values.tolist():
        submitterId=metadata.loc[ind,"matchedNormalSubmitterSampleId"]
        if normalSubmitterToArgo.get(submitterId):
            metadata.loc[ind,"matchedNormalSampleId"]=normalSubmitterToArgo.get(submitterId)
        else:
            warnings.append("No ID associated with %s based on SONG records" % submitterId)
            
    for warning in list(set(warnings)):
        print(warning)

    print("Aggregating Files and IDs...Complete")        
    return(metadata)


if __name__ == "__main__":
    main()
