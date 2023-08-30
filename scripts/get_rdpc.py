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
from datetime import datetime

def main():
    """

    """
    parser = argparse.ArgumentParser(description='Retrieve stats from SONG API and generate plots')
    parser.add_argument('-p', '--project', dest="project", help="projects to query", required=True,type=str,nargs="+")
    parser.add_argument('-u', '--url', dest="rdpc_url", help="RDPC GraphQL URL", required=True,type=str)
    parser.add_argument("-t","--token", dest="token", help="EGO token", required=True,type=str)
    parser.add_argument('-o', '--output_directory', dest="out_dir", help="SONG RDPC URL", default=os.getcwd(),type=str)
    parser.add_argument('-r', '--repo', dest="repo", help="Workflow Github Repo", required=True,type=str,nargs="+")
    parser.add_argument('-x', '--exclude_runs', dest="excluded_runs", help="Runs to exclude",nargs="+",type=str)
    parser.add_argument('-z', '--plot', dest="plot", help="make pretty plots", default=True,type=bool)
    parser.add_argument('-d', '--debug', dest="debug", help="debug", default=False,type=bool)
    parser.add_argument('-s', '--state', dest="state", help="analysis state to query : True False",default=["COMPLETE"],choices=["EXECUTOR_ERROR","SYSTEM_ERROR","COMPLETE"],nargs="+")

    cli_input= parser.parse_args()
    excluded_runs=cli_input.excluded_runs
    debug=cli_input.debug
    rdpc_stats={}
    for repo in cli_input.repo:
        rdpc_stats[repo]={}
        for state in cli_input.state:
            response=rdpc_phone_home(cli_input.rdpc_url,cli_input.token,state,repo)
            run_aggregate_df,task_aggregate_df=generate_rdpc_aggregates(response,debug)

            for project in cli_input.project:
                write_dir="%s/%s_%s_%s" % (cli_input.out_dir,repo.replace(".git","").split("/")[-1],state,project)
                if not os.path.exists(write_dir):
                    os.mkdir(write_dir)

                tsv_dir="%s/%s" % (write_dir,"tsv")
                if not os.path.exists(tsv_dir):
                    os.mkdir(tsv_dir)

                rdpc_stats[repo][project]={}
                rdpc_stats[repo][project][state]={}
                if debug:
                    print(len(run_aggregate_df))
                    print(len(task_aggregate_df))
                
                rdpc_stats[repo][project][state]['runs']=run_aggregate_df.query("study_id==@project and runId!=@excluded_runs")
                rdpc_stats[repo][project][state]['tasks']=task_aggregate_df.query("study_id==@project and runId!=@excluded_runs")
                if debug:
                    print(len(rdpc_stats[repo][project][state]['runs']))
                    print(len(rdpc_stats[repo][project][state]['tasks']))


                rdpc_stats[repo][project][state]['runs'].to_csv("%s/%s_%s_%s_%s.tsv" % (tsv_dir,repo.replace(".git","").split("/")[-1],project,str(state),"runs"),sep="\t")
                rdpc_stats[repo][project][state]['tasks'].to_csv("%s/%s_%s_%s_%s.tsv" % (tsv_dir,repo.replace(".git","").split("/")[-1],project,str(state),"tasks"),sep="\t")

                plots={}
                for ind,item in enumerate(["total_realtime_hrs","max_mem_gb"]):
                    title="%s %s %s" % (repo.replace(".git","").split("/")[-1],project,item)
                    plots["fig.%s.%s.%s" % (1,ind+1,title.replace(" ","_"))]=generate_plot(
                        run_aggregate_df,
                        500,
                        500,
                        [project],
                        [item],
                        title
                    )

                save_pkl_plots(write_dir,plots,cli_input.plot)



    
 
def save_pkl_plots(out_dir,generated_plots,plot):
    svg_dir="%s/%s" % (out_dir,"svg")
    pkl_dir="%s/%s" % (out_dir,"pkl")

    if not os.path.exists(svg_dir):
        os.mkdir(svg_dir)
    if not os.path.exists(pkl_dir):
        os.mkdir(pkl_dir)

    print("Saving plots...")
    for generated_plot in generated_plots.keys():
        file = open("%s/%s.pkl" % (pkl_dir,generated_plot),"wb")
        pickle.dump(generated_plots[generated_plot],file)
        file.close()
    print("Saving plots...Complete")

    if plot:
        print("Saving plots SVGs...")
        for generated_plot in generated_plots.keys():
            generated_plots[generated_plot].write_image("%s/%s.svg" % (svg_dir,generated_plot))
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
                    x=metrics.query("study_id==@col").sort_values(row)['sample_id'].values.tolist(),
                    y=metrics.query("study_id==@col").sort_values(row)[row].values.tolist(),
                    mode='markers+lines',
                    showlegend=False),
                row_ind+1,
                col_ind+1
            )
            
            for val in [25,5,75]:
                fig.append_trace(
                    go.Scatter(
                        x=[
                            metrics.query("study_id==@col").sort_values(row)['sample_id'].values.tolist()[0],
                            metrics.query("study_id==@col").sort_values(row)['sample_id'].values.tolist()[-1]
                        ],
                        y=[np.percentile(metrics.query("study_id==@col").sort_values(row)[row].values.tolist(),val,axis=0)]*2,
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
            
def rdpc_phone_home(rdpc_url,token,state,repo):
    print("Calling RDPC API...")
    
    headers = {
    'Accept-Encoding': 'gzip, deflate, br',
    'Content-Type': 'application/json',
    'Accept': 'application/json',
    'Connection': 'keep-alive',
    'DNT': '1',
    'Origin': '%s',
    'authorization': 'Bearer %s' % (token),
    }

    variables={
        "RunsFilter":{"repository":repo,"state":state},"analysisPage":{"from":0,"size":10000}
    }

    query=\
    """
    query($analysisPage: Page,$RunsFilter: RunsFilter) {
    runs(
        filter: $RunsFilter
        page: $analysisPage
        sorts: { fieldName: startTime, order: asc }
    ) {
        content {
            runId
            startTime
            completeTime
            state
            duration
            repository
            inputAnalyses{
                analysisId
                studyId
                donors{
                    donorId
                    specimens{
                        specimenId
                        samples{
                            sampleId
                            }
                    }
                }
            }
        tasks {
            process
            cpus
            duration
            peakRss
            peakVmem
            readBytes
            writeBytes
            container
            realtime
            vmem
            runId
            memory
            name
            sessionId
            startTime
            completeTime
            realtime
            state
            }
        }
    }}
    """
    response = requests.post(rdpc_url, json={'query': query,"variables":variables},headers=headers)
    if response.status_code!=200:
        sys.exit("Query response failed, return status_code :%s" % response.status_code)

    print("Calling RDPC API...Complete")

    if len(response.json()['data']['runs']['content'])==0:
        sys.exit("Query returned zero responses. Check Repo URL :%s" % repo)        

    return(response.json()['data']['runs']['content'])

def generate_rdpc_aggregates(response,debug):
    print("Aggregating Files and IDs...")  
    run_aggregate_df=pd.DataFrame()
    task_aggregate_df=pd.DataFrame()
    agg_count=0
    task_count=0

    for run in response:
        #print(run['runId'])

        if len(run['inputAnalyses'])==0:
          if debug:
            print("Skipping %s" % (run['runId']))
          continue
        run_aggregate_df.loc[agg_count,"runId"]=run['runId']
        run_aggregate_df.loc[agg_count,"repository"]=run['repository']
        run_aggregate_df.loc[agg_count,"study_id"]=run['inputAnalyses'][0]['studyId']
        run_aggregate_df.loc[agg_count,"analysis_id"]=run['inputAnalyses'][0]['studyId']
        run_aggregate_df.loc[agg_count,"specimen_id"]=run['inputAnalyses'][0]['donors'][0]['specimens'][0]['specimenId']
        run_aggregate_df.loc[agg_count,"sample_id"]=run['inputAnalyses'][0]['donors'][0]['specimens'][0]['samples'][0]['sampleId']
        run_aggregate_df.loc[agg_count,"donor_id"]=run['inputAnalyses'][0]['donors'][0]['donorId']
        run_aggregate_df.loc[agg_count,'startTime']=datetime.fromtimestamp(int(run['startTime'])/1000).strftime('%Y-%m-%dT%H:%M:%S') if run['startTime'] else None
        run_aggregate_df.loc[agg_count,'completeTime']=datetime.fromtimestamp(int(run['completeTime'])/1000).strftime('%Y-%m-%dT%H:%M:%S') if run['completeTime'] else None
        run_aggregate_df.loc[agg_count,"total_realtime_hrs"]=sum([task['realtime'] if task['realtime'] else 0 for task in run['tasks']])/3.6e+6
        run_aggregate_df.loc[agg_count,"total_duration_hrs"]=sum([task['duration'] if task['duration'] else 0  for task in run['tasks']])/3.6e+6
        run_aggregate_df.loc[agg_count,"task_misisng_info_count"]=len([task for task in run['tasks'] if task['duration']==None or task['realtime']==None or task['memory']==None])
        run_aggregate_df.loc[agg_count,"duration_hrs"]=run['duration']/3.6e+6 if run['duration'] else 0
        run_aggregate_df.loc[agg_count,"max_mem_gb"]=np.max([task['memory']/8/1073741824 if task['memory'] else 0 for task in run['tasks']])
        run_aggregate_df.loc[agg_count,"max_cpus"]=np.max([task['cpus'] if task['cpus'] else 0 for task in run['tasks']])
        run_aggregate_df.loc[agg_count,"state"]=run['state']
        if run['startTime'] and run['completeTime']:
            run_aggregate_df.loc[agg_count,"realtime"]=(int(run['completeTime'])/3.6e+6)-(int(run['startTime'])/3.6e+6)
        agg_count+=1
        if len(run['tasks'])==0:
            continue
        for task in run['tasks']:
            task_aggregate_df.loc[task_count,"sessionId"]=task['sessionId']
            task_aggregate_df.loc[task_count,"study_id"]=run['inputAnalyses'][0]['studyId']
            task_aggregate_df.loc[task_count,"cpus"]=task['cpus']
            task_aggregate_df.loc[task_count,'startTime']=datetime.fromtimestamp(int(task['startTime'])/1000).strftime('%Y-%m-%dT%H:%M:%S') if task['startTime'] else None
            task_aggregate_df.loc[task_count,'completeTime']=datetime.fromtimestamp(int(task['completeTime'])/1000).strftime('%Y-%m-%dT%H:%M:%S') if task['completeTime'] else None
            task_aggregate_df.loc[task_count,'realtime']=task['realtime']/3.6e+6 if task['realtime'] else 0
            task_aggregate_df.loc[task_count,'container']=task['container']
            task_aggregate_df.loc[task_count,'memory_gb']=task['memory']/8/1073741824 if task['memory'] else 0
            task_aggregate_df.loc[task_count,"runId"]=run['runId']
            task_aggregate_df.loc[task_count,"repository"]=run['repository']
            task_aggregate_df.loc[task_count,'readBytes']=task['readBytes']
            task_aggregate_df.loc[task_count,'writeBytes']=task['writeBytes']
            task_aggregate_df.loc[task_count,'peakVmem']=task['peakVmem']
            task_aggregate_df.loc[task_count,'peakRss']=task['peakRss']
            task_aggregate_df.loc[task_count,'process']=task['process']
            task_aggregate_df.loc[task_count,'state']=task['state']
            task_count+=1

    print("Aggregating Files and IDs...Complete")        
    return(run_aggregate_df,task_aggregate_df)


if __name__ == "__main__":
    main()
