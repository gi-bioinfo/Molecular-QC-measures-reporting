# Molecular-QC-measures-reporting
## Example run:
```
python get_analysis.py -p APGI-AU -u https://song.rdpc-qa.cancercollaboratory.org -e RNA-Seq
python get_rdpc.py -p APGI-AU -t ${ego_token} -u https://api.rdpc-qa.cancercollaboratory.org/graphql -r https://github.com/icgc-argo-workflows/rna-seq-alignment -o .
```

## Requirements:
[numpy](https://anaconda.org/anaconda/numpy)<Br>
[pandas](https://anaconda.org/anaconda/pandas)<Br>
[plotly](https://anaconda.org/conda-forge/plotly)<Br>
[pickle](https://anaconda.org/conda-forge/pypickle/files)<Br>
[kaleido](https://anaconda.org/conda-forge/python-kaleido)