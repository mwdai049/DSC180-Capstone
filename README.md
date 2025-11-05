# DSC180-Capstone

## About the Project 1 Checkpoint
For our replication, our mentor, Rob Knight, wanted us to focus more on understanding the microbiology field and software tools than replicating the code exactly. To that end, we have all spent weeks 3-5 reading through various clinical trials and statistical analyses to get a better understanding of the field and what we as data scientists can do within it. 

For this reason, the code contained with this repository is focused on QIIME 2, the open-source microbiology data analysis tool. Rather than Python code, we have the step-by-step notebook and artifacts produced by the QIIME 2 pipeline. We are working on replicating this analysis in Python rather than the command line, but for the checkpoint this is what we have accomplished so far. 

## QIIME 2
To upload and process the data, you must install QIIME 2. Instructions can be found here: https://docs.qiime2.org/2024.10/install/#qiime-2-2024-10-distributions

## QIIME 2 Visualizations
In order to view visualizations created in QIIME 2, you will need to upload the output .qza or .qzv file onto the QIIME 2 View website. A link to the website can be found here: https://view.qiime2.org/

## Data
To download the data, visit: https://qiita.ucsd.edu/analysis/description/30622/ (you will need a Qiita account)

## Building phylogenetic and taxonomic trees
greengenes2 is needed for obtaining the phylogenetic and taxanomic trees. The repository can be found here: https://github.com/biocore/q2-greengenes2

Installation requires the following steps: 
```
$ source activate qiime2.2022.8
$ pip install q2-greengenes2
```

To get the trees, you will need to download the following file from the FTP linked here: https://ftp.microbio.me/greengenes_release/current/
- 2024.09.phylogeny.nsv.nwk.qza
