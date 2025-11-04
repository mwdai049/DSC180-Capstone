# DSC180-Capstone

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
