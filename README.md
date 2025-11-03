# DSC180-Capstone

## Uploading and processing the data
To upload and process the data, you must install QIIME 2. Instructions can be found here: https://docs.qiime2.org/2024.10/install/#qiime-2-2024-10-distributions

## QIIME 2 Visualizations
In order to view visualizations created in QIIME 2, you will need to upload the output .qza or .qzv file onto the QIIME 2 View website. A link to the website can be found here: https://view.qiime2.org/

## Building phylogenetic and taxonomic trees
greengenes2 is needed for building the phylogenetic and taxanomic trees. The repository can be found here: https://github.com/biocore/q2-greengenes2

Installation requires the following steps: 
```
$ source activate qiime2.2022.8
$ pip install q2-greengenes2
```

To build the trees, you will need to download the following files from the FTP linked here: https://ftp.microbio.me/greengenes_release/current/
- 2024.09.phylogeny.id.nwk.qza
