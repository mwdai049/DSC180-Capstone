# DSC180-Capstone

## Microbiome Diversity & Power Analysis
For our replication, our mentor, Rob Knight, wanted us to focus more on understanding the microbiology field and software tools than replicating the code exactly. To that end, we have all spent weeks 3-5 reading through various clinical trials and statistical analyses to get a better understanding of the field and what we as data scientists can do within it. 

For this reason, the code contained with this repository is focused on QIIME 2, the open-source microbiology data analysis tool. Rather than Python code, we have the step-by-step notebook and artifacts produced by the QIIME 2 pipeline. We are working on replicating this analysis in Python rather than the command line, but for the checkpoint this is what we have accomplished so far. 

## QIIME 2
To upload and process the data, you must install QIIME 2. Instructions can be found here: https://docs.qiime2.org/2024.10/install/#qiime-2-2024-10-distributions

## QIIME 2 Visualizations
In order to view visualizations created in QIIME 2, you will need to upload the output .qza or .qzv file onto the QIIME 2 View website. A link to the website can be found here: https://view.qiime2.org/

## Data
To download the data, visit: https://qiita.ucsd.edu/analysis/description/25761/ (you will need a Qiita account)

Place the following under the `qiita_artifacts/` directory:
1. insertion_tree.relabelled.tre
    - Newick tree to be imported as a rooted phylogeny.
2. feature_table.qza
    - Feature table (`FeatureTable[Frequency]`) QIIME 2 artifact.
3. mapping file
    - The metadata

The first two can be found in the first triangle node of the processing tree labelled dflt_name (BIOM). The mapping file can be found at the top, next to the title of the analysis.

### Directory Structure

```text
analysis.py          
qiita_artifacts/
    78873_feature-table.qza
    insertion_tree.relabelled.tre
    metadata.txt
artifacts/
figs/  
```

`analysis.py` is the code script, `artifacts/` will contain created artifacts, and `figs/` will contain the visualizations of results of the analysis.

## What the Script Does
1. Load the QIIME2 Artifacts
    - Loads the feature table (feature-table.qza)
    - Imports insertion_tree.relabelled.tre as a rooted phylogenetic tree and saves it as artifacts/insertion_tree.qza
2. Rarefaction
    - Rarefies (a method for subsampling) the feature table to a sequencing depth of 10,000 reads per sample
3. Diversity Metrics
    - Computes Faith’s Phylogenetic Diversity (alpha diversity)
    - Computes unweighted UniFrac distance matrix (beta diversity)
4. Metadata merge & subsetting
    - Merges Faith’s PD data with metadata.txt by Subject ID
    - Subsets samples into:
        - b1 – cd_behavior == "Non-stricturing, non-penetrating (B1)"
        - bother – all other Crohn’s behaviors (excluding not applicable)
5. Power analysis
    - Alpha diversity: two-sample t-test power for a range of effect sizes and total sample sizes
    - Beta diversity: uses within-group UniFrac distances to derive an effect size and computes power over varying sample sizes and significance levels
6. Output
    - Prints mean and standard deviation for alpha and beta diversity in both groups
    - Saves:
        - figs/figure1.pdf – alpha-diversity power curves
        - figs/figure2.pdf – beta-diversity power curves and distance distributions

## Running the Code
Activate the QIIME 2 Python environment, then from the project directory:

```
python analysis.py
```

Running this code will run the analysis and generate the resulting figures for visualization. After it finishes, check the figs/ directory for the generated PDF figures and the console output for summary statistics.
