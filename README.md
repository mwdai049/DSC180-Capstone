# DSC180-Capstone

## Microbiome Diversity & Power Analysis

## Setup
### QIIME 2
To upload and process the data, you must install QIIME 2. The yml files are in the project directory.

### Conda
#### Linux/Windows WSL

```
conda env create \
  --name qiime2-amplicon-2025.10 \
  --file qiime2-environment-windows.yml
```

#### macOS (Apple Silicon)

```
CONDA_SUBDIR=osx-64 conda env create \
  --name qiime2-amplicon-2025.10 \
  --file qiime2-environment-macos.yml
conda activate qiime2-amplicon-2025.10
conda config --env --set subdir osx-64
```

#### macOS (Intel)

```
conda env create \
  --name qiime2-amplicon-2025.10 \
  --file qiime2-environment-macos.yml
```

To test the installation:

```
conda deactivate
conda activate qiime2-amplicon-2025.10
qiime info
```

### Docker

## Data
To download the data, visit: https://qiita.ucsd.edu/analysis/description/25761/ (you will need a Qiita account)

Place the following under the `qiita_artifacts/` directory:
1. insertion_tree.relabelled.tre
    - Newick tree to be imported as a rooted phylogeny.
2. feature_table.qza
    - Feature table (`FeatureTable[Frequency]`) QIIME 2 artifact.
3. metadata.txt
    - The metadata

In Qiita, the first two can be found in the first triangle node of the processing tree labelled dflt_name (BIOM). The mapping file can be found at the top, next to the title of the analysis.

### Directory Structure

```text
src/
    analysis.py
    notebook.ipynb       
qiita_artifacts/
    feature-table.qza
    insertion_tree.relabelled.tre
    metadata.txt
artifacts/
figs/  
qiime2-environment-macos.yml
qiime2-environment-windows.yml
README.md
```

`src/analysis.py` is the code script, `src/notebook.ipynb` contains the code in a notebook, `artifacts/` will contain created artifacts, and `figs/` will contain the visualizations of results of the analysis.

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
6. PCoA Visualization
    - Performs Principal Coordinates Analysis (PCoA) on the unweighted UniFrac distance matrix
    - Generates an interactive plot 
7. Random Forest Classification
    - Trains a random forest classifier on a range of total sample sizes
    - Evaluates performance via test-set ROC AUC and aggregates mean and standard deviation over multiple repetitions
    - Plots mean AUC vs. total sample size
8. Feature Importance Analysis
    - Trains a final RandomForestClassifier on the full training set.
    - Extracts Gini feature importances, sorted from most important to least important
    - Computes permutation importance on the test set using ROC AUC as the scoring metric
        - Estimates the drop in AUC when each feature is permuted.
9. Output
    - figs/figure1.pdf – alpha-diversity power curves
    - figs/figure2.pdf – beta-diversity power curves and distance distributions
    - figs/diversity_summary_stats.csv - summary statistics of alpha and beta diversities
    - figs/figure3.html - interactive PCoA plot
    - figs/figure4.pdf - mean AUC over different sample sizes
    - figs/rf_feature_importances.csv - extracted Gini feature importances
    - figs/rf_feature_importances_permutation.csv - permutation-based importance results

## Running the Code
Activate the QIIME 2 Python environment, then from the project directory:

```
conda activate qiime2-amplicon-2025.7 #if not already activated
python src/analysis.py
```

Running this code will run the analysis and generate the resulting figures for visualization. After it finishes, check the figs/ directory for the generated PDF figures and the console output for summary statistics.
