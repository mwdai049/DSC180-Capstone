# DSC180-Capstone

## Microbiome Diversity & Power Analysis

## Credit
Monica Dai and Katelyn Zhao worked on the code replication part of the study, as seen in Github. Camille Sicat helped with the Github repository formatting and project logistics. Camille Sicat, Nathan Wang, and Sophie Wang wrote the report. 

## Setup
### QIIME 2
To upload and process the data, you must install QIIME 2. We used the amplicon distribution, and the yml files are in the project directory.

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

### Gemelli

To run the RPCA analysis, you will need a separate environment compatible with gemelli. The yml file is provided, and you may install by running the following command in terminal:

```
conda env create -f gemelli-env.yml
```

## Data
To download the data, visit: https://qiita.ucsd.edu/analysis/description/25761/ (you will need a Qiita account)

Place the following under the `qiita_artifacts/` directory:
1. insertion_tree.relabelled.tre
    - Newick tree to be imported as a rooted phylogeny.
2. feature_table.qza
    - Feature table (`FeatureTable[Frequency]`) QIIME 2 artifact.
3. metadata.txt
    - The metadata/mapping file

In Qiita, the first two can be found in the first triangle node of the processing tree labelled dflt_name (BIOM). The mapping file can be found at the top, next to the title of the analysis.

### Directory Structure

```text
src/
    analysis.py
    analysis.ipynb     
    rpca_ordination.py
    rpca_ordination.ipynb  
qiita_artifacts/
    feature-table.qza
    insertion_tree.relabelled.tre
    metadata.txt
out/
    figs/  
qiime2-environment-macos.yml
qiime2-environment-windows.yml
gemelli-env.yml
README.md
```

The `analysis.py` and `analysis.ipynb` files contain the same code, and make up most of the data processing and analysis that we conduct. We just provided different ways to execute. `rpca_ordination.py` and `rpca_ordination.ipynb` exist for the same reason, and are separate from previous analysis due to different environment requirements. `out/` will contain any intermediate outputs, and `out/figs/` will contain the visualizations of the analysis results.

## What the Script Does
1. Load the QIIME 2 Artifacts
    - Loads the feature table (qiita_artifacts/feature-table.qza)
    - Imports insertion_tree.relabelled.tre as a rooted phylogenetic tree and saves it as out/insertion_tree.qza
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
    - SHAP value calculations and summary visualizations
9. Output (paths relative to `./out/`)
    - insertion_tree.qza - phylogenetic tree in compatible format with QIIME 2
    - rarefied_table.qza - feature table after the rarefaction subsampling step
    - rarefied_table.biom - rarefied feature table in BIOM format for RPCA compability
    - alpha_diversity.qza - alpha diversity resulting artifact
    - distance_matrix.qza - beta diversity resulting artifact
    - diversity_summary_stats.csv - summary statistics of alpha and beta diversities
    - random_forest_results.csv - auc values over different sample sizes
    - rf_feature_importances.csv - extracted Gini feature importances

    - figs/figure1_power_vs_sample_size_alpha.pdf – alpha-diversity power curves
    - figs/figure2_power_vs_sample_size_beta.pdf – beta-diversity power curves and distance distributions
    - figs/figure3_unifrac_pcoa.html - interactive PCoA plot
    - figs/figure4_random_forest_learning_curve.pdf - mean AUC over different sample sizes
    - figs/figure5_shap_summary.pdf - SHAP summary for feature importance for class B1
    - figs/figure6_rpca_ordination.pdf - RPCA ordination plot

## Running analysis.py
Activate the QIIME 2 Python environment, then from the project directory:

```
conda activate qiime2-amplicon-2025.7 #if not already activated
python src/analysis.py
```

Running this code will run the analysis and generate the resulting figures for visualization. After it finishes, check the out/figs/ directory for the generated figures.

## Running analysis.ipynb
After selecting the QIIME 2 environment, you may need to edit the env's `kernel.json` file if the Setup cell does not resolve import issues related to `rpy2`.

In a terminal, activate the environment and run the following commands to retrieve relevant paths:

```
R_HOME=$(R RHOME)
echo "R_HOME=$R_HOME"
```

```
ENV_BIN=$(dirname "$(which R)")
echo "ENV_BIN=$ENV_BIN"
```

```
LD_LIB="$R_HOME/lib"
echo "LD_LIBRARY_PATH=$LD_LIB"
```

Then find the path to the environment's `kernel.json` using the terminal command `jupyter kernelspec list`, and add the following:

```
{
 ...
 ...,
 "env": {
  "R_HOME": "[REPLACE WITH R_HOME]",
  "PATH": "[REPLACE WITH ENV_BIN]",
  "LD_LIBRARY_PATH": "[REPLACE WITH LD_LIB]" 
 }
}
```

Then restart the notebook and you should no longer run into `rpy2` import issues.

## Running RPCA Code

To run RPCA analysis, you will need to first convert `rarefied_table.qza` into BIOMV210 format. To do so, activate your QIIME 2 environment and run the following command:

```
qiime tools export --input-path out/rarefied_table.qza --output-path out/
mv out/feature-table.biom out/rarefied_table.biom
```

After this preparation step, you will use the `gemelli-env` environment for the rest of the analysis. You may run the code using either the script or the notebook.

