# gemelli env
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from skbio import DistanceMatrix
from skbio.stats.distance import permanova
from biom import load_table
from biom.table import Table
import seaborn as sns

from gemelli.rpca import rpca
from scipy.spatial import procrustes
from skbio.stats.distance import mantel

from adjustText import adjust_text

table = load_table('../out/rarefied_table.biom')
metadata = pd.read_csv('../qiita_artifacts/metadata.txt', sep='\t')
metadata = metadata.set_index('SampleID')

def run_rpca(table, metadata, target_column):
    ordination, distance = rpca(table)
    samples_df = ordination.samples
    
    target_data = {}

    for sample_id in samples_df.index:
        val = metadata.loc[sample_id, target_column]
        target_data[sample_id] = val

    df = samples_df.copy()

    # Add cohort information to the samples dataframe
    df[target_column] = pd.Series(target_data)

    # Extract variance explained
    variance_explained = ordination.proportion_explained
    
    return ordination, distance, df, variance_explained


df = table.to_dataframe()
df.shape, metadata.shape


metadata['cd_behavior'].value_counts()


metadata['cd_phenotype'] = metadata['cd_behavior'].map({
    'Non-stricturing, non-penetrating (B1)': 'B1', 
    'not applicable': np.nan,
    'Stricturing (B2)': 'B-other',
    'Penetrating (B3)': 'B-other'
})


metadata['cd_phenotype'].value_counts()


target_column = 'cd_phenotype'

ordination, distance, df, variance_explained = run_rpca(table, metadata, target_column)


# Create the RPCA plot with PERMANOVA results annotation
sns.scatterplot(
    data=df,
    x='PC1', 
    y='PC2',
    hue=target_column,
    s=100,
)

plt.title('RPCA Ordination by Crohn\'s Disease Phenotype', fontsize=15)

plt.xlabel(f'PC1 ({variance_explained[0]:.2%} variance explained)')
plt.ylabel(f'PC2 ({variance_explained[1]:.2%} variance explained)')

# Adjust legend
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='cd_phenotype')

# Save plot
plt.savefig('../out/figs/figure6_rpca_ordination.pdf', dpi=300, bbox_inches='tight')





