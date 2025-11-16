import pandas as pd
import qiime2 as q2
import seaborn as sns
import matplotlib.pyplot as plt
import re

from statsmodels.stats.power import tt_ind_solve_power
from statistics import mean, stdev
from skbio import DistanceMatrix
from qiime2.plugins import fragment_insertion, diversity, emperor
from qiime2 import Metadata, Artifact
from qiime2.plugins.feature_table.methods import rarefy
from qiime2.plugins.phylogeny.methods import midpoint_root
from skbio import TreeNode
from q2_types.tree import NewickFormat
from skbio.diversity.alpha import faith_pd
import plotly.express as px
from qiime2 import Artifact, Metadata
import skbio
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Load the feature table and insertion tree. Both are QIIME2 artifacts.
ft = Artifact.load('qiita_artifacts/feature_table.qza')

# Import the relabelled Newick tree as a rooted QIIME 2 phylogeny artifact and save it as insertion_tree.qza
nwk = NewickFormat('qiita_artifacts/insertion_tree.relabelled.tre', mode='r')
tree_qza = Artifact.import_data('Phylogeny[Rooted]', nwk)
tree_qza.save('artifacts/insertion_tree.qza')

insertion_tree = q2.Artifact.load('artifacts/insertion_tree.qza')

# Rarefaction (a kind of random subsampling) on the feature table of counts. 
# We sample down to 10,000 reads per sample to make them comparable, and storing that new standardized table as rarefied
rarefied = rarefy(table=ft, sampling_depth=10000).rarefied_table

# Compute Faith's Phylogenetic Diversity (alpha diversity)
faith_res = diversity.actions.alpha_phylogenetic(
    table=rarefied,
    phylogeny=insertion_tree,
    metric='faith_pd'
)

faith_pd = faith_res.alpha_diversity.view(pd.Series)
faith_res.alpha_diversity.save('artifacts/alpha_diversity.qza')

# Compute unweighted UniFrac distance matrix (beta diversity)
unifrac_res = diversity.actions.beta_phylogenetic(
    table=rarefied,
    phylogeny=insertion_tree,
    metric='unweighted_unifrac'
)
unifrac_dm = unifrac_res.distance_matrix.view(DistanceMatrix)
unifrac_df = pd.DataFrame(unifrac_dm.data, index=unifrac_dm.ids, columns=unifrac_dm.ids)
unifrac_res.distance_matrix.save('artifacts/distance_matrix.qza')

# Load sample metadata and merge with alpha diversity results. Subset data based on Crohn's disease behavior.
metadata = pd.read_csv('qiita_artifacts/metadata.txt', sep='\\t', dtype=str, index_col=0, engine='python')

# Ensure that the index of the metadata matches the sample IDs in the diversity results
metadata['deblur alpha diversity'] = faith_pd
metadata.dropna(subset=['deblur alpha diversity'], inplace=True)

# Subset metadata for Crohn's disease behavior groups
b1 = metadata[metadata['cd_behavior'] == 'Non-stricturing, non-penetrating (B1)']
bother = metadata[(metadata['cd_behavior'] != 'Non-stricturing, non-penetrating (B1)') & 
                  (metadata['cd_behavior'] != 'not applicable')]

b1_dtx = unifrac_dm.filter(ids=b1.index).to_series().values
bother_dtx = unifrac_dm.filter(ids=bother.index).to_series().values

# Power analysis for alpha and beta diversity comparisons
# alpha
sd1 = b1['deblur alpha diversity'].std()
sd2 = bother['deblur alpha diversity'].std()
sd12 = mean([sd1, sd2])

alpha = 0.05
data_alpha = []
for observations in range(10, 155, 5):
    for difference in [2, 3, 4]:
        effect_size = difference/sd12
        x = tt_ind_solve_power(effect_size=effect_size, nobs1=observations, 
                               alpha=alpha, ratio=1.0,
                               alternative='two-sided')
        data_alpha.append({
            'alpha': alpha, 'Power (1-β)': x, 
            'Total sample size (N)': observations * 2, 
            'Difference': '%d (effect size %0.2f)' % (
                difference, effect_size)})
data_alpha = pd.DataFrame(data_alpha) 

# beta
u2u1 = abs(mean(b1_dtx) - mean(bother_dtx))
effect_size = u2u1/mean([stdev(b1_dtx), stdev(bother_dtx)])

data_beta = []
for observations in range(10, 155, 5):
    for alpha in [.001, .01, .05, .1]:
        x = tt_ind_solve_power(effect_size=effect_size, nobs1=observations, 
                               alpha=alpha, ratio=1.0,
                               alternative='two-sided')
        data_beta.append({
            'Significance level (α)': alpha, 'Power (1-β)': x, 
            'Total sample size (N)': observations * 2, 
            'Difference': '%d (effect size %0.2f)' % (
                difference, effect_size)})
data_beta = pd.DataFrame(data_beta) 

print ('Alpha diversity:')
print('B1 Mean Alpha Diversity: ', b1['deblur alpha diversity'].mean())
print('B Other Mean Alpha Diversity: ', bother['deblur alpha diversity'].mean())
print('B1 STD: ', b1['deblur alpha diversity'].std())
print('B Other STD: ', bother['deblur alpha diversity'].std())

print ('Beta diversity:')
print('B1 Mean Beta Diversity: ', mean(b1_dtx))
print('B Other Mean Beta Diversity: ', mean(bother_dtx))
print('B1 STD: ', stdev(b1_dtx))
print('B Other STD: ', stdev(bother_dtx))

# Plotting the alpha diversity power analysis results
fig, ax1 = plt.subplots(figsize=(15, 9))

sns.set(style="whitegrid")
sns.set_context("paper", font_scale=2, 
                rc={'lines.linewidth': 2, 'lines.markersize': 12})

f = sns.lineplot(x='Total sample size (N)', y='Power (1-β)',
                 markers=True, dashes=False, style='Difference',
                 ax=ax1, data=data_alpha)
ax1.xaxis.set_major_locator(plt.MultipleLocator(20))

plt.axhline(0.8, 0, data_alpha['Total sample size (N)'].max())
plt.show()
fig.savefig('figs/figure1.pdf')

# Plotting beta diversity power analysis and distributions
fig = plt.figure(figsize=(20, 8))
grid = plt.GridSpec(ncols=3, nrows=1, hspace=0.2, wspace=0.2)

ax1 = fig.add_subplot(grid[0, :2])
ax2 = fig.add_subplot(grid[0, 2:])

# LEFT: beta power curves
sns.lineplot(x='Total sample size (N)', y='Power (1-β)',
             style='Significance level (α)',
             markers=True, dashes=False,
             ax=ax1, data=data_beta)
ax1.axhline(0.8, 0, data_beta['Total sample size (N)'].max())
ax1.xaxis.set_major_locator(plt.MultipleLocator(20))

# RIGHT: distributions using histplot (+ KDE)
sns.histplot(b1_dtx,
             label="B1 within distances",
             color="red",
             kde=True,
             stat="density",
             element="step",
             ax=ax2)

ax2.axvline(mean(b1_dtx), 0, 6, color="red")

sns.histplot(bother_dtx,
             label="B2-3 within distances",
             color="skyblue",
             kde=True,
             stat="density",
             element="step",
             ax=ax2)

ax2.axvline(mean(bother_dtx), 0, 6, color="skyblue")

ax2.xaxis.set_major_locator(plt.MultipleLocator(.1))
ax2.legend()
plt.show()
fig.savefig('figs/figure2.pdf')

# Perform PCoA on the unweighted UniFrac distance matrix and visualize with Emperor
pcoa = diversity.methods.pcoa(distance_matrix=unifrac_dm)
viz = emperor.visualizers.plot(pcoa=pcoa.pcoa, metadata=metadata)
viz.visualization.save('unweighted-unifrac-pcoa.qza')

pcoa_art = Artifact.load('unweighted-unifrac-pcoa.qza')

ordination: skbio.OrdinationResults = pcoa_art.view(skbio.OrdinationResults)
coords = ordination.samples.copy()
coords.index.name = 'SampleID'

meta_df = metadata.to_dataframe()
df = coords.join(meta_df, how='inner')

x_lab = f"PC1 ({ordination.proportion_explained['PC1'] * 100:.1f}% var)"
y_lab = f"PC2 ({ordination.proportion_explained['PC2'] * 100:.1f}% var)"

fig = px.scatter(
    df,
    x='PC1',
    y='PC2',
    color='cd_behavior',       
    hover_name=df.index,         
    hover_data=df.columns,       
    title="Unweighted UniFrac PCoA"
)

fig.update_layout(
    xaxis_title=x_lab,
    yaxis_title=y_lab
)

fig.savefig('figs/figure3.html')
