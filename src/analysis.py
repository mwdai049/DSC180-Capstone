import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re
import numpy as np

import qiime2 as q2
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

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.inspection import permutation_importance

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Load the feature table and insertion tree. Both are QIIME2 artifacts.
ft = Artifact.load('qiita_artifacts/feature-table.qza')

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
metadata_df = pd.read_csv('qiita_artifacts/metadata.txt', sep='\t', dtype=str, index_col=0, engine='python')

# Ensure that the index of the metadata matches the sample IDs in the diversity results
metadata_df['deblur alpha diversity'] = faith_pd
metadata_df.dropna(subset=['deblur alpha diversity'], inplace=True)

# Subset metadata for Crohn's disease behavior groups
b1 = metadata_df[metadata_df['cd_behavior'] == 'Non-stricturing, non-penetrating (B1)']
bother = metadata_df[(metadata_df['cd_behavior'] != 'Non-stricturing, non-penetrating (B1)') & 
                  (metadata_df['cd_behavior'] != 'not applicable')]

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

alpha_b1_mean = b1['deblur alpha diversity'].mean()
alpha_bother_mean = bother['deblur alpha diversity'].mean()
alpha_b1_std = b1['deblur alpha diversity'].std()
alpha_bother_std = bother['deblur alpha diversity'].std()

beta_b1_mean = mean(b1_dtx)
beta_bother_mean = mean(bother_dtx)
beta_b1_std = stdev(b1_dtx)
beta_bother_std = stdev(bother_dtx)

summary_rows = [
    {"metric": "alpha", "group": "B1",      "mean": alpha_b1_mean,     "std": alpha_b1_std},
    {"metric": "alpha", "group": "B-other", "mean": alpha_bother_mean, "std": alpha_bother_std},
    {"metric": "beta",  "group": "B1",      "mean": beta_b1_mean,      "std": beta_b1_std},
    {"metric": "beta",  "group": "B-other", "mean": beta_bother_mean,  "std": beta_bother_std},
]

summary_df = pd.DataFrame(summary_rows)
summary_df_rounded = summary_df.copy()
summary_df_rounded["mean"] = summary_df_rounded["mean"].round(3)
summary_df_rounded["std"] = summary_df_rounded["std"].round(3)
summary_df_rounded.to_csv("figs/diversity_summary_stats.csv", index=False)

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
pcoa_res = diversity.methods.pcoa(distance_matrix=unifrac_res.distance_matrix)
metadata = Metadata.load('qiita_artifacts/metadata.txt')
viz = emperor.visualizers.plot(
    pcoa=pcoa_res.pcoa,
    metadata=metadata
)

ordination: skbio.OrdinationResults = pcoa_res.pcoa.view(skbio.OrdinationResults)
coords = ordination.samples.copy()
coords.index.name = 'SampleID'

df = coords.join(metadata_df, how='inner')

df = df[df['cd_behavior'] != 'not applicable'].copy()

df['behavior_group'] = df['cd_behavior'].apply(
    lambda x: 'B1' if x == 'Non-stricturing, non-penetrating (B1)' else 'B-other'
)
df['behavior_group'] = pd.Categorical(df['behavior_group'],
                                      categories=['B1', 'B-other'],
                                      ordered=True)

pc1_col = coords.columns[0]
pc2_col = coords.columns[1]

pc1_var = ordination.proportion_explained.iloc[0] * 100
pc2_var = ordination.proportion_explained.iloc[1] * 100

x_lab = f"{pc1_col} ({pc1_var:.1f}% var)"
y_lab = f"{pc2_col} ({pc2_var:.1f}% var)"

fig = px.scatter(
    df,
    x=pc1_col,
    y=pc2_col,
    color='behavior_group',   
    hover_name=df.index,
    hover_data=df.columns,
    title="Unweighted UniFrac PCoA"
)

fig.update_layout(
    xaxis_title=x_lab,
    yaxis_title=y_lab
)

fig.write_html('figs/figure3.html')

#Random Forest Classification
meta_rf = metadata_df.copy()
meta_rf = meta_rf[meta_rf['cd_behavior'] != 'not applicable'].copy()
meta_rf['behavior_group'] = np.where(
    meta_rf['cd_behavior'] == 'Non-stricturing, non-penetrating (B1)',
    'B1',
    'B-other'
)

rarefied_df = rarefied.view(pd.DataFrame)
common_ids = rarefied_df.index.intersection(meta_rf.index)

X_full = rarefied_df.loc[common_ids]
y_full = meta_rf.loc[common_ids, 'behavior_group']

X_train, X_test, y_train, y_test = train_test_split(
    X_full,
    y_full,
    test_size=0.3,
    stratify=y_full,
    random_state=42,
)

def rf_learning_curve(X_train,y_train,X_test,y_test,sample_sizes,n_reps=10,random_state=42):
    rng = np.random.default_rng(random_state)

    Xtr = X_train.values
    ytr = y_train.values
    Xte = X_test.values
    yte = y_test.values

    classes = np.unique(ytr)
    if len(classes) != 2:
        raise ValueError(f"Expected 2 classes, got {classes}")

    class_to_idx = {c: np.where(ytr == c)[0] for c in classes}
    min_class_size = min(len(idx) for idx in class_to_idx.values())
    max_balanced_N = 2 * min_class_size

    sample_sizes = [N for N in sample_sizes if N <= max_balanced_N]

    results = []

    yte_bin = (yte == classes[1]).astype(int)

    for N in sample_sizes:
        n_per_class = N // 2
        aucs = []

        for rep in range(n_reps):
            chosen_idxs = []

            for c in classes:
                idxs = class_to_idx[c]
                chosen = rng.choice(idxs, size=n_per_class, replace=False)
                chosen_idxs.append(chosen)

            chosen_idxs = np.concatenate(chosen_idxs)

            X_sub = Xtr[chosen_idxs]
            y_sub = ytr[chosen_idxs]

            clf = RandomForestClassifier(
                n_estimators=500,
                max_depth=None,
                class_weight="balanced",
                n_jobs=-1,
                random_state=random_state + rep,
            )
            clf.fit(X_sub, y_sub)

            y_prob = clf.predict_proba(Xte)[:, 1]
            auc = roc_auc_score(yte_bin, y_prob)
            aucs.append(auc)

        results.append(
            {
                "N_total": N,
                "mean_auc": np.mean(aucs),
                "std_auc": np.std(aucs),
                "n_reps": len(aucs),
            }
        )

    return pd.DataFrame(results)

sample_sizes = list(range(40, X_train.shape[0] * 2, 20))

rf_results = rf_learning_curve(
    X_train,
    y_train,
    X_test,
    y_test,
    sample_sizes,
    n_reps=10,
)

print(rf_results)

plt.figure(figsize=(8, 6))
sns.lineplot(x="N_total", y="mean_auc", data=rf_results, marker="o")
plt.fill_between(
    rf_results["N_total"],
    rf_results["mean_auc"] - rf_results["std_auc"],
    rf_results["mean_auc"] + rf_results["std_auc"],
    alpha=0.2,
)
plt.axhline(0.5, color="gray", linestyle="--", label="Chance (AUC=0.5)")
plt.xlabel("Total sample size (B1 + B-other)")
plt.ylabel("Test-set ROC AUC")
plt.title("Random Forest learning curve: B1 vs B-other")
plt.legend()
plt.tight_layout()
plt.savefig("figs/figure4.pdf")
plt.show()

#Which features are most important in classification?
rf_final = RandomForestClassifier(
    n_estimators=500,
    max_depth=None,
    class_weight='balanced',
    n_jobs=-1,
    random_state=42,
)
rf_final.fit(X_train, y_train)

importances = rf_final.feature_importances_

feat_names = X_train.columns 
imp_df = pd.DataFrame({
    'feature': feat_names,
    'gini_importance': importances
})

imp_df = imp_df.sort_values('gini_importance', ascending=False)
imp_df.to_csv('figs/rf_feature_importances.csv', index=False)

#Permutation importance
y_prob = rf_final.predict_proba(X_test)[:, 1]
baseline_auc = roc_auc_score((y_test == rf_final.classes_[1]).astype(int), y_prob)

perm = permutation_importance(
    rf_final,
    X_test,
    y_test,
    n_repeats=20,
    random_state=42,
    n_jobs=-1,
    scoring='roc_auc',
)

perm_importances = perm.importances_mean  
perm_df = pd.DataFrame({
    'feature': X_test.columns,
    'perm_importance_auc_drop': perm_importances
})

perm_df = perm_df.sort_values('perm_importance_auc_drop', ascending=False)
perm_df.to_csv('figs/rf_feature_importances_permutation.csv', index=False)
