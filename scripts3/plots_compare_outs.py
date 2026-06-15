import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

"""Requires the tables produced by `compare_outs__gsva_kpca.py`"""

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 13})

palette_box = {'kpca' : "#66c2a5", 'gsva': "skyblue"}
palette_swarm = {'kpca': "green", 'gsva': "cadetblue" }


def add_spatial_params_for_plot(adata):
    adata.uns['spatial'] = {'my_slice_name': {
        'scalefactors': {'fiducial_diameter_fullres': 100,
                         'spot_diameter_fullres': 1,
                         'tissue_hires_scalef': 1,
                         'tissue_lowres_scalef': 1}}}
    return adata


def compute_ordering_tissues(huge: pd.DataFrame, which_values: str):
    """
    huge is a dataframe with columns: pathway_id, <which_values>, tissue
    ranks the tissue names by the median of the values of the
        column <which_values>     across pathway_id s
    """
    assert which_values in ['r', 'abs(r)', 'value'], "Error, column not recognized"
    choo = pd.pivot(huge, columns="pathway_id", index="tissue",
                    values=which_values)
    choo = choo.assign(median_score_per_tissue=choo.median(axis=1).to_numpy())
    choo = choo.assign(tissue_name=choo.index.tolist())
    choo = choo[['median_score_per_tissue', 'tissue_name']]
    choo = choo.sort_values(by=['median_score_per_tissue'], ascending=True)

    return choo.index.tolist()


def count_above_i_cutoff(df, i_thresh:float):
    df['name_w_type'] = df['tissue'].astype(str) + '|' + df[
        'type'].astype(str)
    df = df.loc[df['I'] >= i_thresh,:].copy()
    cool = df.groupby(['name_w_type']).count()
    cool.rename(columns={'pathway_id': 'n_sets'}, inplace=True)
    cool.reset_index(inplace=True)
    cool[['tissue', 'type']] = cool[
        'name_w_type'].str.split("|", expand=True)
    return cool


def mae_between_methods(df, method_a="A", method_b="B"):
    """
    Mean Absolute Error (MAE):
    average absolute difference in n_sets between the two methods across tissues.
    """
    pivot = (
        df.pivot(index="tissue", columns="type", values="n_sets")
          [[method_a, method_b]]
          .dropna()
    )

    return np.mean(np.abs(pivot[method_a] - pivot[method_b]))


def rmse_between_methods(df, method_a="A", method_b="B"):
    """
    Root Mean Squared Error (RMSE):
    typical difference in n_sets between methods, with larger discrepancies weighted more heavily.
    """
    pivot = (
        df.pivot(index="tissue", columns="type", values="n_sets")
          [[method_a, method_b]]
          .dropna()
    )

    return np.sqrt(np.mean((pivot[method_a] - pivot[method_b]) ** 2))


# ******** START **********

scores_avgs_merged = pd.read_csv("scores_avg_merged.csv", index_col=None)

# -- save correlations
# # -- correlations between spot-wise values: corr(gsva, kpca)
correls_merged = pd.read_csv("correls_merged.csv",  index_col=None)

# -- now save the moran I scores of the features for each gsva and kpca
moran_comb = pd.read_csv("moran_merged.csv",  index_col=None)

#      -- compute how many of each category are above an 'I' threshold

i_thresh = 0.1
df_counts = count_above_i_cutoff(moran_comb, i_thresh)

# ## PLOTS
# -- plot correlations
# # -- correlations between spot-wise values: corr(gsva, kpca)
tissues_order_by_correl = compute_ordering_tissues(correls_merged,
                                                       which_values="r")

correls_merged = correls_merged.assign(
    tissue=pd.Categorical(correls_merged['tissue'].tolist(),
                          categories=tissues_order_by_correl))
fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(15, 6))
sns.boxplot(data=correls_merged, x="tissue", y="r", hue="tissue",
            # orient='h',
            palette="PRGn", legend=False,
            flierprops={"marker": "."}, zorder=0)
sns.swarmplot(data=correls_merged, x="tissue", y="r", hue="tissue",
              size=2,
              # color="black", # deprecated !
              palette='dark:black',
              # palette="Set2",
              zorder=1,
              legend=False)
plt.xticks(rotation=90)
plt.tight_layout()
# plt.show()
plt.savefig("figs/correls_comparison.pdf")
plt.close()

# -- just plot the average scores to check their variability
scores_avgs_merged = scores_avgs_merged.assign(
    tissue=pd.Categorical(scores_avgs_merged['tissue'].tolist(),
                          categories=tissues_order_by_correl))
fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(15, 6))
sns.boxplot(data=scores_avgs_merged, x="tissue", y="avg_scores", hue="type",
            palette=palette_box,
            zorder=0,
            flierprops={"marker": "."})

sns.swarmplot(data=scores_avgs_merged, x="tissue", y="avg_scores", hue="type",
              dodge=True, size=2.5,
              palette=palette_swarm,
              zorder=1)
plt.xticks(rotation=90)
plt.tight_layout()
# plt.show()
plt.savefig("figs/avg_scores_comparison.pdf")
plt.close()

# -- now plot the moran I scores of the features for each gsva and kpca
moran_comb = moran_comb.assign(
    tissue=pd.Categorical(moran_comb['tissue'].tolist(),
                          categories=tissues_order_by_correl))
fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(15, 6))
sns.boxplot(data=moran_comb, x="tissue", y="I", hue="type",
            palette=palette_box,
            legend=True,
            flierprops={"marker": "o"}, zorder=0)
sns.stripplot(data=moran_comb, x="tissue", y="I", hue="type",
              size=2.2, dodge=True,
              palette=palette_swarm,
              zorder=1,
              legend=False)
plt.xticks(rotation=90)
plt.tight_layout()
# plt.show()
plt.savefig("figs/moran_i_comparison.pdf")
plt.close()

# -- plot counts of sets over 'I' threshold
df_counts = df_counts.assign(
    tissue=pd.Categorical(df_counts['tissue'].tolist(),
                          categories=tissues_order_by_correl))
fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(15, 6))
sns.barplot(
        data=df_counts,
        x="tissue",
        y="n_sets",
        hue="type",
        dodge=True,
        gap=0.01,
        width=0.6,
        alpha=0.55,
        palette=palette_box
    )

sns.lineplot(
 data=df_counts,
    x="tissue",
    y="n_sets",
    hue="type",
    style="type",
    markers=True,
    palette=palette_swarm,
    alpha=0.9
)
plt.xticks(rotation=90)
plt.tight_layout()
# plt.show()
plt.savefig("figs/n_sets_barplots_lineplot.pdf")
plt.close()

#      - - plot simple boxplot showing global difference
fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(2.6, 4))
sns.boxplot(
    data=df_counts,
    x="type",
    y="n_sets",
    hue="type",
    palette=palette_box,
    zorder=0,
    legend=False,
    #gap=0.2
)
sns.stripplot(data=df_counts, x="type", y="n_sets", hue="type",
              size=4,
              #jitter=True,   # no need of dodge here
              palette=palette_swarm,
              zorder=1,
              legend=False)
plt.xticks(rotation=90)
plt.title("Global difference")
plt.tight_layout()
# plt.show()
plt.savefig("figs/global_diff_n_sets.pdf")
plt.close()


# --  difference between the two curves
mae = mae_between_methods(df_counts, "gsva", "kpca")
rmse = rmse_between_methods(df_counts, "gsva", "kpca")

print(f"Mean absolute error (MAE): {mae:.3f}")
print(f"Root mean squared error (RMSE): {rmse:.3f}")

# ---------------------------------------------------------------------
# # Computed difference:
# ## Mean absolute error (MAE): 1.333
# ## Root mean squared error (RMSE): 2.646
