import os
import anndata as ad
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.spatial import distance
import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from smpath.helpers import generate_coordinates_arr
import smpath.processing.measures_methods as spm
from scipy.stats import zscore

def normalize_scores(scores_df, mode="zscore"):
    if mode == "zscore":
        df_z = pd.DataFrame(
            zscore(scores_df, axis=1, nan_policy="omit"),
            index=scores_df.index,
            columns=scores_df.columns
        )
        return df_z
    else:
        print("not implemented, will return None")
        return None


def open_scores_objects(out_parent_dir, dir_adata, SAMPLE):
    matrix_gsva_kegg = pd.read_csv(os.path.join(
        out_parent_dir, SAMPLE, f"{SAMPLE}-cmps-mx-KEGG-KEGG_gsva.tsv"
    ), sep='\t', index_col=0)
    adata_score_kegg = ad.read_h5ad(os.path.join(dir_adata, SAMPLE,
                                                 f"{SAMPLE}-scores-KEGG.h5ad"))
    kpca_kegg = adata_score_kegg.to_df()

    print(
        f"initial sets numbers: kpca -> {kpca_kegg.shape[1]}, gsva -> {matrix_gsva_kegg.shape[1]}")

    print("  ", kpca_kegg.columns)
    print("  ", matrix_gsva_kegg.columns)

    ok_ids = list(set(kpca_kegg.columns.tolist()
                      ).intersection(set(matrix_gsva_kegg.columns.tolist())))

    ok_spots = list(set(kpca_kegg.index.tolist()
                        ).intersection(set(matrix_gsva_kegg.index.tolist())))

    kpca_kegg = kpca_kegg.loc[ok_spots, ok_ids]
    matrix_gsva_kegg = matrix_gsva_kegg.loc[ok_spots, ok_ids]

    return kpca_kegg, matrix_gsva_kegg


def return_pearsons_all_feats(kpca_kegg,
                              matrix_gsva_kegg) -> pd.DataFrame:
    correls_res = list()

    assert np.all(
        kpca_kegg.columns == matrix_gsva_kegg.columns), "Error, columns do not match, aborting!"
    assert np.all(
        kpca_kegg.index == matrix_gsva_kegg.index), "Error, rows do not match, aborting!"

    for u in kpca_kegg.columns:
        res = pearsonr(matrix_gsva_kegg[u].to_numpy(),
                       kpca_kegg[u].to_numpy())
        correls_res.append(res.statistic)

    return pd.DataFrame({"pathway_id": kpca_kegg.columns, "r": correls_res})


def return_flattened_scores(kpca_kegg: pd.DataFrame,
                            matrix_gsva_kegg: pd.DataFrame,
                            return_spotwise=False):
    assert np.all(
        kpca_kegg.columns == matrix_gsva_kegg.columns), "Error, columns do not match, aborting!"
    assert np.all(
        kpca_kegg.index == matrix_gsva_kegg.index), "Error, rows do not match, aborting!"

    axis_int = 0  # by default, across pathways (faster computation)
    if return_spotwise:
        axis_int = 1

    flattened_kpca = kpca_kegg.to_numpy().var(axis=axis_int)  # .flatten()
    flattened_gsva = matrix_gsva_kegg.to_numpy().var(
        axis=axis_int)  # .flatten()

    types_list = list(np.repeat("kpca", len(flattened_kpca)
                                )) + list(
        np.repeat("gsva", len(flattened_gsva)))

    df = pd.DataFrame(
        {"avg_scores": list(flattened_kpca) + list(flattened_gsva),
         "type": types_list})

    return df


def compute_ordering_tissues(huge: pd.DataFrame, which_values: str):
    """
    huge is a dataframe with columns: pathway_id, r, tissue
    ranks the tissue names by the median of the r values across pathway_id s
    """
    assert which_values in ['r', 'abs(r)'], "Error, must be 'r' or 'abs(r)'"
    choo = pd.pivot(huge, columns="pathway_id", index="tissue",
                    values=which_values)
    choo = choo.assign(median_score_per_tissue=choo.median(axis=1).to_numpy())
    choo = choo.assign(tissue_name=choo.index.tolist())
    choo = choo[['median_score_per_tissue', 'tissue_name']]
    choo = choo.sort_values(by=['median_score_per_tissue'], ascending=True)

    return choo.index.tolist()


def do_binary_mask(v_array, centile_value=0.7):
    assert centile_value <= 1, "Error, centile_value must be between 0 and 1"
    try:
        value = np.percentile(v_array, q=centile_value * 100)
        print(v_array.max(), v_array.min())
        print("centile; ", value)
        mask = np.zeros(shape=v_array.shape, dtype=np.uint32)
        mask[mask >= value] = 1
    except Exception as e:
        print(e)

    return mask


def dice_coefficient(mask1, mask2):
    mask1 = np.asarray(mask1, dtype=bool)
    mask2 = np.asarray(mask2, dtype=bool)

    if mask1.shape != mask2.shape:
        raise ValueError("Masks must have the same shape")

    intersection = np.logical_and(mask1, mask2).sum()

    return 2.0 * intersection / (mask1.sum() + mask2.sum())


def plot_this_pathway(pathway_id, pdseries, modality: str,
                      ax):
    coords_strings_l = pdseries.index.tolist()
    if ('x' in coords_strings_l[0]) and ('y' in coords_strings_l[0]):
        coords_arr = generate_coordinates_arr(coords_strings_l)

    df = pd.DataFrame(data=coords_arr, columns=['x', 'y'])
    df = df.assign(values=pdseries.to_numpy())
    sns.scatterplot(data=df, x='x', y='y', hue='values',
                    palette='Spectral_r', legend=True,
                    s=5,
                    ax=ax)
    ax.set_title(f'{pathway_id} : {modality}')

    return ax


def wrap_comp_moran(df:pd.DataFrame, coords_df:pd.DataFrame,
                    set_type:str, tissue:str):
    adata = ad.AnnData(
        X=df.values,
        obs=coords_df,
        var=pd.DataFrame(index=df.columns)
    )

    adata.obsm['spatial'] = np.array(coords_df)
    adata = spm.compute_moran_i(adata, n_perms=2)

    print(adata.uns.keys())

    moran_df = adata.uns['moranI'].copy()
    moran_df['pathway_id'] = moran_df.index.tolist()
    moran_df['type'] = set_type
    moran_df['tissue'] = tissue

    return moran_df[['pathway_id', 'I', 'type', 'tissue']]



if __name__ == '__main__':

    print(os.listdir("../../"), "%%%%%%")
    out_parent_dir = '../../apollo-gsva/compar_out_2026'  # TODO: modify in server: compar_out_2026
    dir_adata = "../../apollo-data"

    tisues_df = pd.read_csv("../the_data_list.tsv", sep="\t",
                              index_col=None, header=0)
    tissues_list = tisues_df['dataset_name'].tolist()   #  [-4:]
    print(tissues_list)

    correls_dfs = list()
    scores_each_dfs = list()
    morans_dfs = list()

    for SAMPLE in tissues_list:   #  ["LP-brain-mouse-7"]
        try:
            kpca_kegg, matrix_gsva_kegg = open_scores_objects(out_parent_dir,
                                                              dir_adata,
                                                              SAMPLE)
            coords_df = pd.read_csv(os.path.join(out_parent_dir,
                                                 SAMPLE, f'{SAMPLE}-coords.tsv'),
                                    sep='\t', index_col=0, header=0)
        except Exception as e:
            print(e, "failed opening")

        kpca_I_df = wrap_comp_moran(kpca_kegg, coords_df, "kpca",
                                    tissue=SAMPLE)
        gsva_I_df = wrap_comp_moran(matrix_gsva_kegg, coords_df, "gsva", tissue=SAMPLE)

        morans_dfs.append(pd.concat([kpca_I_df, gsva_I_df], ignore_index=True))

        try:
            kpca_kegg = normalize_scores(kpca_kegg)
            matrix_gsva_kegg = normalize_scores(matrix_gsva_kegg)
            if len(kpca_kegg.columns) > 1:
                tmp_df = return_pearsons_all_feats(kpca_kegg, matrix_gsva_kegg)
                tmp_df = tmp_df.assign(**{'abs(r)': tmp_df['r'].abs()})
                tmp_df['tissue'] = SAMPLE
                correls_dfs.append(tmp_df)

                tmp_flat = return_flattened_scores(kpca_kegg, matrix_gsva_kegg)
                tmp_flat['tissue'] = SAMPLE
                scores_each_dfs.append(tmp_flat)
            # end if
        except:
            print(SAMPLE, "=== ? what happened with this tissue ? ")


    megaflat = pd.concat(scores_each_dfs, axis=0, ignore_index=True)

    # -- plot correlations
    # # -- correlations between spot-wise values: corr(gsva, kpca)

    correls_merged = pd.concat(correls_dfs, axis=0, ignore_index=True)

    tissues_order_by_correl = compute_ordering_tissues(correls_merged,
                                                       which_values="abs(r)")

    correls_merged = correls_merged.assign(
        tissue=pd.Categorical(correls_merged['tissue'].tolist(),
                              categories=tissues_order_by_correl))

    sns.boxplot(data=correls_merged, x="tissue", y="abs(r)", hue="tissue",
                # orient='h',
                palette="PuBu_r", legend=False,
                flierprops={"marker": "."}, zorder=0)
    sns.swarmplot(data=correls_merged, x="tissue", y="abs(r)", hue="tissue",
                  size=2,
                  # color="black", # deprecated !
                  palette='dark:black',
                  # palette="Set2",
                  zorder=1,
                  legend=False)
    plt.xticks(rotation=90)
    plt.show()
    plt.close()

    # -- just plot the average scores to check their variability
    megaflat = megaflat.assign(
        tissue=pd.Categorical(megaflat['tissue'].tolist(),
                              categories=tissues_order_by_correl))

    palette_box = {'kpca' : "#66c2a5", 'gsva': "skyblue"}
    palette_swarm = {'kpca': "green", 'gsva': "cadetblue" }

    sns.boxplot(data=megaflat, x="tissue", y="avg_scores", hue="type",
                palette=palette_box,
                zorder=0,
                flierprops={"marker": "."} )


    sns.swarmplot(data=megaflat, x="tissue", y="avg_scores", hue="type",
                  dodge=True, size=2.5,
                  palette=palette_swarm,
                  zorder=1)
    plt.xticks(rotation=90)
    plt.show()
    plt.close()

    # -- now plot the moran I scores of the features for each gsva and kpca
    moran_comb = pd.concat(morans_dfs, ignore_index=True)

    moran_comb = moran_comb.assign(
        tissue=pd.Categorical(moran_comb['tissue'].tolist(),
                              categories=tissues_order_by_correl))

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
    plt.show()
    plt.close()


