import os
import anndata as ad
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.stats import pearsonr
import seaborn as sns
import numpy as np
import pandas as pd
import scanpy as sc
from smpath.helpers import generate_coordinates_arr
import smpath.processing.measures_methods as spm
from scipy.stats import zscore
import squidpy as sq

np.random.seed(42)

"""
Important:
Only the SpacePath's k-PCA implementation allows "mz aware" set enrichment.

The GSVA tool has no functionality to make the scoring enrichment
 "mz aware" (at least 2 different IDs from at least 2 distinct mz);
so we take the pathways that were both enriched by kPCA and GSVA.
"""



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


def open_scores_objects_full(out_parent_dir, dir_adata, SAMPLE):
    """
    opens both gsva and kpca matrices (dataframes)
    the columns' names are the set identifiers, the rows' names the spots.
    returns both separated matrices, sets ids are the column names.
    """
    matrix_gsva_kegg = pd.read_csv(os.path.join(
        out_parent_dir, SAMPLE, f"{SAMPLE}-cmps-mx-KEGG-KEGG_gsva.tsv"
    ), sep='\t', index_col=0)
    adata_score_kegg = ad.read_h5ad(os.path.join(dir_adata, SAMPLE,
                                                 f"{SAMPLE}-scores-KEGG.h5ad"))

    coords_df = pd.DataFrame(adata_score_kegg.obsm['spatial'],
                             index=adata_score_kegg.obs.index.tolist(), columns=['x', 'y'])
    kpca_kegg = adata_score_kegg.to_df()
  

    ok_spots = coords_df.index.tolist()
    kpca_kegg = kpca_kegg.loc[ok_spots, :]
    matrix_gsva_kegg = matrix_gsva_kegg.loc[ok_spots, kpca_kegg.columns.tolist()]


    return kpca_kegg, matrix_gsva_kegg, coords_df


def subset_intersection_scores_objects(kpca_kegg:pd.DataFrame, matrix_gsva_kegg:pd.DataFrame):
    """
    opens both matrices (dataframes)
    drops NaN columns
    returns both separately, with the intersection of the columns' names
    """
    kpca_kegg = kpca_kegg.dropna(axis=1, how="all")
    matrix_gsva_kegg = matrix_gsva_kegg.dropna(axis=1, how="all")

    ok_ids = list(set(kpca_kegg.columns.tolist()
                      ).intersection(set(matrix_gsva_kegg.columns.tolist())))

    kpca_kegg = kpca_kegg.loc[:, ok_ids]
    matrix_gsva_kegg = matrix_gsva_kegg.loc[:, ok_ids]

    return kpca_kegg, matrix_gsva_kegg


def return_pearsons_all_feats(kpca_kegg,
                              matrix_gsva_kegg) -> pd.DataFrame:
    kpca_kegg, matrix_gsva_kegg = subset_intersection_scores_objects(
        kpca_kegg, matrix_gsva_kegg
    )
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


def return_variances_of_scores(kpca_kegg: pd.DataFrame,
                               matrix_gsva_kegg: pd.DataFrame,
                               return_spotwise=False):
    assert np.all(
        kpca_kegg.columns == matrix_gsva_kegg.columns), "Error, columns do not match, aborting!"
    assert np.all(
        kpca_kegg.index == matrix_gsva_kegg.index), "Error, rows do not match, aborting!"

    axis_int = 0  # by default, across pathways (faster computation)
    if return_spotwise:
        axis_int = 1   # across spots (very slow when too many spots)

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

    # shift values to be all positive so avoids errors
    min_val_all = df.min(skipna=True).min(skipna=True)

    df = df + abs(min_val_all) + 1e-20
    coords_df['obs'] = coords_df.index.tolist()
    coords_df.index = coords_df.index.astype(str).tolist()

    adata = ad.AnnData(
        X=csr_matrix(df.to_numpy()),
        var=pd.DataFrame({'pathway_id': np.array(df.columns.tolist(), dtype=str)},
                         index=df.columns.tolist(), dtype=str),
        obs=coords_df,
    )
    adata.obsm['spatial'] = coords_df[['x', 'y']].to_numpy()

    adata = spm.compute_moran_i(adata, n_perms=2) # fast, no need pvalues here

    moran_df = adata.uns['moranI'].copy()
    moran_df['pathway_id'] = moran_df.index.tolist()
    moran_df['type'] = set_type
    moran_df['tissue'] = tissue

    return moran_df[['pathway_id', 'I', 'type', 'tissue']]



if __name__ == '__main__':

    out_parent_dir = '../../compar_out_2026'  # TODO: modify in server: compar_out_2026
    dir_adata = "../../data" # TODO: modify in server: data

    tisues_df = pd.read_csv("../the_data_list.tsv", sep="\t",
                              index_col=None, header=0)
    tissues_list = tisues_df['dataset_name'].tolist()      # [23:25]
    print(tissues_list)

    correls_dfs = list()
    scores_each_dfs = list()
    morans_dfs = list()

    for SAMPLE in tissues_list:
        try:
            kpca_kegg, matrix_gsva_kegg, coords_df = open_scores_objects_full(
                out_parent_dir,  dir_adata,  SAMPLE)

        except Exception as e:
            print(e, "failed opening")

        try:
            kpca_I_df = wrap_comp_moran(kpca_kegg, coords_df, "kpca",
                                        tissue=SAMPLE)
            gsva_I_df = wrap_comp_moran(matrix_gsva_kegg, coords_df, "gsva", tissue=SAMPLE)

            morans_dfs.append(pd.concat([kpca_I_df, gsva_I_df], ignore_index=True))
        except Exception as e:
            print(e, "failed moran of features for this tissue")

        try:
            kpca_kegg = normalize_scores(kpca_kegg)
            matrix_gsva_kegg = normalize_scores(matrix_gsva_kegg)
            if len(kpca_kegg.columns) > 1:
                tmp_df = return_pearsons_all_feats(kpca_kegg, matrix_gsva_kegg)
                #tmp_df = tmp_df.assign(**{'abs(r)': tmp_df['r'].abs()})   # TODO: stop using abs r, use r correl only
                tmp_df['tissue'] = SAMPLE
                correls_dfs.append(tmp_df)

                tmp_variances_df = return_variances_of_scores(
                    kpca_kegg, matrix_gsva_kegg)
                tmp_variances_df['tissue'] = SAMPLE
                scores_each_dfs.append(tmp_variances_df)
            # end if
        except Exception as e:
            print(SAMPLE, "failed pearson")


    # -- save average scores
    scores_avgs_merged = pd.concat(scores_each_dfs, axis=0, ignore_index=True)
    scores_avgs_merged.to_csv("scores_avg_merged.csv", index=False)

    # -- save correlations
    # # -- correlations between spot-wise values: corr(gsva, kpca)
    correls_merged = pd.concat(correls_dfs, axis=0, ignore_index=True)
    correls_merged.to_csv("correls_merged.csv", index=False)

    # -- now save the moran I scores of the features for each gsva and kpca
    moran_comb = pd.concat(morans_dfs, ignore_index=True)
    moran_comb.to_csv("moran_merged.csv", index=False)

