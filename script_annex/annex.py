import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq
import seaborn as sns
import anndata as ad
import smpath.processing.measures_methods as spm


SAMPLE = "Stomach-rat"

print(os.listdir(os.path.join( "../../compar_out_2026")))

#   f'{SAMPLE}-{pw_type}-CleanAnnotations.tsv')


scores_oh = ad.read_h5ad(os.path.join("../", SAMPLE, f"{SAMPLE}-scores-KEGG.h5ad"))
print(scores_oh.uns.keys())
print(scores_oh.uns['current_pathways_df'].head())

oo = scores_oh.uns['current_pathways_df']['pathway_id'].unique().tolist()
df_i = pd.DataFrame({
    'pathway_id' : oo,
    'U' : np.random.gamma(shape=0.2, scale=0.8, size=len(oo))
})

foo = pd.merge(df_i, scores_oh.var, how='left', on=['pathway_id'])

print(foo)

print(foo.shape)


print(scores_oh)
scores_oh = spm.compute_moran_i(scores_oh, n_perms=2)

df_i = scores_oh.uns['moranI'].copy()
df_i['pathway_id'] = scores_oh.uns['moranI'].index.tolist()
df_i = df_i.sort_values(by=['I'], ascending=False)
print(df_i[['pathway_id', 'I']].head(8))
print("=====================")
# moran
print(scores_oh.var.head(2))

pws = ['hsa00120', 'hsa02010']
for po in pws:
    sq.pl.spatial_scatter(
        scores_oh, color=[po],
        img=False, size=1
    )
    plt.show()


## --
sco_gsva = pd.read_csv(os.path.join(
         "../../compar_out_2026",
    SAMPLE, f"{SAMPLE}-cmps-mx-KEGG-KEGG_gsva.tsv"
    ), sep='\t', index_col=0)
coords_df = pd.read_csv(os.path.join( "../../compar_out_2026",
                                                 SAMPLE, f'{SAMPLE}-coords.tsv'),
 sep='\t', index_col=0)

from scipy.sparse import csr_matrix
ob_gsva = ad.AnnData(
X = csr_matrix(sco_gsva.values),
var = pd.DataFrame({'pathway_id': np.array(sco_gsva.columns.tolist(), dtype=str)},
                   index=sco_gsva.columns.tolist(), dtype=str),
obs = coords_df,
)
ob_gsva.obsm['spatial'] = np.array(coords_df)
ob_gsva.uns['spatial'] = {'my_slice_name': {
    'scalefactors': {'fiducial_diameter_fullres': 100,
                     'spot_diameter_fullres': 1,
                     'tissue_hires_scalef': 1,
                     'tissue_lowres_scalef': 1}}}

pws = ['hsa00120', 'hsa02010']
for po in pws:
    sq.pl.spatial_scatter(
        ob_gsva, color=[po],
        img=False, size=1
    )
    plt.title("gsva")
    plt.show()