
import os
import pandas as pd
import numpy as np
import anndata as ad

import seaborn as sns
import matplotlib.pyplot as plt


palette_methods = {'kpca': "green", 'gsva': "cadetblue" }


def select_line_log_kpca(out_kpca_logfile,
                         tissue, set_type='KEGG'):
    kpca_time = 'not found'
    with open(out_kpca_logfile, 'r') as f:
        log_kpca_run = f.readlines()
    for l in log_kpca_run:
        if (set_type in l) and (tissue in l) and (
                'computed in seconds' in l):
            kpca_time = l.split(':')[-1].strip()

    return kpca_time


out_gsva_reps = "../out_reports_scripts2"
out_kpca_logfiles = ["copy_log_remote/26-03-21_18-30-11-328845.log",
                     "copy_log_remote/26-03-21_19-07-52-280348.log"
                     ]
data_orig = "../../apollo-data"

input_file1 = "../the_data_list.tsv"
df_data = pd.read_csv(input_file1,
        header=0, index_col=None, sep='\t')

df_time = df_data.copy()
df_time = df_time.assign(gsva_time=np.repeat(np.nan, df_time.shape[0]))
df_time = df_time.assign(kpca_time=np.repeat(np.nan, df_time.shape[0]))

for tissue in df_data['dataset_name'].tolist():
    print(tissue)
    with open(os.path.join(out_gsva_reps, f'{tissue}-gsva.txt'), 'r') as f:
        # gsva time
        gsva_report_lines = f.readlines()
        line_time = [x for x in gsva_report_lines if 'Execution time in seconds' in x]
        gsva_time = line_time[0].split(":")[-1].strip()
        print(gsva_time, "gsva!!")
        df_time.loc[df_time['dataset_name'] == tissue, 'gsva_time'] = float(gsva_time)


        # kpca time
        for out_kpca_logfile in out_kpca_logfiles:
            try:
                kpca_time = select_line_log_kpca(out_kpca_logfile, tissue)
                print(kpca_time, "kpca!!")
                df_time.loc[
                    df_time['dataset_name'] == tissue, 'kpca_time'] = float(kpca_time)
            except:
                continue


# add spots number
df_time = df_time.assign(spots_n=np.nan)

for tissue in df_data['dataset_name'].tolist():
    adata_mz = ad.read_h5ad(os.path.join(data_orig, tissue, f'{tissue}.h5ad'))

    df_time.loc[
        df_time['dataset_name'] == tissue, 'spots_n'] = adata_mz.obs.shape[0]

print(df_time)

#

foo = pd.melt(df_time, id_vars='spots_n', value_vars=['gsva_time', 'kpca_time'],
              var_name='method')

foo.replace('gsva_time', 'gsva', inplace=True)
foo.replace('kpca_time', 'kpca', inplace=True)
foo = foo.assign(method=pd.Categorical(foo['method'].tolist(),
                                       categories=['kpca', 'gsva']))
print(foo.columns)
foo = foo.assign(value=foo['value'].astype(float).to_numpy())
foo['spots_n'] = foo['spots_n'].astype(float).to_numpy()
foo = foo.sort_values(by=['value'] , ascending=True)
sns.lineplot(foo, x='spots_n', y='value', hue='method', style='method',
             markers=True, palette=palette_methods)
plt.ylabel("time")
plt.show()



# trash
# df = df_time[['gsva_time', 'kpca_time', 'spots_n']]
# sns.lineplot(data=df, x='spots_n', y='gsva_time', color='purple')
# sns.lineplot(data=df, x='spots_n', y='kpca_time', color='cadetblue')
# plt.show()