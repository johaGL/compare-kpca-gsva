# Compare k-PCA vs GSVA
**Across spatially resolved metabolomics/lipidomics datasets**

This comparison was performed using the datasets specified in:

```
the_data_list.tsv
```

The steps were performed in order, always from `compare-kpca-gsva`:
* `cd scripts0` then `./copy_files.sh`: to copy the base .h5ad files and their annotations, plus the **KEGG** pathways
* `cd scripts1` then `module load smpath` then `python3 run.py`: to obtain compound resolved matrix per dataset
* `cd scripts2` then `sbatch slurm_launch.sh`: to run **GSVA**
* `cd scripts3` then  `module load smpath` and finally: 

  
```
python3 compare_outs__gsva_kpca.py
python3 viz_runtime.py
```

to run **k-PCA** (via SpacePath which in turn uses the sspa implementation of k-PCA)

The input matrices (compound resolved) as well as the reference pathways were the same for both methods.
