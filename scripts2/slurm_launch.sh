#!/bin/bash
########################## Slurm options
#SBATCH --job-name=gsva2026
#SBATCH --output=/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/slurmout/gsva2026_%j.out
#SBATCH --chdir=/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/
#SBATCH --mail-user=juana7@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=180G
#SBATCH --ntasks-per-node=20
#SBATCH --time=30:00:00
##################################################

# run: sbatch -p compute_long <SCRIPT>
scontrol show job $SLURM_JOB_ID

eval "$(conda shell.bash hook)"  # this solves commandnotfounderror when conda activate


module load r_gsva

cd /mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/compare_2026/scripts2

samples=$(cat ../the_data_list.tsv | cut -f 1 )

for SAMPLE in ${samples[@]}; do
    if [ $SAMPLE != "dataset_name" ]; then
        echo $SAMPLE
        Rscript gsva-spata2.R "../../compar_out_2026" $SAMPLE "KEGG" > ../out_reports_scripts2/$SAMPLE-gsva.txt
        echo ""
    fi
        
done

