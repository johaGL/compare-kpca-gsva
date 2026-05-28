#!/bin/bash

samples=$(cat ../the_data_list.tsv | cut -f 1 )

for SAMPLE in ${samples[@]}; do
    if [ $SAMPLE != "dataset_name" ]; then
        echo $SAMPLE
        Rscript gsva-spata2.R "../../compar_out_2026" $SAMPLE "KEGG" > ../out_reports_scripts2/$SAMPLE-gsva.txt
        echo ""
    fi
        
done

echo "when slurm job, make sure head is no longer when all"

#
