#!/bin/bash
orig_dir="../../data"
tissues=$(cat ../the_data_list.tsv | cut -f 1 )

echo ${tissues[@]}
destination="../../compar_out_2026"

for u in ${tissues[@]}; do
  if [ $u != "dataset_name" ]; then
    mkdir -p $destination/$u;
    cp $orig_dir/$u/$u".h5ad" $destination/$u;
    cp "$orig_dir/$u/$u-"*.csv $destination/$u;
    cp "$orig_dir/$u/$u"*CleanAnnotations.tsv $destination/$u;
  fi
done

cp -r "$orig_dir/metabolites-mappings/" $destination
