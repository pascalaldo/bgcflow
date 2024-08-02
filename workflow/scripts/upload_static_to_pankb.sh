#!/bin/bash

for d in data/processed/*/ ; do
    SPECIES=`basename $d`
    ALLELEOME_DIR="${d}alleleome/"
    echo $SPECIES
    # final_dn_ds_count_per_gene.csv
    for f in dn_ds.json step_line.json ; do 
        azcopy copy ${ALLELEOME_DIR}Pan/$f https://pankb.blob.core.windows.net/data/PanKB/web_data/species/$SPECIES/panalleleome/$PANKB_WRITE_SAS
    done
    # azcopy copy --recursive ${ALLELEOME_DIR}gene_data https://pankb.blob.core.windows.net/data/PanKB/web_data/species/$SPECIES/panalleleome/$PANKB_WRITE_SAS
done
