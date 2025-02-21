rule phylons:
    input:
        mash_distances="data/processed/{stage}/{name}/mash/df_mash.csv",
        gene_presence="data/interim/{stage}/roary/{name}/df_gene_presence_binary.csv",
        pangenome_summary="data/interim/{stage}/alleleome/{name}/pangene_v2.csv",
    output:
        L="data/processed/{stage}/{name}/phylons/NMF_L.csv",
        L_binarized="data/processed/{stage}/{name}/phylons/NMF_L_binarized.csv",
        A="data/processed/{stage}/{name}/phylons/NMF_A.csv",
        A_binarized="data/processed/{stage}/{name}/phylons/NMF_A_binarized.csv",
    log:
        "logs/{stage}/{species}/compute_phylons.log"
    conda:
        "../envs/..."
    shell:
        """
        python3 workflow/scripts/compute_phylons.py \
            --mash_distances {input.mash_distances} \
            --gene_presence {input.gene_presence} \
            --pangenome_summary {input.pangenome_summary} \
            &> {log}
        """
