rule phylons:
    input:
        mash_distances="data/processed/{stage}/{name}/mash/df_mash.csv",
        gene_presence="data/interim/{stage}/roary/{name}/df_gene_presence_binary.csv",
        pangenome_summary="data/interim/{stage}/alleleome/{name}/pangene_v2.csv",
    output:
        "data/processed/{stage}/{name}/phylons/NMF_L.csv",
        "data/processed/{stage}/{name}/phylons/NMF_L_binarized.csv",
        "data/processed/{stage}/{name}/phylons/NMF_A.csv",
        "data/processed/{stage}/{name}/phylons/NMF_A_binarized.csv",
        out_dir=directory("data/processed/{stage}/{name}/phylons/"),
    log:
        "logs/{stage}/{name}/compute_phylons.log"
    conda:
        "../envs/pyphylon.yaml"
    shell:
        """
        python3 workflow/scripts/compute_phylons.py \
            --mash_distances {input.mash_distances} \
            --gene_presence {input.gene_presence} \
            --pangenome_summary {input.pangenome_summary} \
            --output_dir {output.out_dir} \
            &> {log}
        """

rule pankb_phylons:
    input:
        L="data/processed/{stage}/{name}/phylons/NMF_L.csv",
        L_binarized="data/processed/{stage}/{name}/phylons/NMF_L_binarized.csv",
        A="data/processed/{stage}/{name}/phylons/NMF_A.csv",
        A_binarized="data/processed/{stage}/{name}/phylons/NMF_A_binarized.csv",
    output:
        genome_to_phylon="data/processed/{stage}/pankb/web_data/species/{name}/phylons/genome_to_phylon.json",
        phylon_to_genome="data/processed/{stage}/pankb/web_data/species/{name}/phylons/phylon_to_genome.json",
        gene_to_phylons="data/processed/{stage}/pankb/web_data/species/{name}/phylons/gene_to_phylons.json",
        phylon_to_genes="data/processed/{stage}/pankb/web_data/species/{name}/phylons/phylon_to_genes.json",
        phylon_to_gene_weights="data/processed/{stage}/pankb/web_data/species/{name}/phylons/phylon_to_gene_weights.json",
        gene_to_phylon_weights="data/processed/{stage}/pankb/web_data/species/{name}/phylons/gene_to_phylon_weights.json",
    conda:
        "../envs/pyphylon.yaml"
    shell:
        """
        python3 workflow/scripts/make_phylons_json.py \
            --L {input.L} \
            --L_binarized {input.L_binarized} \
            --A {input.A} \
            --A_binarized {input.A_binarized} \
            --genome_to_phylon {output.genome_to_phylon} \
            --phylon_to_genome {output.phylon_to_genome} \
            --phylon_to_genes {output.phylon_to_genes} \
            --phylon_to_gene_weights {output.phylon_to_gene_weights} \
            --gene_to_phylons {output.gene_to_phylons} \
            --gene_to_phylon_weights {output.gene_to_phylon_weights}
        """
