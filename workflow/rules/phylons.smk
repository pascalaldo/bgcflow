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
        output="data/processed/{stage}/pankb/web_data/species/{name}/phylons.json",
    conda:
        "../envs/pyphylon.yaml"
    log:
        "logs/{stage}/{name}/make_phylons_json.log"
    shell:
        """
        python3 workflow/scripts/make_phylons_json.py \
            --L {input.L} \
            --L_binarized {input.L_binarized} \
            --A {input.A} \
            --A_binarized {input.A_binarized} \
            --output {output.output} \
            &> {log}
        """

rule phylons_all:
    input:
        fexpand(
            ["data/processed/{stage}/pankb/web_data/species/{name}/phylons.json"],
            stage=RULE_FUNCTIONS["pankb_data_prep"]["stages"],
            name=RULE_FUNCTIONS["pankb_data_prep"]["projects"],
        )
