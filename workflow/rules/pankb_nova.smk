# rule eggnog_roary:
#     input:
#         faa="data/interim/{stage}/roary/{name}/pan_genome_reference.fa",
#         eggnog_db="resources/eggnog_db",
#         dmnd="resources/eggnog_db/bacteria.dmnd",
#     output:
#         annotations="data/interim/{stage}/eggnog_roary/{name}/{name}.emapper.annotations",
#         xlsx="data/interim/{stage}/eggnog_roary/{name}/{name}.emapper.annotations.xlsx",
#         hits="data/interim/{stage}/eggnog_roary/{name}/{name}.emapper.hits",
#         seed_orthologs="data/interim/{stage}/eggnog_roary/{name}/{name}.emapper.seed_orthologs",
#         temp_dir=temp(directory("data/interim/{stage}/eggnog_roary/tmp/{name}"))
#     params:
#         eggnog_dir="data/interim/{stage}/eggnog_roary/{name}/",
#     conda:
#         "../envs/eggnog.yaml"
#     threads: 8
#     log:
#         "logs/{stage}/eggnog-roary/eggnog-{name}.log",
#     shell:
#         """
#         mkdir -p {output.temp_dir}
#         emapper.py -i {input.faa} --translate --itype "CDS" --excel --cpu {threads} -o {wildcards.name} --output_dir {params.eggnog_dir} --data_dir {input.eggnog_db} --temp_dir {output.temp_dir} &>> {log}
#         """

# rule eggnog_roary_result_copy:
#     input:
#         annotations="data/interim/{stage}/eggnog_roary/{name}/{name}.emapper.annotations",
#         xlsx="data/interim/{stage}/eggnog_roary/{name}/{name}.emapper.annotations.xlsx",
#     output:
#         eggnog_xlsx="data/processed/{stage}/{name}/eggnog_roary/eggnog_roary.xlsx",
#         eggnog_annotations="data/processed/{stage}/{name}/eggnog_roary/emapper.annotations"
#     conda:
#         "../envs/bgc_analytics.yaml"
#     log:
#         "logs/{stage}/eggnog-roary/eggnog-result-copy-{name}.log",
#     shell:
#         """
#         cp {input.xlsx} {output.eggnog_xlsx}
#         cp {input.annotations} {output.eggnog_annotations}
#         """

rule pankb_nova_organism:
    input:
        gp_binary="data/interim/{stage}/roary/{name}/df_gene_presence_binary.csv",
        summary_v2="data/interim/{stage}/alleleome/{name}/pangene_v2.csv",
        gtdb_meta="data/processed/{stage}/{name}/tables/df_gtdb_meta.csv",
        filt_norm="data/processed/{stage}/{name}/alleleome/Pan/final_pan_aa_thresh_core_genes_dom_var_genome_count_pos_normalized.csv",
        sel_genes="data/interim/{stage}/alleleome/{name}/sel_genes.csv",
    output:
        organism="data/processed/{stage}/pankb/web_data/species/{name}/nova/organism.json",
    log:
        "logs/{stage}/pankb_nova/pankb_nova_organism_{name}.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_nova organism {wildcards.name} \
            --gp_binary {input.gp_binary} \
            --summary {input.summary_v2} \
            --gtdb_meta {input.gtdb_meta} \
            --filt_norm {input.filt_norm} \
            --sel_genes {input.sel_genes} \
            -o {output.organism} > {log} 2>&1
        """

rule pankb_nova_pangene:
    input:
        gp_binary="data/interim/{stage}/roary/{name}/df_gene_presence_binary.csv",
        eggnog_table="data/processed/{stage}/{name}/eggnog_roary/emapper.annotations",
        summary_v2="data/interim/{stage}/alleleome/{name}/pangene_v2.csv",
        gtdb_meta="data/processed/{stage}/{name}/tables/df_gtdb_meta.csv",
    output:
        pangene="data/processed/{stage}/pankb/web_data/species/{name}/nova/pangene.json",
        pathway="data/processed/{stage}/pankb/web_data/species/{name}/nova/pathway.json",
    log:
        "logs/{stage}/pankb_nova/pankb_nova_pangene_{name}.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_nova pangene {wildcards.name} \
            --gp_binary {input.gp_binary} \
            --summary {input.summary_v2} \
            --gtdb_meta {input.gtdb_meta} \
            --eggnog_table {input.eggnog_table} \
            -o {output.pangene} \
            -p {output.pathway} > {log} 2>&1
        """

rule pankb_nova_gene:
    input:
        gp_locustag="data/interim/{stage}/roary/{name}/df_gene_presence_locustag.csv",
        all_locustag="data/interim/{stage}/alleleome/{name}/all_locustag.csv",
        gtdb_meta="data/processed/{stage}/{name}/tables/df_gtdb_meta.csv",
        dummy="data/interim/{stage}/alleleome/{name}/pangenome_alignments/fasta_dummy_Pan",
    output:
        gene="data/processed/{stage}/pankb/web_data/species/{name}/nova/gene.jsonl.gz",
    params:
        fasta_dir="data/interim/{stage}/alleleome/{name}/pangenome_alignments/input/"
    log:
        "logs/{stage}/pankb_nova/pankb_nova_gene_{name}.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_nova gene {wildcards.name} \
            --gp_locustag {input.gp_locustag} \
            --all_locustag {input.all_locustag} \
            --fasta_dir {params.fasta_dir} \
            --gtdb_meta {input.gtdb_meta} \
            -o {output.gene} > {log} 2>&1
        """


rule pankb_nova_all:
    input:
        fexpand(
            [
                "data/processed/{stage}/pankb/web_data/species/{name}/nova/organism.json",
                "data/processed/{stage}/pankb/web_data/species/{name}/nova/pangene.json",
                "data/processed/{stage}/pankb/web_data/species/{name}/nova/pathway.json",
                "data/processed/{stage}/pankb/web_data/species/{name}/nova/gene.jsonl.gz",
            ],
            stage=RULE_FUNCTIONS["pankb_data_prep"]["stages"],
            name=RULE_FUNCTIONS["pankb_data_prep"]["projects"],
        ),
