rule prepare_alleleome:
    input:
        roary="data/interim/roary/{name}"
    output:
        pangene_v2="data/processed/{name}/alleleome/pangene_v2.csv",
    log:
        "logs/prepare_alleleome/{name}.log"
    conda:
        "../envs/alleleome.yaml"
    shell:
        """
        python workflow/scripts/alleleome_reassign_pangene_categories.py \
            --data-dir {input.roary} \
            --output-file {output.pangene_v2} 2>> {log}
        """

rule prepare_alleleome_fasta:
    input:
        roary_path="data/interim/roary/{name}",
        gbk_folder=lambda wildcards: expand("data/interim/processed-genbank/{strains}.gbk",
            name=wildcards.name,
            strains=[s for s in PEP_PROJECTS[wildcards.name].sample_table.genome_id.unique()],
        ),
        pangene_summary_path="data/processed/{name}/alleleome/pangene_v2.csv"
    output:
        fasta=directory("data/processed/{name}/alleleome/pangenome_alignments")
    log:
        "logs/prepare_alleleome_fasta/{name}.log"
    conda:
        "../envs/alleleome.yaml"
    shell:
        """
        python workflow/scripts/alleleome_get_core_genes_fasta.py \
            --roary_path {input.roary_path} \
            --gbk_folder data/interim/processed-genbank \
            --pangene_summary_path {input.pangene_summary_path} \
            --output_folder data/processed/{wildcards.name}/alleleome 2>> {log}
        """

rule pan_alleleome:
    input:
        pangenome_aligments=directory("data/processed/{name}/alleleome/pangenome_alignments"),
        alleleome=directory("data/processed/{name}/alleleome/alleleome")
    output:
        pan_alleleome=directory("data/processed/{name}/alleleome/pan_alleleome")
    conda:
        "../envs/alleleome.yaml"
    log:
        "logs/pan_alleleome/{name}.log"
    shell:
        """
        PAN_Alleleome --path1 {input.pangenome_aligments} \
            --path2 {input.alleleome} 2>> {log}
        """
