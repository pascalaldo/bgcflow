rule roary:
    input:
        gff=lambda wildcards: get_prokka_outputs(wildcards.name, filter_samples_qc(wildcards, DF_SAMPLES)),
    output:
        core_alignment_header="data/interim/roary/{name}/core_alignment_header.embl",
        pan_genome_reference="data/interim/roary/{name}/pan_genome_reference.fa",
        accessory_header="data/interim/roary/{name}/accessory.header.embl",
        clustered_proteins="data/interim/roary/{name}/clustered_proteins",
        core_gene_alignment="data/interim/roary/{name}/core_gene_alignment.aln",
        gene_presence_absence="data/interim/roary/{name}/gene_presence_absence.csv",
        summary_statistics="data/interim/roary/{name}/summary_statistics.txt",
        accessory_tab="data/interim/roary/{name}/accessory.tab",
        accessory_binary_genes="data/interim/roary/{name}/accessory_binary_genes.fa",
        core_accessory_header="data/interim/roary/{name}/core_accessory.header.embl",
        numer_of_genes_in_pan_genome="data/interim/roary/{name}/number_of_genes_in_pan_genome.Rtab",
        accessory_binary_genes_newick="data/interim/roary/{name}/accessory_binary_genes.fa.newick",
        core_accessory="data/interim/roary/{name}/core_accessory.tab",
        number_of_new_genes="data/interim/roary/{name}/number_of_new_genes.Rtab",
        accessory_graph="data/interim/roary/{name}/accessory_graph.dot",
        core_accessory_graph="data/interim/roary/{name}/core_accessory_graph.dot",
        number_of_unique_genes="data/interim/roary/{name}/number_of_unique_genes.Rtab",
    conda:
        "../envs/roary.yaml"
    params:
        i=80,
        g=80000,
        roary_dir="data/interim/roary/{name}/",
    threads: 16
    log:
        "logs/roary/roary-{name}.log",
    shell:
        """
        rm -rf {params.roary_dir}
        roary -p {threads} -f {params.roary_dir} -i {params.i} -g {params.g} -e -n -r -v {input.gff} &>> {log}
        """

checkpoint roary_reassign_pangene_categories:
    input:
        gene_presence_binary="data/interim/roary/{name}/df_gene_presence_binary.csv",
        pangene_summary="data/interim/roary/{name}/df_pangene_summary.csv"
    output:
        pangene="data/processed/{name}/tables/df_roary_pangene_summary_reassigned.csv",
    log:
        "logs/roary_reassign_pangene_categories/{name}.log"
    conda:
        "../envs/bgc_analytics.yaml"
    shell:
        """
        python workflow/scripts/alleleome_reassign_pangene_categories.py \
            --gp_binary {input.gene_presence_binary} \
            --summary {input.pangene_summary} \
            --output-file {output.pangene} 2>> {log}
        """

rule eggnog_roary:
    input:
        faa="data/interim/roary/{name}/pan_genome_reference.fa",
        eggnog_db="resources/eggnog_db",
        dmnd="resources/eggnog_db/bacteria.dmnd",
    output:
        annotations="data/interim/eggnog_roary/{name}/{name}.emapper.annotations",
        xlsx="data/interim/eggnog_roary/{name}/{name}.emapper.annotations.xlsx",
        hits="data/interim/eggnog_roary/{name}/{name}.emapper.hits",
        seed_orthologs="data/interim/eggnog_roary/{name}/{name}.emapper.seed_orthologs",
        temp_dir=temp(directory("data/interim/eggnog_roary/tmp/{name}"))
    params:
        eggnog_dir="data/interim/eggnog_roary/{name}/",
    conda:
        "../envs/eggnog.yaml"
    threads: 8
    log:
        "logs/eggnog-roary/eggnog-{name}.log",
    shell:
        """
        mkdir -p {output.temp_dir}
        emapper.py -i {input.faa} --translate --itype "CDS" --excel --cpu {threads} -o {wildcards.name} --output_dir {params.eggnog_dir} --data_dir {input.eggnog_db} --temp_dir {output.temp_dir} &>> {log}
        """

rule eggnog_roary_result_copy:
    input:
        annotations="data/interim/eggnog_roary/{name}/{name}.emapper.annotations",
        xlsx="data/interim/eggnog_roary/{name}/{name}.emapper.annotations.xlsx",
        pangene="data/processed/{name}/tables/df_roary_pangene_summary_reassigned.csv",
    output:
        eggnog_xlsx="data/processed/{name}/eggnog_roary/eggnog_roary.xlsx",
        eggnog_annotations="data/processed/{name}/eggnog_roary/emapper.annotations"
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/eggnog-roary/eggnog-result-copy-{name}.log",
    shell:
        """
        cp {input.xlsx} {output.eggnog_xlsx}
        cp {input.annotations} {output.eggnog_annotations}
        """

rule deeptfactor_roary:
    input:
        pan_genome_reference="data/interim/roary/{name}/pan_genome_reference.fa",
        resource="resources/deeptfactor/",
    output:
        df_deeptfactor="data/processed/{name}/tables/df_deeptfactor_roary.tsv",
    conda:
        "../envs/deeptfactor.yaml"
    threads: 8
    log:
        "logs/deeptfactor-roary/deeptfactor-roary-{name}.log",
    params:
        outdir="data/interim/deeptfactor_roary/{name}",
    shell:
        """
        workdir=$PWD
        mkdir -p data/processed/{wildcards.name} 2>> {log}
        (cd {input.resource} && python tf_running.py \
            -i $workdir/{input.pan_genome_reference.fa} -o $workdir/{params.outdir} \
            -g cpu -cpu {threads}) 2>> {log}
        cp {params.outdir}/prediction_result.txt {output.df_deeptfactor} &>> {log}
        """

rule diamond_roary:
    input:
        faa="data/interim/roary/{name}/pan_genome_reference.fa",
        resource="resources/deeptfactor/",
    output:
        diamond_interim="data/interim/diamond/{name}/{name}_pangenome.dmnd",
        diamond_processed="data/processed/{name}/diamond/{name}_pangenome.dmnd",
    conda:
        "../envs/antismash.yaml"
    threads: 8
    log:
        "logs/diamond-roary/diamond-roary-{name}.log",
    shell:
        """
        diamond makedb --in {input.faa} -d {output.diamond_interim} -p {threads} &>> {log}
        cp {output.diamond_interim} {output.diamond_processed} &>> {log}
        """


rule roary_out:
    input:
        summary_statistics="data/interim/roary/{name}/summary_statistics.txt",
        genomes_tree="data/processed/{name}/automlst_wrapper/df_genomes_tree.csv",
        final_tree="data/processed/{name}/automlst_wrapper/final.newick",
    output:
        gene_presence_binary="data/processed/{name}/roary/df_gene_presence_binary.csv",
        gene_presence_locustag="data/processed/{name}/roary/df_gene_presence_locustag.csv",
        summary="data/processed/{name}/roary/df_pangene_summary.csv",
        summary_interim="data/interim/roary/{name}/df_pangene_summary.csv",
        gene_presence_binary_interim="data/interim/roary/{name}/df_gene_presence_binary.csv",
        gene_presence_locustag_interim="data/interim/roary/{name}/df_gene_presence_locustag.csv",
    params:
        roary_interim_dir="data/interim/roary/{name}/",
        roary_processed_dir="data/processed/{name}/roary",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/roary/roary-out-{name}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/make_pangenome_dataset.py {params.roary_interim_dir} {params.roary_processed_dir} {input.final_tree} {input.genomes_tree} 2>> {log}
        """
