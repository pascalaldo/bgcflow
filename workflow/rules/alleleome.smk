import pandas as pd

rule alleleome_prepare:
    input:
        gp_binary="data/interim/roary/{name}/df_gene_presence_binary.csv",
        gp_locustag="data/interim/roary/{name}/df_gene_presence_locustag.csv",
        summary="data/interim/roary/{name}/df_pangene_summary.csv",
        gbk_files=lambda wildcards: get_prokka_outputs(wildcards.name, filter_samples_qc(wildcards, DF_SAMPLES), ext="gbk", path="processed-genbank"),
    output:
        summary_v2="data/processed/{name}/alleleome/pangene_v2.csv",
        all_locustag="data/processed/{name}/alleleome/all_locustag.csv",
        all_genes="data/processed/{name}/alleleome/all_genes.csv",
        sel_locustag="data/processed/{name}/alleleome/sel_locustag.csv",
        sel_genes="data/processed/{name}/alleleome/sel_genes.csv",
    params:
        gbk_folder="data/interim/processed-genbank/",
    log:
        "logs/alleleome/prepare_{name}.log"
    conda:
        "../envs/alleleome.yaml"
    shell:
        """
        alleleome prepare \
            --gp_binary {input.gp_binary} \
            --gp_locustag {input.gp_locustag} \
            --summary {input.summary} \
            --gbk_folder {params.gbk_folder} \
            --summary_v2 {output.summary_v2} \
            --all_locustag {output.all_locustag} \
            --all_genes {output.all_genes} \
            --sel_locustag {output.sel_locustag} \
            --sel_genes {output.sel_genes} > {log} 2>&1
        """

rule alleleome_fasta:
    input:
        all_locustag="data/processed/{name}/alleleome/all_locustag.csv",
        all_genes="data/processed/{name}/alleleome/all_genes.csv",
        sel_locustag="data/processed/{name}/alleleome/sel_locustag.csv",
        sel_genes="data/processed/{name}/alleleome/sel_genes.csv",
    output:
        dummy="data/processed/{name}/alleleome/pangenome_alignments/fasta_dummy_{pan_core}",
        gene_list="data/processed/{name}/alleleome/{pan_core}/gene_list.txt",
    params:
        out_dir="data/processed/{name}/alleleome/pangenome_alignments/",
        pan_core_flag=lambda wildcards: ("--pan" if wildcards.pan_core == "Pan" else "--no-pan"),
    log:
        "logs/alleleome/fasta_{name}_{pan_core}.log"
    conda:
        "../envs/alleleome.yaml"
    shell:
        """
        alleleome fasta \
            --all_locustag {input.all_locustag} \
            --all_genes {input.all_genes} \
            --sel_locustag {input.sel_locustag} \
            --sel_genes {input.sel_genes} \
            --gene_list {output.gene_list} \
            --out_dir {params.out_dir} \
            {params.pan_core_flag} > {log} 2>&1
        touch {output.dummy}
        """

rule alleleome_process:
    input:
        gene_list="data/processed/{name}/alleleome/{pan_core}/gene_list.txt",
        dummy="data/processed/{name}/alleleome/pangenome_alignments/fasta_dummy_{pan_core}",
    output:
        dummy="data/processed/{name}/alleleome/pangenome_alignments/process_dummy_{pan_core}",
    params:
        out_dir="data/processed/{name}/alleleome/pangenome_alignments/",
    threads: workflow.cores
    log:
        "logs/alleleome/process_{name}_{pan_core}.log"
    conda:
        "../envs/alleleome.yaml"
    shell:
        """
        alleleome process \
            --gene_list {input.gene_list} \
            --out_dir {params.out_dir} \
            -p {threads} > {log} 2>&1
        touch {output.dummy}
        """

rule alleleome_analyze:
    input:
        gene_list="data/processed/{name}/alleleome/{pan_core}/gene_list.txt",
        dummy="data/processed/{name}/alleleome/pangenome_alignments/process_dummy_{pan_core}",
    output:
        aa_vars="data/processed/{name}/alleleome/{pan_core}/pan_amino_acid_vars_df.csv",
        codon_muts="data/processed/{name}/alleleome/{pan_core}/pan_gene_syno_non_syno_df.csv",
    params:
        out_dir="data/processed/{name}/alleleome/pangenome_alignments/",
    log:
        "logs/alleleome/analyze_{name}_{pan_core}.log"
    conda:
        "../envs/alleleome.yaml"
    shell:
        """
        alleleome analyze \
            --gene_list {input.gene_list} \
            --out_dir {params.out_dir} \
            --aa_vars {output.aa_vars} \
            --codon_muts {output.codon_muts} > {log} 2>&1
        """

rule alleleome_preplot:
    input:
        gene_list="data/processed/{name}/alleleome/{pan_core}/gene_list.txt",
        aa_vars="data/processed/{name}/alleleome/{pan_core}/pan_amino_acid_vars_df.csv",
    output:
        dominant_aa="data/processed/{name}/alleleome/{pan_core}/final_core_consensus_dominant_aa_count_df.csv",
        variable_aa="data/processed/{name}/alleleome/{pan_core}/final_core_pan_aa_thresh_vars_all_substitutions_sep_df.csv",
        dom_var="data/processed/{name}/alleleome/{pan_core}/final_pan_aa_thresh_core_genes_dominant_variant_genome_count_pos.csv",
        gaps="data/processed/{name}/alleleome/{pan_core}/pan_aa_thresh_core_genes_aa_pos_with_gaps.csv",
        filt_norm="data/processed/{name}/alleleome/{pan_core}/final_pan_aa_thresh_core_genes_dom_var_genome_count_pos_normalized.csv",
        dom_var_out_dir=directory("data/processed/{name}/alleleome/{pan_core}/dom_var/"),
    params:
        out_dir="data/processed/{name}/alleleome/pangenome_alignments/",
    log:
        "logs/alleleome/preplot_{name}_{pan_core}.log"
    conda:
        "../envs/alleleome.yaml"
    shell:
        """
        alleleome preplot \
            --gene_list {input.gene_list} \
            --aa_vars {input.aa_vars} \
            --out_dir {params.out_dir} \
            --dominant_aa {output.dominant_aa} \
            --variable_aa {output.variable_aa} \
            --dom_var {output.dom_var} \
            --gaps {output.gaps} \
            --filt_norm {output.filt_norm} \
            --dom_var_out_dir {output.dom_var_out_dir} \
            > {log} 2>&1
        """
