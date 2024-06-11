import pandas as pd

def get_genes(sel_genes_path, which="Core"):
    sel_genes = pd.read_csv(sel_genes_path, index_col=0)
    if which == "All":
        gene_list = sel_genes.index.to_list()
    elif which == "Pan":
        gene_list = sel_genes.index[sel_genes["Pan"]].to_list()
    else:
        gene_list = sel_genes.index[sel_genes["Core"]].to_list()
    gene_list = [str(i).replace("/", "") for i in gene_list]
    return gene_list

checkpoint alleleome_prepare:
    input:
        gp_binary="data/interim/roary/{name}/df_gene_presence_binary.csv",
        gp_locustag="data/interim/roary/{name}/df_gene_presence_locustag.csv",
        summary="data/interim/roary/{name}/df_pangene_summary.csv"
        gbk_folder="data/interim/processed-genbank/"
    output:
        summary_v2="data/processed/{name}/alleleome/pangene_v2.csv"
        all_locustag="data/processed/{name}/alleleome/all_locustag.csv"
        all_genes="data/processed/{name}/alleleome/all_genes.csv"
        sel_locustag="data/processed/{name}/alleleome/sel_locustag.csv"
        sel_genes="data/processed/{name}/alleleome/sel_genes.csv"
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
            --gbk_folder {input.gbk_folder} \
            --summary_v2 {output.summary_v2} \
            --all_locustag {output.all_locustag} \
            --all_genes {output.all_genes} \
            --sel_locustag {output.sel_locustag} \
            --sel_genes {output.sel_genes} > {log} 2>&1
        """


rule alleleome_fasta:
    input:
        all_locustag="data/processed/{name}/alleleome/all_locustag.csv"
        all_genes="data/processed/{name}/alleleome/all_genes.csv"
        sel_locustag="data/processed/{name}/alleleome/sel_locustag.csv"
        sel_genes="data/processed/{name}/alleleome/sel_genes.csv"
    output:
        fna=lambda wildcards: expand("data/processed/{{name}}/alleleome/pangenome_alignments/{gene}/input/pan_genes.fna", gene=get_genes(checkpoints.alleleome_prepare.get(name=wildcards.name).output.sel_genes, which=wildcards.pan_core))
        faa=lambda wildcards: expand("data/processed/{{name}}/alleleome/pangenome_alignments/{gene}/input/pan_genes.faa", gene=get_genes(checkpoints.alleleome_prepare.get(name=wildcards.name).output.sel_genes, which=wildcards.pan_core))
        dummy="data/processed/{name}/alleleome/pangenome_alignments/dummy_{pan_core}"
    params:
        out_dir="data/processed/{name}/alleleome/pangenome_alignments/"
        pan_core_flag=lambda wildcards: ("--pan" if wildcards.pan_core == "Pan" else "--no-pan")
    log:
        "logs/alleleome/fasta_{name}.log"
    conda:
        "../envs/alleleome.yaml"
    shell:
        """
        alleleome fasta \
            --all_locustag {input.all_locustag} \
            --all_genes {input.all_genes} \
            --sel_locustag {input.sel_locustag} \
            --sel_genes {input.sel_genes} \
            --out_dir {params.out_dir} \
            {params.pan_core_flag} > {log} 2>&1
        touch {output.dummy}
        """

rule alleleome_process:
    input:
        fna="data/processed/{name}/alleleome/pangenome_alignments/{gene}/input/pan_genes.fna"
        faa="data/processed/{name}/alleleome/pangenome_alignments/{gene}/input/pan_genes.faa"
    output:
        mafft_na="data/processed/{name}/alleleome/pangenome_alignments/{gene}/output/mafft_nucleotide_{gene}.fasta"
        mafft_aa="data/processed/{name}/alleleome/pangenome_alignments/{gene}/output/mafft_amino_acid_{gene}.fasta"
        consensus_na="data/processed/{name}/alleleome/pangenome_alignments/{gene}/output/nucleotide_consensus_{gene}.fna"
        consensus_aa="data/processed/{name}/alleleome/pangenome_alignments/{gene}/output/amino_acid_consensus_{gene}.faa"
        blast_na="data/processed/{name}/alleleome/pangenome_alignments/{gene}/output/nucleotide_blast_out_{gene}.xml"
        blast_aa="data/processed/{name}/alleleome/pangenome_alignments/{gene}/output/amino_acid_blast_out_{gene}.xml"
    params:
        out_dir="data/processed/{name}/alleleome/pangenome_alignments/{gene}/"
    log:
        "logs/alleleome/process_{name}_{gene}.log"
    conda:
        "../envs/alleleome.yaml"
    shell:
        """
        alleleome process_gene \
            --out_dir {params.out_dir} \
            --gene_id {wildcards.gene} > {log} 2>&1
        """

rule alleleome_analyze:
    input:
        all_genes="data/processed/{name}/alleleome/all_genes.csv"
        sel_locustag="data/processed/{name}/alleleome/sel_locustag.csv"
        sel_genes="data/processed/{name}/alleleome/sel_genes.csv"
        blast_na=lambda wildcards: expand("data/processed/{{name}}/alleleome/pangenome_alignments/{gene}/output/nucleotide_blast_out_{gene}.xml", gene=get_genes(checkpoints.alleleome_prepare.get(name=wildcards.name).output.sel_genes, which=wildcards.pan_core))
        blast_aa=lambda wildcards: expand("data/processed/{{name}}/alleleome/pangenome_alignments/{gene}/output/amino_acid_blast_out_{gene}.xml", gene=get_genes(checkpoints.alleleome_prepare.get(name=wildcards.name).output.sel_genes, which=wildcards.pan_core))
    output:
        aa_vars="data/processed/{name}/alleleome/pan_amino_acid_vars_df.csv"
        codon_muts="data/processed/{name}/alleleome/pan_gene_syno_non_syno_df.csv"
    params:
        out_dir="data/processed/{name}/alleleome/pangenome_alignments/"
        pan_core_flag=lambda wildcards: ("--pan" if wildcards.pan_core == "Pan" else "--no-pan")
    log:
        "logs/alleleome/analyze_{name}.log"
    conda:
        "../envs/alleleome.yaml"
    shell:
        """
        alleleome analyze \
            --all_genes {input.all_genes} \
            --sel_locustag {input.sel_locustag} \
            --sel_genes {input.sel_genes} \
            --out_dir {params.out_dir} \
            --aa_vars {output.aa_vars} \
            --codon_muts {output.codon_muts} \
            {params.pan_core_flag} > {log} 2>&1
        touch {output.dummy}
        """
