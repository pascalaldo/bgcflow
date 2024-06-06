import pandas as pd

def get_genes(pangene_summary_path, which="Core"):
    """
    Get list of genes.

    Parameters:
    pangene_summary_path (Path): Path to pangenome summary file
    which (string): "Core", "Accessory", "All"/"Pan"

    Returns:
    list: List of genes
    """
    gene_class_table = pd.read_csv(pangene_summary_path, index_col=0)
    if which in ["All", "Pan"]:
        gene_list = list(gene_class_table.index)
    else:
        gene_list = list(
            gene_class_table.loc[gene_class_table["pangenome_class_2"] == which, :].index
        )
    gene_list = [str(i).replace("/", "") for i in gene_list]
    return gene_list

checkpoint prepare_alleleome_summary:
    input:
        gene_presence_binary="data/interim/roary/{name}/df_gene_presence_binary.csv",
        pangene_summary="data/interim/roary/{name}/df_pangene_summary.csv"
    output:
        pangene_v2="data/processed/{name}/alleleome/pangene_v2.csv",
    log:
        "logs/prepare_alleleome/{name}.log"
    conda:
        "../envs/alleleome.yaml"
    shell:
        """
        python workflow/scripts/alleleome_reassign_pangene_categories.py \
            --gp_binary {input.gene_presence_binary} \
            --summary {input.pangene_summary} \
            --output-file {output.pangene_v2} 2>> {log}
        """

rule prepare_alleleome_genes:
    input:
        gene_presence_binary="data/interim/roary/{name}/df_gene_presence_binary.csv",
        gene_presence_locustag="data/interim/roary/{name}/df_gene_presence_locustag.csv",
        gbk_folder="data/interim/processed-genbank/"
        # gbk_files=lambda wildcards: expand("data/interim/processed-genbank/{strains}.gbk",
        #     name=wildcards.name,
        #     strains=[s for s in PEP_PROJECTS[wildcards.name].sample_table.genome_id.unique()],
        # )
    output:
        all_locustags="data/interim/alleleome/{name}/df_all_locustag.csv"
    log:
        "logs/prepare_alleleome_fasta/{name}_locustags.log"
    conda:
        "../envs/alleleome.yaml"
    shell:
        """
        python workflow/scripts/alleleome_prepare.py locustags \
            --gp_binary {input.gene_presence_binary} \
            --gp_locustag {input.gene_presence_locustag} \
            --all_locustag {output.all_locustags} \
            --gbk_folder {input.gbk_folder} 2>> {log}
        """

rule alleleome_collect_pangene:
    input:
        gene_presence_binary="data/interim/roary/{name}/df_gene_presence_binary.csv",
        gene_presence_locustag="data/interim/roary/{name}/df_gene_presence_locustag.csv",
        all_locustags="data/interim/alleleome/{name}/df_all_locustag.csv"
    output:
        fna="data/processed/{name}/alleleome/pangenome_alignments/{gene}/input/pangenes.fna",
        faa="data/processed/{name}/alleleome/pangenome_alignments/{gene}/input/pangenes.faa"
    log:
        "logs/prepare_alleleome_fasta/{name}_{gene}_collect.log"
    conda:
        "../envs/alleleome.yaml"
    shell:
        """
        python workflow/scripts/alleleome_prepare.py collect \
            --gp_binary {input.gene_presence_binary} \
            --gp_locustag {input.gene_presence_locustag} \
            --all_locustag {input.all_locustags} \
            --fna {output.fna} \
            --faa {output.faa} \
            --gene_id {wildcards.gene} 2>> {log}
        """

# rule prepare_alleleome_fasta:
#     input:
#         roary_path="data/interim/roary/{name}",
#         gbk_folder=lambda wildcards: expand("data/interim/processed-genbank/{strains}.gbk",
#             name=wildcards.name,
#             strains=[s for s in PEP_PROJECTS[wildcards.name].sample_table.genome_id.unique()],
#         ),
#         pangene_summary_path="data/processed/{name}/alleleome/pangene_v2.csv"
#     output:
#         fasta=directory("data/processed/{name}/alleleome/pangenome_alignments")
#     log:
#         "logs/prepare_alleleome_fasta/{name}.log"
#     conda:
#         "../envs/alleleome.yaml"
#     shell:
#         """
#         python workflow/scripts/alleleome_get_core_genes_fasta.py \
#             --roary_path {input.roary_path} \
#             --gbk_folder data/interim/processed-genbank \
#             --pangene_summary_path {input.pangene_summary_path} \
#             --which Core \
#             --output_folder data/processed/{wildcards.name}/alleleome 2>> {log}
#         """

rule alleleome:
    input:
        pangenome_aligments_fna=lambda wildcards: expand("data/processed/{{name}}/alleleome/pangenome_alignments/{gene}/input/pangenes.fna", gene=get_genes(checkpoints.prepare_alleleome_summary.get(name=wildcards.name).output.pangene_v2, which=wildcards.pan_core)),
        pangenome_aligments_faa=lambda wildcards: expand("data/processed/{{name}}/alleleome/pangenome_alignments/{gene}/input/pangenes.faa", gene=get_genes(checkpoints.prepare_alleleome_summary.get(name=wildcards.name).output.pangene_v2, which=wildcards.pan_core)),
        pangene_summary_path="data/processed/{name}/alleleome/pangene_v2.csv"
    output:
        alleleome_directory_path="data/processed/{name}/alleleome/{pan_core}/pan_gene_syno_non_syno_df.csv"
    conda:
        "../envs/alleleome.yaml"
    log:
        "logs/alleleome/{name}_{pan_core}.log"
    shell:
        """
        Alleleome {wildcards.pan_core} --path1 data/processed/{wildcards.name}/alleleome/pangenome_alignments/ \
            --path2 data/processed/{wildcards.name}/alleleome/{wildcards.pan_core}/ \
            --table {input.pangene_summary_path} \
            --log_to_terminal 2>> {log}
        """
