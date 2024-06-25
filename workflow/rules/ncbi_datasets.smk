#datasets summary genome taxon Lactobacillaceae  --assembly-source refseq --as-json-lines | dataformat tsv genome > Lactobacillaceae_refseq.tsv

rule ncbi_dataset_for_taxon:
    output:
        tsv="data/interim/ncbi_datasets/taxon/{taxon}.tsv",
    params:
        taxon=lambda wildcards: wildcards.taxon.replace('_', ' '),
        fields=",".join(["accession", "organism-name", "organism-tax-id", "source_database", "assminfo-name", "assminfo-level"])
    conda:
        "../envs/ncbi_datasets.yaml"
    log:
        "logs/ncbi_datasets/ncbi_dataset_for_taxon_{taxon}.log",
    shell:
        """
            datasets summary genome taxon {params.taxon} --assembly-version latest --assembly-source RefSeq --fields {params.fields} --as-json-lines | dataformat tsv genome > {output.tsv} 2> {log}
        """

rule ncbi_dataset_download_genome_for_taxon_dehydrated:
    input:
        genome_list="data/interim/ncbi_datasets/taxon/{taxon}.genome_list",
    output:
        dehydrated_dataset="data/interim/ncbi_datasets/taxon/{taxon}.dehydrated.zip",
        readme="data/interim/ncbi_datasets/datasets/{taxon}/README.md",
        fetch="data/interim/ncbi_datasets/datasets/{taxon}/ncbi_dataset/fetch.txt",
        assembly_report="data/interim/ncbi_datasets/datasets/{taxon}/ncbi_dataset/data/assembly_data_report.jsonl",
        dataset_catalog="data/interim/ncbi_datasets/datasets/{taxon}/ncbi_dataset/data/dataset_catalog.json",
    params:
        taxon_out_dir="data/interim/ncbi_datasets/datasets/{taxon}/"
    conda:
        "../envs/ncbi_datasets.yaml"
    log:
        "logs/ncbi_datasets/ncbi_dataset_download_genome_for_taxon_dehydrated_{taxon}.log",
    shell:
        """
            datasets download genome accession --inputfile {input.genome_list} --dehydrated --filename {output.dehydrated_dataset} > {log} 2>&1
            unzip -o -d {params.taxon_out_dir} {output.dehydrated_dataset}
        """

rule ncbi_dataset_rehydrate:
    input:
        genome_list="data/interim/ncbi_datasets/taxon/{taxon}.genome_list",
        fetch="data/interim/ncbi_datasets/datasets/{taxon}/ncbi_dataset/fetch.txt",
        assembly_report="data/interim/ncbi_datasets/datasets/{taxon}/ncbi_dataset/data/assembly_data_report.jsonl",
    output:
        dataset_dir=directory("data/interim/ncbi_datasets/datasets/{taxon}/"),
        dummy="data/interim/ncbi_datasets/taxon/{taxon}.dummy",
    conda:
        "../envs/ncbi_datasets.yaml"
    log:
        "logs/ncbi_datasets/ncbi_dataset_redydrate_{taxon}.log",
    shell:
        """
            datasets rehydrate --directory {output.dataset_dir} > {log} 2>&1
        """

rule ncbi_dataset_collect:
    input:
        fna=lambda wildcards: expand("data/interim/ncbi_datasets/datasets/{taxon}/ncbi_dataset/data/{{accession}}/{{accession}}.{{assembly}}_genomic.fna", taxon=get_taxon_for_accession(wildcards.accession))
        jsonl_report=lambda wildcards: expand("data/interim/ncbi_datasets/datasets/{taxon}/ncbi_dataset/data/assembly_data_report.jsonl", taxon=get_taxon_for_accession(wildcards.accession))
    output:
        fna="data/interim/fasta/{accession}.fna",
        json_report="data/interim/assembly_report/{accession}.json",
    conda:
        "../envs/ncbi_datasets.yaml"
    log:
        "logs/ncbi_datasets/ncbi_dataset_collect_{accession}.log",
    shell:
        """
            mv {input.fna} {output.fna}
            grep '^{"accession":"$p"' {input.jsonl_report} > {output.json_report}
        """