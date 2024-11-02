def ncbi_genome_list_dummy_input(wildcards):
    RULE_FUNCTIONS["ncbi_annotation_datasets"]["accessions_for_project"](wildcards)
    return []

checkpoint select_genomes_for_annotation_download:
    input:
        dummy=ncbi_genome_list_dummy_input
    output:
        genome_list="data/interim/{stage}/ncbi_annotation_datasets/project/{name}.genome_list",
    run:
        import pandas as pd
        genomes = RULE_FUNCTIONS["ncbi_annotation_datasets"]["accessions_for_project"](wildcards)
        with open(output.genome_list, "w") as f:
            for genome in genomes:
                f.write(f"{genome}\n")

rule ncbi_annotation_dataset_download_dehydrated:
    input:
        genome_list="data/interim/{stage}/ncbi_annotation_datasets/project/{name}.genome_list",
    output:
        dehydrated_dataset="data/interim/{stage}/ncbi_annotation_datasets/projects/{name}.dehydrated.zip",
        readme="data/interim/{stage}/ncbi_annotation_datasets/datasets/{name}/README.md",
        fetch="data/interim/{stage}/ncbi_annotation_datasets/datasets/{name}/ncbi_dataset/fetch.txt",
        assembly_report="data/interim/{stage}/ncbi_annotation_datasets/datasets/{name}/ncbi_dataset/data/assembly_data_report.jsonl",
        dataset_catalog="data/interim/{stage}/ncbi_annotation_datasets/datasets/{name}/ncbi_dataset/data/dataset_catalog.json",
    params:
        project_out_dir="data/interim/{stage}/ncbi_annotation_datasets/datasets/{name}/"
    conda:
        "../envs/ncbi_datasets.yaml"
    log:
        "logs/{stage}/ncbi_annotation_datasets/ncbi_annotation_dataset_download_dehydrated_{name}.log",
    shell:
        """
            datasets download genome accession --inputfile {input.genome_list} --dehydrated --filename {output.dehydrated_dataset} --include gff3 > {log} 2>&1
            unzip -o -d {params.project_out_dir} {output.dehydrated_dataset} >> {log} 2>&1
        """

rule ncbi_annotation_dataset_rehydrate:
    input:
        genome_list="data/interim/{stage}/ncbi_annotation_datasets/project/{name}.genome_list",
        fetch="data/interim/{stage}/ncbi_annotation_datasets/datasets/{name}/ncbi_dataset/fetch.txt",
        assembly_report="data/interim/{stage}/ncbi_annotation_datasets/datasets/{name}/ncbi_dataset/data/assembly_data_report.jsonl",
    output:
        dummy="data/interim/{stage}/ncbi_annotation_datasets/project/{name}.dummy",
    params:
        project_out_dir="data/interim/{stage}/ncbi_annotation_datasets/datasets/{name}/"
    conda:
        "../envs/ncbi_datasets.yaml"
    log:
        "logs/{stage}/ncbi_annotation_datasets/ncbi_annotation_dataset_redydrate_{name}.log",
    shell:
        """
            datasets rehydrate --directory {params.project_out_dir} > {log} 2>&1
            while IFS="" read -r p || [ -n "$p" ]
            do
                if [ -f {params.project_out_dir}/ncbi_annotation_dataset/data/$p/genomic.gff ]
                then
                    ls -w 1 {params.project_out_dir}/ncbi_annotation_dataset/data/$p/genomic.gff | wc -l | grep "^1$" > /dev/null
                    mv {params.project_out_dir}/ncbi_annotation_dataset/data/$p/genomic.gff {params.project_out_dir}/ncbi_annotation_dataset/data/$p/$p.gff
                else
                    echo "File for $p not found."
                fi
            done < {input.genome_list}
            touch {output.dummy}
        """

rule ncbi_annotation_dataset_collect:
    input:
        dummy=fexpand("data/interim/{stage}/ncbi_annotation_datasets/project/{name}.dummy", name=RULE_FUNCTIONS["ncbi_annotation_datasets"]["project_for_accession"], stage=RULE_FUNCTIONS["ncbi_annotation_datasets"]["stages"]),
    output:
        gff="data/interim/all/gff/{accession}.gff",
    params:
        gff=fexpand("data/interim/{stage}/ncbi_annotation_datasets/datasets/{name}/ncbi_dataset/data/{accession}/{accession}.gff", name=RULE_FUNCTIONS["ncbi_annotation_datasets"]["project_for_accession"], accession=(lambda wildcards: wildcards.accession), stage=RULE_FUNCTIONS["ncbi_annotation_datasets"]["stages"]),
    conda:
        "../envs/data_processing.yaml"
    log:
        "logs/ncbi_annotation_datasets/ncbi_dataset_collect_{accession}.log",
    shell:
        """
            if [ -f {params.gff} ]
            then
                if [ ! -f {output.gff} ]
                then
                    RELDIR=`dirname {output.gff}`
                    LINKPATH=`realpath -s --relative-to=$RELDIR "{params.gff}"`
                    ln -s $LINKPATH {output.gff}
                fi
            else
                echo "File {params.gff} not found."
            fi
            touch {output.gff}
        """