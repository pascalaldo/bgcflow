checkpoint pankb_imodulon_organism_data:
    input:
        samples="data/processed/species/samples/{name}.csv",
    output:
        genome_list="data/interim/{stage}/pankb_imodulon/{name}/genome_list.csv",
        output_dir=directory("data/interim/{stage}/pankb_imodulon/{name}/data/"),
    log:
        "logs/{stage}/pankb_imodulon/pankb_imodulon_organism_data_{name}.log"
    conda:
        "../envs/pankb_data_prep.yaml"
    shell:
        """
        pankb_imodulon organism_data {wildcards.name} \
            --samples {input.samples} \
            --genome_list {output.genome_list} \
            -o {output.output_dir} > {log} 2>&1
        """

def get_imodulon_organism_genomes(wildcards):
    genome_list_path = checkpoints.pankb_imodulon_organism_data.get(stage=wildcards.stage, name=wildcards.name).output["genome_list"]
    try:
        df_genomes = pd.read_csv(genome_list_path)
    except pd.errors.EmptyDataError:
        return []
    genomes = df_genomes["genome_id"].to_list()
    return genomes

rule imodulon_locus_tags:
    input:
        imodulon_table="data/interim/{stage}/pankb_imodulon/{name}/data/",
        prokka_gbk=fexpand("data/interim/{{stage}}/processed-genbank/{accession}.gbk", accession=get_imodulon_organism_genomes),
        samples="data/interim/{stage}/pankb_imodulon/{name}/genome_list.csv",
    output:
        csv="data/interim/{stage}/pankb_imodulon/{name}/locus_tag_mapping.csv",
    params:
        prokka_gbk="data/interim/{stage}/processed-genbank/",
    conda:
        "../envs/locus_tag_mapping.yaml"
    log:
        "logs/{stage}/imodulon_locus_tags/imodulon_locus_tags-{name}.log"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/map_locus_tags.py \
            --samples {input.samples} \
            --imodulon_table {input.imodulon_table} \
            --prokka_gbk {params.prokka_gbk} \
            -o {output.csv} > {log}
        """