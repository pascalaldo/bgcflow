import pandas as pd

checkpoint species_split:
    input:
        gtdb="data/processed/{taxon}/tables/df_gtdb_meta.csv",
        qc="data/processed/qc/qc_passed.csv",
        overview="data/interim/ncbi_datasets/taxon/{taxon}.csv",
    output:
        classification="data/interim/split/{taxon}/classification.csv",
    params:
        min_clas_size=30,
    conda:
        "../envs/data_processing.yaml"
    log: "log/split/species_split_{taxon}.log",
    shell:
        """
            python workflow/scripts/species_split.py \
                --gtdb_meta {input.gtdb} \
                --qc_passed {input.qc} \
                --overview {input.overview} \
                --classification {output.classification} \
                --min {params.min_clas_size} > {log} 2>&1
        """

checkpoint extract_species_split_samples:
    input:
        classification=lambda wildcards: f"data/interim/split/{get_taxon_for_species_project(wildcards.name)}/classification.csv",
    output:
        samples="data/processed/samples/{name}.csv",
    run:
        import pandas as pd
        taxon = get_taxon_for_species_project(wildcards.name)
        df_samples = get_species_projects_samples_df_for_taxon(taxon)
        df_samples = df_samples.loc[df_samples["name"] == wildcards.name, :]
        df_samples.to_csv(output.samples)
        # df_samples.to_csv(output.samples_2)

def get_samples_for_species_project(name):
    df = pd.read_csv(checkpoints.extract_species_split_samples.get(name=name).output.samples, header=0, index_col=0, low_memory=False)
    return df.index.to_list()
def get_samples_df_for_species_project(name):
    df = pd.read_csv(checkpoints.extract_species_split_samples.get(name=name).output.samples, header=0, index_col=0, low_memory=False)
    return df

rule extract_species_split_all_samples:
    input:
        classification=expand("data/interim/split/{taxon}/classification.csv", taxon=TAXONS.index.to_list()),
    output:
        samples="data/interim/bgcflow_utils/samples.csv",
    run:
        df = get_species_projects_samples_df()
        df.to_csv(output.samples)

rule extract_species_split_gtdb_meta:
    input:
        gtdb=lambda wildcards: f"data/processed/{get_taxon_for_species_project(wildcards.name)}/tables/df_gtdb_meta.csv",
        samples="data/processed/samples/{name}.csv",
    output:
        gtdb="data/processed/{name,[A-Za-z]+_[A-Za-z]+}/tables/df_gtdb_meta.csv",
    run:
        import pandas as pd
        df_gtdb = pd.read_csv(input.gtdb, header=0, index_col=0, low_memory=False)
        df_samples = pd.read_csv(input.samples, header=0, index_col=0)
        df_gtdb = df_gtdb.loc[df_samples.index.to_list(), :]
        df_gtdb.to_csv(output.gtdb)

rule extract_species_split_overview:
    input:
        overview=lambda wildcards: f"data/interim/ncbi_datasets/taxon/{get_taxon_for_species_project(wildcards.name)}.csv",
        samples="data/processed/samples/{name}.csv",
    output:
        overview="data/processed/samples/overview_{name, [A-Za-z]+_[A-Za-z]+}.csv",
    run:
        import pandas as pd
        df_overview = pd.read_csv(input.overview, header=0, index_col=0, low_memory=False)
        df_samples = pd.read_csv(input.samples, header=0, index_col=0)
        overview = overview.loc[df_samples.index.to_list(), :]
        df_overview.to_csv(output.overview)

def get_species_projects_for_taxon(taxon):
    df = pd.read_csv(checkpoints.species_split.get(taxon=taxon).output.classification, index_col=0, header=0)
    return df.index.to_list()
def get_species_projects():
    projects = []
    for taxon in TAXONS.index.to_list():
        projects.extend(get_species_projects_for_taxon(taxon))
    return projects

def get_species_projects_samples_df_for_taxon(taxon):
    df = pd.read_csv(checkpoints.species_split.get(taxon=taxon).output.classification, index_col=0, header=0)
    df["genome_id"] = df["genome_id"].apply(lambda x: x.split(","))
    df.reset_index(drop=False, inplace=True)
    df = df.explode("genome_id")
    df.set_index("genome_id", drop=True, inplace=True)
    df.rename(columns={"Organism [_]": "name"}, inplace=True)
    df["source"] = "ncbi"
    df["taxon"] = taxon
    return df
def get_species_projects_samples_df():
    df = None
    for taxon in TAXONS.index.to_list():
        df = pd.concat([df, get_species_projects_samples_df_for_taxon(taxon)])
    return df

def get_taxon_for_species_project(name):
    for taxon in TAXONS.index.to_list():
        projects = get_species_projects_for_taxon(taxon)
        if name in projects:
            return taxon
    return None
