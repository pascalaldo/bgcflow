checkpoint select_genomes:
    input:
        tsv="data/interim/{stage}/ncbi_datasets/{taxon}.tsv",
    output:
        genome_list="data/interim/{stage}/ncbi_datasets/taxon/{taxon}.genome_list",
    run:
        import pandas as pd
        df = pd.read_csv(input.tsv, sep="\t", header=0, index_col=0, low_memory=False)
        genomes = df.index.to_list()
        with open(output.genome_list, "w") as f:
            for genome in genomes:
                f.write(f"{genome}\n")

# def get_accessions_for_taxon(taxon):
#     genome_list = checkpoints.select_genomes.get(taxon=taxon).output.genome_list
#     with open(genome_list, "r") as f:
#         accessions = [g.strip() for g in f.readlines()]
#     return accessions
# def get_all_accessions():
#     accessions = []
#     for taxon in TAXONS.index.to_list():
#         accessions.extend(get_accessions_for_taxon(taxon))
#     return accessions
# def get_taxon_for_accession(accession):
#     for taxon in TAXONS.index.to_list():
#         genome_list = checkpoints.select_genomes.get(taxon=taxon).output.genome_list
#         with open(genome_list, "r") as f:
#             for genome in f:
#                 if genome.strip() == accession:
#                     return taxon
#     raise ValueError(f"No taxon found for {accession}")
