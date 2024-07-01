checkpoint preselect_genomes:
    input:
        tsv="data/interim/ncbi_datasets/taxon/{taxon}_original.tsv",
    output:
        tsv="data/interim/ncbi_datasets/taxon/{taxon}.tsv",
    shell:
        """
            cp {input.tsv} {output.tsv}
        """
