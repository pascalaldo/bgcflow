checkpoint merge_taxonomy:
    input:
        gtdb="data/interim/gtdb/{taxon}/tables/df_gtdb_meta.csv",
    output:
        meta=report(
            "data/processed/{taxon}/tables/df_gtdb_meta.csv",
            caption="../report/table-gtdb.rst",
            category="{taxon}",
            subcategory="Taxonomic Placement",
        ),
    shell:
        """
            cp {input.gtdb} {output.meta}
        """