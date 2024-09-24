checkpoint merge_taxonomy:
    input:
        gtdb="data/interim/{stage}/gtdb/{taxon}/tables/df_gtdb_meta.csv",
        gtdbtk="data/processed/{stage}/{taxon}/tables/gtdbtk.bac120.summary.tsv",
        # gtdbtk="data/interim/{stage}/gtdbtk/{taxon}/result/classify/gtdbtk.bac120.summary.tsv",
    output:
        meta=report(
            "data/processed/{stage}/{taxon}/tables/df_gtdb_meta.csv",
            caption="../report/table-gtdb.rst",
            category="{taxon}",
            subcategory="Taxonomic Placement",
        ),
    run:
        import pandas as pd
        gtdb_df = pd.read_csv(input.gtdb, header=0, index_col=0, low_memory=False)
        gtdbtk_df = pd.read_csv(input.gtdbtk, header=0, index_col=0, low_memory=False, sep="\t")

        clas_levels = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
        gtdb_merge_df = gtdb_df[clas_levels + ["Organism"]].copy()
        gtdb_merge_df["classification_source"] = "gtdb"

        gtdbtk_merge_df = pd.DataFrame(index=gtdbtk_df.index)
        gtdbtk_merge_df.index.name = "genome_id"

        def get_clas(values, c):
            l = values.split(";")
            prefix = f"{c[0].lower()}__"
            for v in l:
                if v.startswith(prefix):
                    return v
            return None
        
        for c in clas_levels:
            gtdbtk_merge_df[c] = gtdbtk_df["classification"].apply(lambda x: get_clas(x, c))
        gtdbtk_merge_df["Organism"] = gtdbtk_merge_df["Species"]
        gtdbtk_merge_df["Species"] = gtdbtk_merge_df["Species"].str.split(" ", n=1).str[1]
        gtdbtk_merge_df.fillna("", inplace=True)
        # gtdbtk_merge_df["Species"] = gtdbtk_merge_df["Species"].str.removeprefix("s__")
        # gtdbtk_merge_df["Organism"] = "s__" + gtdbtk_merge_df["Genus"].str.removeprefix("g__") + " " + gtdbtk_merge_df["Species"]
        gtdbtk_merge_df["classification_source"] = "gtdbtk"

        gtdb_merge_df.update(gtdbtk_merge_df)
        gtdb_merge_df.to_csv(output.meta)