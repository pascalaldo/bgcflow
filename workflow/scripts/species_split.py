import pandas as pd
import argparse

def species_split(gtdb_meta_path, qc_passed_path, overview_path, classification_path, min_clas_size):
    df_gtdb = pd.read_csv(gtdb_meta_path, low_memory=False, index_col=0, header=0)
    df_qc = pd.read_csv(qc_passed_path, low_memory=False, index_col=0, header=0)
    df_overview = pd.read_csv(overview_path, low_memory=False, index_col="genome_id", header=0)

    sel_gtdb = set(df_gtdb.index[df_gtdb["Organism"] != "s__"].to_list())
    sel_qc = set(df_qc.index[df_qc["passed"]].to_list())
    sel_overview = set(df_overview.index.to_list())

    sel = list(sel_gtdb & sel_qc & sel_overview)

    df_full = df_overview.loc[sel, :].merge(df_gtdb.loc[sel, :], how='outer', left_index=True, right_index=True)
    df_full["genome_id"] = df_full.index
    df_full["Organism [ ]"] = df_full["Organism"].str.removeprefix("s__")
    df_full["Organism [_]"] = df_full["Organism [ ]"].str.replace(" ", "_")

    df_clas = df_full[["Organism [_]", "genome_id"]].groupby("Organism [_]").agg(list)
    df_clas = df_clas.loc[df_clas["genome_id"].apply(len) >= min_clas_size, :]
    df_clas["genome_id"] = df_clas["genome_id"].apply(lambda x: ",".join(x))
    df_clas.to_csv(classification_path)

def main():
    parser = argparse.ArgumentParser(description="Determine projects based on species info.")
    parser.add_argument(
        "--gtdb_meta", type=str, required=True, help="df_gtdb_meta.csv file"
    )
    parser.add_argument(
        "--qc_passed", type=str, required=True, help="qc_passed.csv file"
    )
    parser.add_argument(
        "--overview", type=str, required=True, help="<taxon>.csv file"
    )
    parser.add_argument(
        "--classification", type=str, required=True, help="classification.csv file"
    )
    parser.add_argument(
        "--min", type=int, required=True, help="Minimum number of genomes in a project."
    )

    args = parser.parse_args()
    species_split(
        gtdb_meta_path=args.gtdb_meta,
        qc_passed_path=args.qc_passed,
        overview_path=args.overview,
        classification_path=args.classification,
        min_clas_size=args.min,
    )

if __name__ == "__main__":
    main()
