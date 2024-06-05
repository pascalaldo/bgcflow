import argparse
import logging
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def remove_special_char(s):
    if "/" in str(s):
        return s.replace("/", "_")
    elif "'" in str(s):
        return s.replace("'", "_variant")
    elif "(" in str(s):
        s = s.replace("(", "_")
        return s.replace(")", "")
    else:
        return s


def load_presence_data(gene_presence_binary_path, gene_presence_locustag_path):
    """
    Load data from Roary output.

    Parameters:
    roary_path (Path): Path to Roary output

    Returns:
    tuple: A tuple containing two pandas DataFrames
    """
    df_gene_presence_binary = pd.read_csv(
        gene_presence_binary_path, index_col="Gene", low_memory=False
    )
    df_gene_presence_locustag = pd.read_csv(
        gene_presence_locustag_path, index_col="Gene", low_memory=False
    )
    df_gene_presence_locustag.index = [
        remove_special_char(str(i)) for i in list(df_gene_presence_locustag.index)
    ]
    return df_gene_presence_binary, df_gene_presence_locustag

def load_locustag_data(all_locustag_path):
    df_all_locustag = pd.read_csv(
        all_locustag_path, index_col=0, low_memory=False
    )
    return df_all_locustag

def parse_genbank_files(df_gene_presence_locustag, gbk_folder):
    """
    Parse GenBank files.

    Parameters:
    df_gene_presence_locustag (DataFrame): DataFrame with gene presence information
    gbk_folder (Path): Path to folder containing GenBank files

    Returns:
    DataFrame: DataFrame with parsed GenBank data
    """
    all_locustag_list = []
    for genome_id in list(df_gene_presence_locustag.columns):
        genbank_file_path = gbk_folder / f"{genome_id}.gbk"
        for record in SeqIO.parse(genbank_file_path, "genbank"):
            for feature in record.features:
                genome_data_list = []
                tag = feature.qualifiers.get("locus_tag")
                if tag:
                    genome_data_list.append(tag[0])  # Locus tag
                    genome_data_list.append(genome_id)  # Genome ID
                    if "product" in feature.qualifiers.keys():
                        genome_data_list.append(
                            feature.qualifiers["product"][0]
                        )  # Prokka annotation
                    else:
                        genome_data_list.append("")
                    genome_data_list.append(
                        int(feature.location.start)
                    )  # Start position
                    genome_data_list.append(int(feature.location.end))  # End position
                    genome_data_list.append(
                        str(feature.extract(record.seq))
                    )  # Nucleotide Seq
                    if "translation" in feature.qualifiers.keys():
                        genome_data_list.append(
                            feature.qualifiers["translation"][0]
                        )  # Amino Acid Seq
                    else:
                        genome_data_list.append("")
                    all_locustag_list.append(genome_data_list)
    all_locustag_df = pd.DataFrame(
        all_locustag_list,
        columns=[
            "Locus_Tag",
            "Genome_ID",
            "Prokka_Annotation",
            "Start_Position",
            "End_Position",
            "Nucleotide_Seq",
            "Amino_Acid_Seq",
        ],
    )
    all_locustag_df.index = all_locustag_df["Locus_Tag"]
    return all_locustag_df

def process_gene(
    gene_id, df_gene_presence_locustag, all_locustag_df, fna_path, faa_path
):
    """
    Process genes and save output files.

    Parameters:
    gene_id (str): Gene name
    df_gene_presence_locustag (DataFrame): DataFrame with gene presence information
    all_locustag_df (DataFrame): DataFrame with all locus tags
    fna_path (Path): .fna output file
    faa_path (Path): .faa output file
    """

    logging.info(f"   Processing gene: {gene_id}")
    
    gene_locustag = []
    for locus_tag_str in df_gene_presence_locustag.loc[gene_id, :].dropna():
        gene_locustag.extend(
            locus_tag_str.split("\t") if "\t" in locus_tag_str else [locus_tag_str]
        )
    nucleotide_records = []
    amino_acid_records = []
    gene_locustag_set = set(gene_locustag)
    filtered_df = all_locustag_df[all_locustag_df.index.isin(gene_locustag_set)]
    nucleotide_records = [
        SeqRecord(
            Seq(row["Nucleotide_Seq"]),
            id=locustag,
            description=f"{row.get('Prokka_Annotation', 'Unknown')} | {row['Genome_ID']}",
        )
        for locustag, row in filtered_df.iterrows()
    ]
    amino_acid_records = [
        SeqRecord(
            Seq(row["Amino_Acid_Seq"]),
            id=locustag,
            description=f"{row.get('Prokka_Annotation', 'Unknown')} | {row['Genome_ID']}",
        )
        for locustag, row in filtered_df.iterrows()
    ]
    with open(fna_path, "w") as fasta_n_file:
        SeqIO.write(nucleotide_records, fasta_n_file, "fasta")
    with open(faa_path, "w") as fasta_aa_file:
        SeqIO.write(amino_acid_records, fasta_aa_file, "fasta")
    logging.info(f"   Finished gene: {gene_id}")

def main():
    parser = argparse.ArgumentParser(description="Process some files.")
    parser.add_argument('mode', type=str, choices=["locustags", "collect"])
    parser.add_argument(
        "--gp_binary", type=str, required=True, help="Path to gene_presence_binary csv file."
    )
    parser.add_argument(
        "--gp_locustag", type=str, required=True, help="Path to gene_presence_locustag csv file."
    )
    parser.add_argument(
        "--gbk_folder", type=str, required=False, help="Folder containing GenBank files (use with mode 'locustags')"
    )
    parser.add_argument(
        "--all_locustag", type=str, required=True, help="Path to all_locustags csv file."
    )
    parser.add_argument(
        "--fna", type=str, required=False, help="Path to pangenes.fna (use with mode 'collect')"
    )
    parser.add_argument(
        "--faa", type=str, required=False, help="Path to pangenes.faa (use with mode 'collect')"
    )
    parser.add_argument(
        "--gene_id", type=str, required=False, help="Gene id to process (use with mode 'collect')"
    )
    parser.add_argument(
        "--which", type=str, default="Core", required=False, help="Which genes to process (Core, Accessory, Pan/All)"
    )

    args = parser.parse_args()

    gp_binary_path = Path(args.gp_binary)
    gp_locustag_path = Path(args.gp_locustag)
    all_locustag_path = Path(args.all_locustag)

    df_gene_presence_binary, df_gene_presence_locustag = load_presence_data(gp_binary_path, gp_locustag_path)
    if args.mode == "locustags":
        gbk_folder = Path(args.gbk_folder)
        df_all_locustag = parse_genbank_files(df_gene_presence_locustag, gbk_folder)
        df_all_locustag.write_csv(all_locustag_path)
    else:
        df_all_locustag = load_locustag_data(all_locustag_path)
        if mode == "collect":
            faa_path = Path(args.faa)
            fna_path = Path(args.fna)
            process_gene(gene_id, df_gene_presence_locustag, df_all_locustag, fna_path, faa_path)
    # gene_list = get_genes(pangene_summary_path, which=which_genes)
    # process_genes(
    #     gene_list, df_gene_presence_locustag, df_all_locustag, output_folder
    # )


if __name__ == "__main__":
    main()
