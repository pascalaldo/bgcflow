import os
import json
from pathlib import Path
import argparse
from BCBio import GFF
from Bio import SeqIO
import pandas as pd

def initialize_parser(parser):
    parser.description = "Process data required for genome pages."
    parser.add_argument(
        "--samples",
        type=str,
        help="Table with genome accessions under the column 'genome_id'.",
    )
    parser.add_argument(
        "--ncbi_gff",
        type=str,
        required=False,
        default=None,
        help="GFF file from NCBI (with the original locus_tags).",
    )
    parser.add_argument(
        "--imodulon_table",
        type=str,
        required=False,
        default=None,
        help="iModulon gene annotation file."
    )
    parser.add_argument(
        "--prokka_gbk",
        type=str,
        required=True,
        help="Prokka GBK file (with the prokka locus_tags).",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Output file with mapping.",
    )

def touch(path):
    with open(path, 'a'):
        os.utime(path, None)

def locus_tag_mapping(
    samples_path,
    ncbi_gff_path,
    imodulon_table_path,
    prokka_gbk_path,
    output_path,
):
    try:
        accessions = pd.read_csv(samples_path, usecols=["genome_id"])
        accessions = accessions["genome_id"].to_list()
    except pd.errors.EmptyDataError:
        accessions = []

    prokka_gbk_path = Path(prokka_gbk_path)

    if not ncbi_gff_path is None and imodulon_table_path is None:
        input_mode = "ncbi_gff"
        ncbi_gff_path = Path(ncbi_gff_path)
    elif not imodulon_table_path is None and ncbi_gff_path is None:
        input_mode = "imodulon_table"
        imodulon_table_path = Path(imodulon_table_path)
    else:
        raise Exception("Specify one of ncbi_gff and imodulon_table")

    df = None

    for accession in accessions:
        mapping = []
        feature_dict = {}
        if input_mode == "ncbi_gff":
            with open(ncbi_gff_path / f"{accession}.gff", "r") as f:
                for record in GFF.parse(f):
                    for feature in record.features:
                        try:
                            feature_type = feature.type
                            if not feature_type in ["gene", "pseudogene"]:
                                continue
                            start_pos = int(feature.location.start)
                            end_pos = int(feature.location.end)
                            strand = feature.location.strand
                            locus_tag = feature.qualifiers["locus_tag"][0]
                            gene_key = "gene" if "gene" in feature.qualifiers else "Name"
                            gene = feature.qualifiers[gene_key][0]
                            synonyms = set()
                            if len(feature.qualifiers[gene_key]) > 1:
                                for synonym in feature.qualifiers[gene_key][1:]:
                                    synonyms.add(synonym)
                            for synonym in feature.qualifiers.get("gene_synonym", []):
                                synonyms.add(synonym)
                            synonyms = list(synonyms)
                            feature_key = (start_pos, end_pos, strand)
                            if feature_key in feature_dict:
                                print(f"Feature already present for key {feature_key}")
                            else:
                                feature_dict[feature_key] = (locus_tag, feature_type, gene, synonyms)
                        except:
                            print(f"Could not properly parse feature {feature.id}")
                            # print(feature)
                            # raise
        elif input_mode == "imodulon_table":
            imod_table = pd.read_csv(imodulon_table_path / accession / "gene_info.csv", index_col=0)
            imod_table.index.name = "locus_tag"
            imod_table.rename(columns={"stop": "end"}, inplace=True)
            for locus_tag, row in imod_table.iterrows():
                try:
                    feature_type = "gene"
                    start_pos = int(row["start"]) -1
                    end_pos = int(row["end"])
                    strand = int(str(row["strand"]).replace("+", "1").replace("-", "-1"))
                    gene = row["gene_name"]
                    synonyms = []
                    feature_key = (start_pos, end_pos, strand)
                    if feature_key in feature_dict:
                        print(f"Feature already present for key {feature_key}")
                    else:
                        # print(f"Adding feature: {feature_key}: {(locus_tag, feature_type, gene, synonyms)}")
                        feature_dict[feature_key] = (locus_tag, feature_type, gene, synonyms)
                except:
                    print(f"Error when parsing locus_tag {locus_tag}.")

        
        for record in SeqIO.parse(prokka_gbk_path / f"{accession}.gbk", "genbank"):
            for i, feature in enumerate(record.features):
                feature_type = feature.type
                if not feature_type == "CDS":
                    continue
                start_pos = int(feature.location.start)
                end_pos = int(feature.location.end)
                strand = feature.location.strand
                feature_key = (start_pos, end_pos, strand)
                locus_tag = feature.qualifiers["locus_tag"][0]
                gene = feature.qualifiers.get("gene", [None])[0]
                
                original_info = None
                if feature_key in feature_dict:
                    original_info = feature_dict[feature_key]
                    exact_match = True
                elif not gene is None:
                    # gene_basename = gene.split("_", 1)[0]
                    for k, v in feature_dict.items():
                        # if v[2] != gene_basename and not gene_basename in v[3]:
                        #     continue
                        start_diff = abs(k[0] - start_pos)
                        end_diff = abs(k[1] - end_pos)
                        if k[2] == strand and (start_diff % 3) == 0 and (end_diff % 3) == 0 and ((start_diff < 210 and end_diff < 210) or (min(start_diff, end_diff) < 21 and max(start_diff, end_diff) < 900)):
                            feature_key = k
                            original_info = feature_dict[feature_key]
                            exact_match = False
                            break

                if original_info is None:
                    # print(feature_key)
                    # print(f"Could not find feature {feature}")
                    # raise Exception()
                    pass
                else:
                    mapping.append(
                        {
                            "prokka_locus_tag": locus_tag,
                            "genome": accession,
                            "prokka_gene": gene,
                            "prokka_start_pos": start_pos,
                            "prokka_end_pos": end_pos,
                            "prokka_strand": strand,
                            "original_locus_tag": original_info[0],
                            "original_type": original_info[1],
                            "original_gene": original_info[2],
                            "original_synonyms": ",".join(original_info[3]),
                            "original_start_pos": feature_key[0],
                            "original_end_pos": feature_key[1],
                            "original_strand": feature_key[2],
                            "exact_match": exact_match,
                        }
                    )
        
        df = pd.concat([df, pd.DataFrame.from_records(mapping, index=["prokka_locus_tag", "genome"])])
    if df is None:
        touch(output_path)
    else:
        df.to_csv(output_path)

def run(args):
    locus_tag_mapping(
        args.samples,
        args.ncbi_gff,
        args.imodulon_table,
        args.prokka_gbk,
        args.output,
    )


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
