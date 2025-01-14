import os
import json
from pathlib import Path
import argparse
from BCBio import GFF
from Bio import SeqIO
import pandas as pd
from sortedcontainers import SortedDict

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
    with open(output_path, "w") as f_output:
        for i_accession, accession in enumerate(accessions):
            print(f"Accession {i_accession + 1}/{len(accessions)} ({accession})", flush=True)
            mapping = []
            feature_dict = SortedDict()
            if input_mode == "ncbi_gff":
                accession_file = ncbi_gff_path / f"{accession}.gff"
                if os.stat(accession_file).st_size > 0:
                    with open(accession_file, "r") as f:
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
                                    # if locus_tag == "BSU_33810":
                                    #     print("### BSU_33810")
                                    #     print(start_pos)
                                    #     print(end_pos)
                                    #     print(strand)
                                    #     print(locus_tag)
                                    #     print(gene_key)
                                    #     print(gene)
                                    #     print(synonyms)
                                    #     print(feature_key)
                                    feature_present = False
                                    if feature_key[0] in feature_dict:
                                        for f in feature_dict[feature_key[0]]:
                                            if f[0] == feature_key:
                                                feature_present = True
                                                break
                                    if feature_present:
                                        print(f"Feature already present for key {feature_key}")
                                    else:
                                        if not feature_key[0] in feature_dict:
                                            feature_dict[feature_key[0]] = []
                                        feature_dict[feature_key[0]].append((feature_key, (locus_tag, feature_type, gene, synonyms)))
                                except:
                                    print(f"Could not properly parse feature {feature.id}")
                                    # print(feature)
                                    # raise
                else:
                    print(f"Skipping, because file {accession_file} is empty.")
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

                        feature_present = False
                        if feature_key[0] in feature_dict:
                            for f in feature_dict[feature_key[0]]:
                                if f[0] == feature_key:
                                    feature_present = True
                                    break
                        if feature_present:
                            print(f"Feature already present for key {feature_key}")
                        else:
                            if not feature_key[0] in feature_dict:
                                feature_dict[feature_key[0]] = []
                            feature_dict[feature_key[0]].append((feature_key, (locus_tag, feature_type, gene, synonyms)))
                    except:
                        print(f"Error when parsing locus_tag {locus_tag}.")

            print("Built feature_dict", flush=True)
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
                    exact_match = False
                    for f in feature_dict.get(feature_key[0], []):
                        if f[0] == feature_key:
                            original_info = f[1]
                            exact_match = True
                            break
                    if not exact_match:
                        # gene_basename = gene.split("_", 1)[0]
                        min_diff = None
                        for sp in feature_dict.irange(start_pos - 300, start_pos + 300):
                            for f in feature_dict[sp]:
                                k, v = f
                                if k[2] != strand:
                                    continue
                                # if v[2] != gene_basename and not gene_basename in v[3]:
                                #     continue
                                f_start_pos = k[0]
                                f_end_pos = k[1]
                                start_diff = abs(f_start_pos - start_pos)
                                end_diff = abs(f_end_pos - end_pos)
                                if (end_pos - f_start_pos) > 30 and (f_end_pos - f_start_pos) > 30 and (start_diff % 3) == 0 and (end_diff % 3) == 0 and ((start_diff < 210 and end_diff < 210) or (min(start_diff, end_diff) < 21 and max(start_diff, end_diff) < 300)):
                                    if min_diff is None or min_diff >= (start_diff + end_diff):
                                        feature_key = k
                                        original_info = v
                                    # feature_key = k
                                    # original_info = feature_dict[feature_key]
                                    # exact_match = False

                    # if abs(start_pos - 3468253) < 3 or abs(end_pos - 3468253) < 3:
                    #     print("$$$ 3468253")
                    #     print(start_pos)
                    #     print(end_pos)
                    #     print(strand)
                    #     print(feature_key)
                    #     print(locus_tag)
                    #     print(gene)
                    #     print(original_info)

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
            first_write = (df is None)

            df = pd.DataFrame.from_records(mapping)
            if first_write:
                df.to_csv(f_output, index=False)
            else:
                df.to_csv(f_output, index=False, header=False)

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
