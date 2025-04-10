"""
Generate JSON output for PanKB from the phylons results
"""

import argparse
import json
import pandas as pd
import os
import logging

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--L", required=True, help="Path of L matrix CSV")
    parser.add_argument("--L_binarized", required=True, help="Path of L binarized matrix CSV")
    parser.add_argument("--A", required=True, help="Path of A matrix CSV")
    parser.add_argument("--A_binarized", required=True, help="Path of A binarized matrix CSV")
    parser.add_argument("--output", required=True, help="Path of the output JSON file")

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(asctime)s - %(levelname)s] - %(message)s")
    logger = logging.getLogger()
    logger.info("Converting phylons to JSON data")
    logger.info(f"Args: {args}")

    N_ROUNDING_DIGITS = 5

    # Load input files
    logger.info("Reading NMF matrices")
    L = pd.read_csv(args.L, index_col=0)
    L_binarized = pd.read_csv(args.L_binarized, index_col=0)
    A = pd.read_csv(args.A, index_col=0)
    A_binarized = pd.read_csv(args.A_binarized, index_col=0)

    logger.info("Making genome to phylons mapping")
    genome_to_phylons = {}
    phylons_all = A_binarized.index
    for genome in A_binarized.columns:
        A_binarized_column = A_binarized[genome]
        phylon_mask = A_binarized_column == 1
        if phylon_mask.sum() < 1:
            genome_to_phylons[genome] = None
            continue

        genome_phylons = [int(phylon) for phylon in phylons_all[phylon_mask].values]
        genome_to_phylons[genome] = genome_phylons

    logger.info("Making genome to phylon weights mapping")
    genome_to_phylon_weights = {}
    for genome in A.columns:
        A_column = A[genome]
        genome_to_phylon_weights[genome] = {
            phylon: round(weight, N_ROUNDING_DIGITS) for phylon, weight in A_column.items()
        }

    logger.info("Making gene to phylons mapping")
    gene_to_phylons = {}
    gene_to_phylon_weights = {}
    for (_, row), (_, row_bin) in zip(L.iterrows(), L_binarized.iterrows()):
        gene = row.name

        phylon_mask = row_bin == 1
        phylons = row.index[phylon_mask]

        for phylon in row.index:
            gene_to_phylon_weights[gene] = {
                phylon: round(float(weight), N_ROUNDING_DIGITS) for phylon, weight in row.items()
            }
            gene_to_phylons[gene] = phylons.tolist()

    output_data = {
        "gene_phylons": gene_to_phylons,
        "gene_phylon_weights": gene_to_phylon_weights,
        "genome_phylons": genome_to_phylons,
        "genome_phylon_weights": genome_to_phylon_weights,
    }

    # Save output files
    logger.info("Saving output file")
    out_dir = os.path.dirname(args.output)
    os.makedirs(out_dir, exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(output_data, f, indent=2)

    logger.info("Done")
