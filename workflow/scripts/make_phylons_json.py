"""
Generate JSON output for PanKB from the phylons results
"""

import argparse
import json
from collections import defaultdict
import pandas as pd
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--L", required=True, help="Path of L matrix CSV")
    parser.add_argument("--L_binarized", required=True, help="Path of L binarized matrix CSV")
    parser.add_argument("--A", required=True, help="Path of A matrix CSV")
    parser.add_argument("--A_binarized", required=True, help="Path of A binarized matrix CSV")
    parser.add_argument("--genome_to_phylons", required=True, help="Path of output genome to phylons JSON")
    parser.add_argument("--phylon_to_genomes", required=True, help="Path of output phylon to genome JSON")
    parser.add_argument("--phylon_to_genes", required=True, help="Path of output phylon to genes JSON")
    parser.add_argument("--phylon_to_gene_weights", required=True, help="Path of output phylon to gene weights JSON")
    parser.add_argument("--gene_to_phylons", required=True, help="Path of output gene to phylons JSON")
    parser.add_argument("--gene_to_phylon_weights", required=True, help="Path of output gene to phylon weights JSON")

    args = parser.parse_args()

    # Load input files
    L = pd.read_csv(args.L, index_col=0)
    L_binarized = pd.read_csv(args.L_binarized, index_col=0)
    # A = pd.read_csv(args.A, index_col=0)
    A_binarized = pd.read_csv(args.A_binarized, index_col=0)

    genome_to_phylons = {}
    phylons = A_binarized.index
    for genome in A_binarized.columns:
        A_column = A_binarized[genome]
        phylon_mask = A_column == 1
        if phylon_mask.sum() < 1:
            genome_to_phylons[genome] = None
            continue

        phylons = [int(phylon) for phylon in phylons[phylon_mask].values]
        genome_to_phylons[genome] = phylons

    phylon_to_genomes = defaultdict(list)
    for genome, phylons in genome_to_phylons.items():
        for phylon in phylons:
            phylon_to_genomes[phylon].append(genome)
    phylon_to_genomes = dict(phylon_to_genomes)

    phylon_to_gene_weights = {}
    phylon_to_genes = {}
    genes = L.index
    for phylon in L.columns:
        column = L[phylon]
        column_bin = L_binarized[phylon]

        genes_mask = column_bin == 1
        phylon_genes = genes[genes_mask]

        phylon_to_gene_weights[phylon] = {gene: round(weight, 5) for gene, weight in column.items()}
        phylon_to_genes[phylon] = list(phylon_genes)

    gene_to_phylons = {}
    gene_to_phylon_weights = {}
    for (_, row), (_, row_bin) in zip(L.iterrows(), L_binarized.iterrows()):
        gene = row.name

        phylon_mask = row_bin == 1
        phylons = row.index[phylon_mask]

        for phylon in row.index:
            gene_to_phylon_weights[gene] = {phylon: round(float(weight), 5) for phylon, weight in row.items()}
            gene_to_phylons[gene] = phylons.tolist()

    # Save output files
    output_files = [
        args.genome_to_phylon,
        args.phylon_to_genomes,
        args.phylon_to_genes,
        args.phylon_to_gene_weights,
        args.gene_to_phylons,
        args.gene_to_phylon_weights,
    ]
    for otuput_fp in output_files:
        os.makedirs(os.path.dirname(otuput_fp), exist_ok=True)

    with open(args.genome_to_phylon, "w", encoding="utf-8") as f:
        json.dump(genome_to_phylons, f, indent=2)

    with open(args.phylon_to_genomes, "w", encoding="utf-8") as f:
        json.dump(phylon_to_genomes, f, indent=2)

    with open(args.phylon_to_genes, "w", encoding="utf-8") as f:
        json.dump(phylon_to_genes, f, indent=2)

    with open(args.phylon_to_gene_weights, "w", encoding="utf-8") as f:
        json.dump(phylon_to_gene_weights, f, indent=2)

    with open(args.gene_to_phylons, "w", encoding="utf-8") as f:
        json.dump(gene_to_phylons, f, indent=2)

    with open(args.gene_to_phylon_weights, "w", encoding="utf-8") as f:
        json.dump(gene_to_phylon_weights, f, indent=2)
