import json
import os
import sys
from shutil import copyfile

sys.setrecursionlimit(3000) # TODO: Temporary fix. Biopython uses recursion, which fails for large trees

import pandas as pd
from Bio import Phylo


def copy_automlst_out(automlst_interim_folder, automlst_processed_folder):
    """
    Copy important files from autoMLST interim to proecessed directory

    Parameters
    ----------
    1. automlst_interim_folder : str / path
        Location of the output from autoMLST in the interim directory

    Returns
    -------
    1. automlst_processed_folder : str / path
        Location of the processed output directory for important autoMLST results
    """

    files_to_copy = [
        "raxmlpart.txt",
        "raxmlpart.txt.bionj",
        "raxmlpart.txt.contree",
        "raxmlpart.txt.iqtree",
        "raxmlpart.txt.log",
        "raxmlpart.txt.mldist",
        "raxmlpart.txt.treefile",
    ]

    for tree_file in files_to_copy:
        from_path = os.path.join(automlst_interim_folder, tree_file)
        to_path = os.path.join(automlst_processed_folder, tree_file)

        if os.path.isfile(from_path):
            copyfile(from_path, to_path)

    # Selected tree newick file for downstream analysis by default automlst_wrapper
    treefile = os.path.join(automlst_interim_folder, "raxmlpart.txt.treefile")
    out_newick_processed = os.path.join(automlst_processed_folder, "final.newick")
    out_newick_interim = os.path.join(automlst_interim_folder, "final.newick")

    if os.path.isfile(treefile):
        copyfile(treefile, out_newick_interim)
        copyfile(treefile, out_newick_processed)

    # Save list of MLST genes in pandas table (currently with only index with TIGR IDs)
    aligned_dir = os.path.join(automlst_interim_folder, "aligned")
    mlst_out_path_processed = os.path.join(
        automlst_processed_folder, "df_mlst_genes.csv"
    )
    mlst_out_path_interim = os.path.join(automlst_interim_folder, "df_mlst_genes.csv")

    mlst_gene_list = os.listdir(aligned_dir)

    df_mlst = pd.DataFrame(index=mlst_gene_list, columns=["aligned_path"])
    for mlst_gene in df_mlst.index:
        df_mlst.loc[mlst_gene, "aligned_path"] = os.path.join(
            aligned_dir, mlst_gene + ".faa"
        )
    df_mlst.index.name = "mlst_gene"

    df_mlst.to_csv(mlst_out_path_interim)
    df_mlst.to_csv(mlst_out_path_processed)

    return None


def get_genome_tree_table(
    automlst_processed_folder, prokka_interim_folder, gtdb_table_path, organism_info_folder
):
    """
    Copy important files from autoMLST interim to proecessed directory

    Parameters
    ----------
    1. automlst_processed_folder : str / path
        Location of the processed output directory for important autoMLST results
    2. prokka_interim_folder : str/ path
        Path to get organism information for each genome as used in prokka
    3. gtdb_table_path : str/ path
        Path to get table with gtdb information

    Returns
    -------
    1. df_genomes_tree : pd.DataFrame
        Combined datadrame of organism info ordered in phylogenetic tree (final.newick file)
    """

    newick_path = os.path.join(automlst_processed_folder, "final.newick")

    t = Phylo.read(newick_path, "newick")

    # Max distance to create better plot
    # mdist = max([t.distance(t.root, x) for x in t.get_terminals()])

    # Sort the matrix according to tip labels in the tree
    genome_ids_list = [
        x.name
        if x.name in os.listdir(prokka_interim_folder)
        else x.name[:-1] + "." + x.name[-1]
        for x in t.get_terminals()
    ]
    # The above is required in the default version of autmlst_wrapper
    # To be fixed in the orginal repo

    columns_list = [
        "genus_original",
        "species_original",
        "strain",
        "domain",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]
    df_genomes_tree = pd.DataFrame(index=genome_ids_list, columns=columns_list)
    df_genomes_tree.index.name = "genome_id"

    df_gtdb_meta = pd.read_csv(gtdb_table_path, index_col=0)

    for genome_id in df_genomes_tree.index:
        # Reading organism infor used for prokka run including strain ID
        org_info_path = os.path.join(
            organism_info_folder, genome_id, "organism_info.txt"
        )
        with open(org_info_path, "r") as file_obj:
            org_info = file_obj.readlines()[0]
            df_genomes_tree.loc[genome_id, "genus_original"] = org_info.split(",")[0]
            df_genomes_tree.loc[genome_id, "species_original"] = org_info.split(",")[1]
            df_genomes_tree.loc[genome_id, "strain"] = org_info.split(",")[2]

        # Reading gtdb metadata from JSON files

        df_genomes_tree.loc[genome_id, "domain"] = df_gtdb_meta.loc[genome_id, "Domain"]
        df_genomes_tree.loc[genome_id, "phylum"] = df_gtdb_meta.loc[genome_id, "Phylum"]
        df_genomes_tree.loc[genome_id, "class"] = df_gtdb_meta.loc[genome_id, "Class"]
        df_genomes_tree.loc[genome_id, "order"] = df_gtdb_meta.loc[genome_id, "Order"]
        df_genomes_tree.loc[genome_id, "family"] = df_gtdb_meta.loc[genome_id, "Family"]
        df_genomes_tree.loc[genome_id, "genus"] = df_gtdb_meta.loc[genome_id, "Genus"]
        df_genomes_tree.loc[genome_id, "organism"] = df_gtdb_meta.loc[genome_id, "Species"]
        try:
            df_genomes_tree.loc[genome_id, "species"] = df_gtdb_meta.loc[genome_id, "Species"].split(" ")[1]
        except IndexError:  # leave blank for empty taxonomy
            df_genomes_tree.loc[genome_id, "species"] = ""

    genomes_tree_path_interim = os.path.join(
        automlst_processed_folder, "df_genomes_tree.csv"
    )
    genomes_tree_path_processed = os.path.join(
        automlst_processed_folder, "df_genomes_tree.csv"
    )

    df_genomes_tree.to_csv(genomes_tree_path_interim)
    df_genomes_tree.to_csv(genomes_tree_path_processed)

    return df_genomes_tree


if __name__ == "__main__":
    automlst_interim_folder = sys.argv[1]
    automlst_processed_folder = sys.argv[2]
    prokka_interim_folder = sys.argv[3]
    gtdb_table_path = sys.argv[4]
    organism_info_folder = sys.argv[5]
    copy_automlst_out(automlst_interim_folder, automlst_processed_folder)
    df_genomes_tree = get_genome_tree_table(
        automlst_processed_folder, prokka_interim_folder, gtdb_table_path, organism_info_folder
    )
