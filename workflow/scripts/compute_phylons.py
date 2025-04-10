"""
This script computes phylons automatically by following the process in the pyphylon examples:
https://github.com/SBRG/pyphylon/tree/c0f732b170593719272d94edcc32f9d33a87ba30/examples

And more broadly the protocol described in the E. coli phylons paper: https://doi.org/10.1128/msphere.00532-24

1. Run clustering on mash correlation distances. The amount of large enough clusters is used to inform the NMF rank
2. Run MCA to further inform the NMF rank
3. Run NMF with the ranks informed by steps 1 and 2
4. Compute metrics of reconstruction accuracy
5. Run NMF with additional ranks around the best rank according to step 4
6. Save the L and A phylon matrices for the best NMF rank

For a broad explanation on the concept of phylons, see the 2 links above.
"""
import argparse
import logging
from math import ceil
from os import makedirs

import pandas as pd
import prince
import scipy
import scipy.sparse
from pyphylon.mash import cluster_corr_dist, sensitivity_analysis
from pyphylon.models import (binarize_nmf_outputs,
                             calculate_nmf_reconstruction_metrics,
                             generate_nmf_reconstructions,
                             normalize_nmf_outputs, run_nmf)

if __name__ == "__main__":
    # Percentile of total genomes that is considered to be a "big enough" cluster of sequences (as per the mash 
    # clustering that informs NMF rank)
    PERCENTILE_CLUSTERS_THRESHOLD = 0.01

    # Maximum rank to run NMF with (the algorithm seems to be O(n²) w.r.t. rank and the number of genomes can get 
    # pretty big on popular species)
    MAX_NMF_RANK = 75

    parser = argparse.ArgumentParser(description="Compute phylons")
    parser.add_argument("--mash_distances", type=str, required=True, help="Path to the mash distances CSV")
    parser.add_argument("--gene_presence", type=str, required=True, help="Path of gene presence/absence CSV")
    parser.add_argument("--pangenome_summary", type=str, required=True, help="Path of the pangenome summary CSV")
    parser.add_argument("--output_dir", type=str, required=True, help="Path of the output directory")

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(asctime)s - %(levelname)s] - %(message)s")
    logger = logging.getLogger()
    logger.info("Starting the computation of phylons")
    logger.info(f"Args: {args}")

    mash_distances_df = pd.read_csv(args.mash_distances, index_col=0)
    gene_presence_df = pd.read_csv(args.gene_presence, index_col=0)
    pangenome_summary_df = pd.read_csv(args.pangenome_summary, index_col=0)

    n_genomes = mash_distances_df.shape[0]

    ## mash clustering ##
    logger.info("Clustering on mash distances")
    mash_corr_df = mash_distances_df.corr()
    mash_corr_distances_df = 1 - mash_corr_df

    thresh_n_clusters_df, elbow_idx, elbow_threshold = sensitivity_analysis(mash_corr_distances_df)
    n_clusters_elbow = thresh_n_clusters_df["num_clusters"][elbow_idx]
    logger.info(f"Num clusters decelerates at {n_clusters_elbow} clusters (threshold: {elbow_threshold:.5f})")

    link, dist, mash_clusters_df = cluster_corr_dist(mash_corr_distances_df, thresh=elbow_threshold)

    small_cluster_threshold = ceil(n_genomes * PERCENTILE_CLUSTERS_THRESHOLD)
    logger.info(
        f"Threshold for small clusters is {small_cluster_threshold} genomes ({n_genomes=} | {PERCENTILE_CLUSTERS_THRESHOLD=})"
    )

    cluster_sizes_df = pd.DataFrame(mash_clusters_df["cluster"].value_counts()).reset_index()
    n_big_enough_mash_clusters = (cluster_sizes_df["count"] >= small_cluster_threshold).sum()

    ## MCA ##
    gene_presence_sparse_df = gene_presence_df.astype(pd.SparseDtype("int8", 0))
    gene_presence_coo: scipy.sparse.coo_matrix = gene_presence_sparse_df.sparse.to_coo()
    gene_presence_csr = scipy.sparse.csr_matrix(gene_presence_coo)

    accessory_genes_mask = pangenome_summary_df["pangenome_class_2"] == "Accessory"
    accessory_genes = pangenome_summary_df.index[accessory_genes_mask]

    accessory_gene_presence_df = gene_presence_df.loc[accessory_genes]
    logger.info(f"{len(accessory_genes)} accessory genes")

    mca = prince.MCA(
        n_components=accessory_gene_presence_df.shape[1],
        n_iter=3,
        copy=True,
        check_input=True,
        engine="sklearn",
        random_state=42,
    )

    mca = mca.fit(accessory_gene_presence_df)
    explained_variance_percentage = mca.percentage_of_variance_
    cumulative_variance = pd.Series(explained_variance_percentage).cumsum()
    cumulative_variance: dict[int, float] = {
        num: cumulative_variance[cumulative_variance >= num].index[0] for num in range(1, 99)
    }

    rank_list = [
        2,
        5,
        7,
        10,
        n_big_enough_mash_clusters,
        cumulative_variance[70],
        cumulative_variance[75],
        cumulative_variance[80],
        cumulative_variance[85],
        cumulative_variance[90],
    ]
    rank_list = sorted(set(min(n, MAX_NMF_RANK) for n in rank_list))

    ## NMF ##
    logger.info(f"Running NMF with ranks: {rank_list}")
    W_dict, H_dict = run_nmf(data=accessory_gene_presence_df, ranks=rank_list, max_iter=40_000)

    L_norm_dict, A_norm_dict = normalize_nmf_outputs(accessory_gene_presence_df, W_dict, H_dict)
    L_binarized_dict, A_binarized_dict = binarize_nmf_outputs(L_norm_dict, A_norm_dict)
    P_reconstructed_dict, P_error_dict, P_confusion_dict = generate_nmf_reconstructions(
        accessory_gene_presence_df, L_binarized_dict, A_binarized_dict
    )
    reconstruction_metrics_df = calculate_nmf_reconstruction_metrics(P_reconstructed_dict, P_confusion_dict)
    reconstruction_metrics_df.sort_values(by="AIC")

    best_rank = reconstruction_metrics_df["AIC"].idxmin(axis=0)
    
    # re-run with additional ranks around the best rank
    extra_ranks = list(set(max(1, i) for i in range(best_rank - 4, best_rank + 5)))
    logger.info(f"Best rank: {best_rank}. Running extra ranks: {extra_ranks}")

    W_dict, H_dict = run_nmf(data=accessory_gene_presence_df, ranks=extra_ranks, max_iter=40_000)

    L_norm_dict, A_norm_dict = normalize_nmf_outputs(accessory_gene_presence_df, W_dict, H_dict)
    L_binarized_dict, A_binarized_dict = binarize_nmf_outputs(L_norm_dict, A_norm_dict)
    P_reconstructed_dict, P_error_dict, P_confusion_dict = generate_nmf_reconstructions(
        accessory_gene_presence_df, L_binarized_dict, A_binarized_dict
    )
    reconstruction_metrics_extra_df = calculate_nmf_reconstruction_metrics(P_reconstructed_dict, P_confusion_dict)

    reconstruction_metrics_all_df = pd.concat([reconstruction_metrics_df, reconstruction_metrics_extra_df]).sort_values(by='AIC')
    
    best_rank = int(reconstruction_metrics_all_df.reset_index().loc[0, 'rank'])
    logger.info(f"Best rank after extra ranks: {best_rank}")

    L_best = L_norm_dict[best_rank]
    L_binarized_best = L_binarized_dict[best_rank]
    A_best = A_norm_dict[best_rank]
    A_binarized_best = A_binarized_dict[best_rank]

    logger.info("Saving L, L(binarized), A and A (binarized) matrices for best rank")
    makedirs(args.output_dir, exist_ok=True)
    L_best.to_csv(f'{args.output_dir}/NMF_L.csv')
    L_binarized_best.to_csv(f'{args.output_dir}/NMF_L_binarized.csv')
    A_best.to_csv(f'{args.output_dir}/NMF_A.csv')
    A_binarized_best.to_csv(f'{args.output_dir}/NMF_A_binarized.csv')

    logger.info("Done")
    