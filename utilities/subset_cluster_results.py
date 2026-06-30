#!/usr/bin/env python3

import os
import argparse
import numpy as np
import pandas as pd
import scipy.sparse as sp
from Bio import SeqIO
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform

__version__ = "0.9.0b"

def load_targets(accessions_file):
    with open(accessions_file, 'r') as f:
        return set(line.strip() for line in f if line.strip())

def linkage_to_newick(linkage_matrix, labels):
    tree, _ = to_tree(linkage_matrix, rd=True)
    def build_newick(node):
        if node.is_leaf(): return f"{labels[node.id]}:{node.dist:.6f}"
        return f"({build_newick(node.left)},{build_newick(node.right)}):{node.dist:.6f}"
    return build_newick(tree) + ";"

def calculate_jaccard(A_sparse):
    """Calculates the exact Jaccard distance matrix."""
    num_genomes = A_sparse.shape[0]
    row_sums = A_sparse.sum(axis=1).A1 
    square_dist = np.zeros((num_genomes, num_genomes), dtype=np.float32)
    chunk_size = 2000
    
    for i in tqdm(range(0, num_genomes, chunk_size), desc="Calculating Jaccard batches"):
        end = min(i + chunk_size, num_genomes)
        int_chunk = A_sparse[i:end].dot(A_sparse.T).toarray()
        union_chunk = row_sums[i:end, None] + row_sums - int_chunk
        
        with np.errstate(divide='ignore', invalid='ignore'):
            jaccard_chunk = int_chunk / union_chunk
            jaccard_chunk[np.isnan(jaccard_chunk)] = 0.0
        square_dist[i:end, :] = 1.0 - jaccard_chunk

    np.fill_diagonal(square_dist, 0.0) 
    return np.clip(square_dist, 0.0, 1.0)

def main():
    parser = argparse.ArgumentParser(description="VORTEX Utility Script: Subset clustering results and metadata by a list of accessions.")
    
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}', help="Show the script's version number and exit.")
    
    parser.add_argument("--matrix_prefix", required=True, help="Prefix of the binary matrix files (e.g., results/protein_cluster_matrix_binary_seqid_0_7)")
    parser.add_argument("--clusters_dir", required=True, help="Path to the original clusters directory containing .faa files")
    parser.add_argument("--accessions", required=True, help="Text file containing genome accessions to keep (one per line)")
    parser.add_argument("--out_dir", required=True, help="Directory to save subset results")
    parser.add_argument("--taxonomy_csv", default=None, help="Path to phage_taxonomy.csv (optional, for downsampling metadata)")
    parser.add_argument("--signatures_csv", default=None, help="Path to phage_signatures.csv (optional, for downsampling metadata)")
    parser.add_argument("--fig_format", default="png", choices=["png", "pdf", "svg"], help="Format for output plots")
    parser.add_argument("--dpi", type=int, default=300, help="DPI for static plots")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    subset_clusters_dir = os.path.join(args.out_dir, "clusters")
    os.makedirs(subset_clusters_dir, exist_ok=True)

    print("Loading target accessions...")
    target_accessions = load_targets(args.accessions)

    print("Loading original sparse matrix and metadata...")
    mat_path = f"{args.matrix_prefix}_matrix.npz"
    cols_path = f"{args.matrix_prefix}_cols.txt"
    rows_path = f"{args.matrix_prefix}_rows.csv"

    with open(cols_path, "r") as f:
        all_genomes = [line.strip() for line in f if line.strip()]
    
    rows_df = pd.read_csv(rows_path)
    mat = sp.load_npz(mat_path).tocsr()

    # 1. Map target accessions to matrix column indices
    col_indices = [i for i, g in enumerate(all_genomes) if g in target_accessions]
    subset_genomes = [all_genomes[i] for i in col_indices]
    
    if not subset_genomes:
        print("Error: None of the provided accessions were found in the dataset.")
        return

    print(f"Found {len(subset_genomes)} out of {len(target_accessions)} requested genomes.")

    # 2. Subset the matrix columns
    sub_mat = mat[:, col_indices]

    # 3. Filter out rows (clusters) that have 0 occurrences in the new subset
    row_sums = sub_mat.sum(axis=1).A1
    active_row_indices = np.where(row_sums > 0)[0]
    
    final_mat = sub_mat[active_row_indices, :]
    final_rows_df = rows_df.iloc[active_row_indices].reset_index(drop=True)
    active_cluster_ids = set(final_rows_df["cluster_id"].tolist())

    print(f"Retained {len(active_cluster_ids)} clusters active in the subset.")

    # 4. Save the subset matrix as a single dense CSV
    out_csv = os.path.join(args.out_dir, "subset_matrix.csv")
    print(f"Saving subset matrix as a single dense CSV to {out_csv}...")
    
    dense_mat = final_mat.toarray()
    mat_df = pd.DataFrame(dense_mat, columns=subset_genomes)
    final_csv_df = pd.concat([final_rows_df, mat_df], axis=1)
    final_csv_df.to_csv(out_csv, index=False)

    # 5. Process Taxonomy and Signatures (If provided)
    if args.taxonomy_csv and args.signatures_csv:
        print("Downsampling taxonomy and signature mappings...")
        try:
            tax_df = pd.read_csv(args.taxonomy_csv)
            sig_df = pd.read_csv(args.signatures_csv)

            # Filter taxonomy to just our successfully mapped subset genomes
            subset_tax_df = tax_df[tax_df['Genome'].isin(subset_genomes)]
            
            # Extract all valid lineage identifiers across all ranks
            valid_lineages = set()
            if 'Family' in subset_tax_df.columns: valid_lineages.update(subset_tax_df['Family'].dropna().unique())
            if 'Subfamily' in subset_tax_df.columns: valid_lineages.update(subset_tax_df['Subfamily'].dropna().unique())
            if 'Genus' in subset_tax_df.columns: valid_lineages.update(subset_tax_df['Genus'].dropna().unique())

            # Filter signatures strictly to those lineages
            subset_sig_df = sig_df[sig_df['Lineage'].isin(valid_lineages)]

            # Save the downsampled metadata
            subset_tax_df.to_csv(os.path.join(args.out_dir, "subset_phage_taxonomy.csv"), index=False)
            subset_sig_df.to_csv(os.path.join(args.out_dir, "subset_phage_signatures.csv"), index=False)
            print(f"Saved subset taxonomy ({len(subset_tax_df)} rows) and signatures ({len(subset_sig_df)} rows).")
        except Exception as e:
            print(f"Warning: Failed to downsample metadata: {e}")

    # 6. Extract multi-FASTA files for the active clusters
    print("Extracting subset FASTA files...")
    genome_prefixes = tuple(f"{g}_" for g in subset_genomes)
    
    for cluster_id in tqdm(active_cluster_ids, desc="Parsing cluster fastas"):
        orig_fasta = os.path.join(args.clusters_dir, f"{cluster_id}.faa")
        if not os.path.exists(orig_fasta):
            continue
            
        subset_records = []
        for record in SeqIO.parse(orig_fasta, "fasta"):
            if record.id.startswith(genome_prefixes):
                subset_records.append(record)
                
        if subset_records:
            SeqIO.write(subset_records, os.path.join(subset_clusters_dir, f"{cluster_id}.faa"), "fasta")

    # 7. Calculate Jaccard distances and regenerate Visualizations
    print("Calculating Jaccard distances for subset...")
    A = final_mat.transpose().tocsr().astype(np.float32)
    square_dist = calculate_jaccard(A)
    condensed_dist = squareform(square_dist, checks=False)

    print("Generating Visualizations...")
    row_linkage = linkage(condensed_dist, method="complete")
    
    with open(os.path.join(args.out_dir, "subset_dendrogram.newick"), "w") as f: 
        f.write(linkage_to_newick(row_linkage, subset_genomes) + "\n")

    if final_mat.shape[1] <= 5000:
        sns.set_theme(font_scale=1)
        
        bin_df = pd.DataFrame(A.toarray(), index=subset_genomes)
        ax1 = sns.clustermap(bin_df, cmap="Blues", row_linkage=row_linkage, col_cluster=True, yticklabels=True, xticklabels=False, figsize=(20, 10))
        ax1.ax_col_dendrogram.set_visible(False)
        ax1.ax_heatmap.set_yticklabels(ax1.ax_heatmap.get_yticklabels(), rotation=0, fontsize=max(1.5, min(10, 100/len(subset_genomes))))
        plt.savefig(os.path.join(args.out_dir, f"subset_binary_heatmap.{args.fig_format}"), format=args.fig_format, dpi=args.dpi, bbox_inches="tight")
        plt.close("all")

        shared_pct = pd.DataFrame((1.0 - square_dist) * 100.0, index=subset_genomes, columns=subset_genomes)
        ax2 = sns.clustermap(shared_pct, cmap="Reds", row_linkage=row_linkage, col_linkage=row_linkage, yticklabels=True, xticklabels=False, figsize=(12, 10))
        ax2.ax_col_dendrogram.set_visible(False)
        ax2.ax_heatmap.set_yticklabels(ax2.ax_heatmap.get_yticklabels(), rotation=0, fontsize=max(1.5, min(10, 100/len(subset_genomes))))
        plt.savefig(os.path.join(args.out_dir, f"subset_shared_pct_heatmap.{args.fig_format}"), format=args.fig_format, dpi=args.dpi, bbox_inches="tight")
        plt.close("all")
        
        shared_pct.to_csv(os.path.join(args.out_dir, "subset_shared_pct.csv"))

    print(f"Done! Subset results saved to: {args.out_dir}")

if __name__ == "__main__":
    main()