#!/usr/bin/env python3

import os
import argparse
import numpy as np
import pandas as pd
import scipy.sparse as sp

def load_targets(accessions_file):
    with open(accessions_file, 'r', encoding='utf-8-sig') as f:
        return set(line.strip() for line in f if line.strip())

def main():
    parser = argparse.ArgumentParser(description="VORTEX utility script: Extract a network of genomes sharing clusters with target accessions.")
    parser.add_argument("--matrix_prefix", required=True, help="Prefix of the binary matrix files (e.g., results/protein_cluster_matrix_binary_seqid_0_7)")
    parser.add_argument("--summary_csv", required=True, help="Path to the protein cluster summary CSV (e.g., results/protein_cluster_summary_seqid_0_7.csv)")
    parser.add_argument("--accessions", required=True, help="Text file containing genome accessions to query (one per line)")
    parser.add_argument("--out_dir", required=True, help="Directory to save the resulting CSVs")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    print("Loading target accessions...")
    target_accessions = load_targets(args.accessions)

    print("Loading original sparse matrix and metadata...")
    mat_path = f"{args.matrix_prefix}_matrix.npz"
    cols_path = f"{args.matrix_prefix}_cols.txt"
    rows_path = f"{args.matrix_prefix}_rows.csv"

    with open(cols_path, "r", encoding='utf-8-sig') as f:
        all_genomes = np.array([line.strip() for line in f if line.strip()])
    
    rows_df = pd.read_csv(rows_path)
    mat = sp.load_npz(mat_path).tocsr()

    # 1. Map target accessions to matrix column indices
    target_indices = [i for i, g in enumerate(all_genomes) if g in target_accessions]
    
    if not target_indices:
        print("Error: None of the provided accessions were found in the dataset.")
        return

    print(f"Found {len(target_indices)} out of {len(target_accessions)} requested target genomes.")

    # DEBUG: Find the exact strings that are in your file but NOT in the matrix
    missing_genomes = target_accessions - set(all_genomes)
    print(f"\nDEBUG: The following genomes failed the exact match:\n{missing_genomes}\n")

    # 2. Identify protein clusters where a protein is present in the target genomes
    # Slice matrix to target genomes, sum rows to find which clusters they possess
    target_mat = mat[:, target_indices]
    active_cluster_mask = target_mat.sum(axis=1).A1 > 0
    active_cluster_indices = np.where(active_cluster_mask)[0]
    
    active_cluster_ids = rows_df.iloc[active_cluster_indices]["cluster_id"].values
    print(f"Identified {len(active_cluster_indices)} unique protein clusters present in the target genomes.")

    # 3. Identify all genomes that have a representative protein within the clusters
    # Slice original matrix to only the active clusters, sum columns to find all genomes that have them
    clusters_mat = mat[active_cluster_indices, :]
    connected_genome_mask = clusters_mat.sum(axis=0).A1 > 0
    connected_genome_indices = np.where(connected_genome_mask)[0]
    
    connected_genomes = all_genomes[connected_genome_indices]
    print(f"Expanded network: Found {len(connected_genomes)} total genomes sharing at least one of these clusters.")

    # 4. Generate Output 1: Shared Network Matrix CSV
    out_matrix_csv = os.path.join(args.out_dir, "shared_network_matrix.csv")
    print(f"Saving shared network matrix to {out_matrix_csv}...")
    
    # Extract the exact cross-section of Active Clusters x Connected Genomes
    final_mat = mat[np.ix_(active_cluster_indices, connected_genome_indices)]
    
    # Convert to dense DataFrame and combine with row metadata
    mat_df = pd.DataFrame(final_mat.toarray(), columns=connected_genomes)
    final_rows_df = rows_df.iloc[active_cluster_indices].reset_index(drop=True)
    
    final_csv_df = pd.concat([final_rows_df, mat_df], axis=1)
    final_csv_df.to_csv(out_matrix_csv, index=False)

    # 5. Generate Output 2: Protein IDs (genome_cds_number) mapping
    out_proteins_csv = os.path.join(args.out_dir, "shared_network_proteins.csv")
    print(f"Extracting specific protein CDS numbers to {out_proteins_csv}...")
    
    summary_df = pd.read_csv(args.summary_csv)
    
    # Filter the summary to only include the clusters we care about
    filtered_summary = summary_df[summary_df["cluster_id"].isin(active_cluster_ids)].copy()
    
    # The protein IDs are formatted as GenomeID_CDSNumber. Extract the GenomeID to filter properly.
    # We use rsplit with n=1 to handle genome accessions that might naturally contain underscores
    filtered_summary["genome_id"] = filtered_summary["protein"].str.rsplit("_", n=1).str[0]
    
    # Keep only the proteins belonging to the connected genomes
    connected_set = set(connected_genomes)
    final_proteins_df = filtered_summary[filtered_summary["genome_id"].isin(connected_set)]
    
    # Clean up and save
    final_proteins_df = final_proteins_df.drop(columns=["genome_id"])
    final_proteins_df.to_csv(out_proteins_csv, index=False)

    print("Network extraction complete!")

if __name__ == "__main__":
    main()