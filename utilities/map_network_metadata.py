#!/usr/bin/env python3
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Map new metadata to an existing VORTEX network and filter by edge weight.")
    parser.add_argument("-n", "--network", required=True, help="Path to the existing phage_global_network_edges.csv")
    parser.add_argument("-m", "--metadata", required=True, help="Path to the new metadata CSV")
    parser.add_argument("-o", "--output", required=True, help="Path to save the updated network edges CSV")
    parser.add_argument("-g", "--genome_col", default="Genome", help="Column name in the metadata file containing genome accessions")
    parser.add_argument("-c", "--clean", action="store_true", help="Remove existing metadata from the network before mapping")
    parser.add_argument("-w", "--min_weight", type=float, default=0.0, help="Minimum jaccard threshold to retain an edge, e.g. -w 0.05: remove all edges <= 0.05")
    
    args = parser.parse_args()

    # Load the existing edge list
    edges_df = pd.read_csv(args.network)
    
    # Apply the weight threshold filter to filter connections below distance x
    if args.min_weight > 0.0:
        initial_count = len(edges_df)
        edges_df = edges_df[edges_df["Weight"] >= args.min_weight]
        print(f"Filtered edges from {initial_count} to {len(edges_df)} using weight threshold >= {args.min_weight}")

    # Load the new metadata
    meta_df = pd.read_csv(args.metadata, dtype=str)

    if args.genome_col not in meta_df.columns:
        print(f"Error: '{args.genome_col}' not found in metadata.")
        return

    # Strip old metadata if requested, keeping only columns pertaining to jaccard distance and shared protein clusters
    if args.clean:
        core_cols = ["Source", "Target", "Weight", "Shared_Clusters"]
        edges_df = edges_df[[c for c in core_cols if c in edges_df.columns]]

    meta_df.set_index(args.genome_col, inplace=True)
    meta_dict = meta_df.to_dict(orient="index")

    # Apply mapping for both source and target nodes
    for col in meta_df.columns:
        edges_df[f"Source_{col}"] = edges_df["Source"].map(lambda x: meta_dict.get(x, {}).get(col, "Unknown"))
        edges_df[f"Target_{col}"] = edges_df["Target"].map(lambda x: meta_dict.get(x, {}).get(col, "Unknown"))

    edges_df.to_csv(args.output, index=False)
    print(f"Saved updated network to {args.output}")

if __name__ == "__main__":
    main()