#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys


def calculate_shared_proteins(matrix):
    num_genomes = matrix.shape[1]
    shared_proteins = np.zeros((num_genomes, num_genomes), dtype=int)

    for i in range(num_genomes):
        for j in range(i, num_genomes):
            intersection = np.sum(np.logical_and(matrix[:, i], matrix[:, j]))
            shared_proteins[i, j] = intersection
            shared_proteins[j, i] = intersection  # Symmetric matrix

    return shared_proteins


def calculate_percentage_shared_proteins(shared_proteins, binary_data):
    num_genomes = shared_proteins.shape[0]
    percentages = np.zeros((num_genomes, num_genomes))

    for i in range(num_genomes):
        for j in range(i, num_genomes):
            total_proteins = np.sum(np.logical_or(binary_data[:, i], binary_data[:, j]))
            if total_proteins > 0:
                percentage = (shared_proteins[i, j] / total_proteins) * 100
                percentages[i, j] = percentage
                percentages[j, i] = percentage  # Symmetric matrix

    return percentages


def process_taxonomic_level(level_name, full_df, metadata, args):
    """
    Helper function to process a specific taxonomic level (family, genus, or subfamily).
    Returns a DataFrame ready for plotting, or None if no valid data found.
    """
    print(f"\n--- Processing level: {level_name} ---")

    # Identify valid groups (more than 1 member)
    counts = metadata[level_name].value_counts()
    valid_groups = counts[counts > 1].index

    if len(valid_groups) == 0:
        print(f"No {level_name} groups found with >1 member. Skipping.")
        return None

    col_1 = f"{level_name}_1"
    col_2 = f"{level_name}_2"

    comparisons_to_concat = []

    # Process Intra-group similarities (Group A vs Group A)
    intra_sim = full_df[
        (full_df[col_1] == full_df[col_2]) & (full_df[col_1].isin(valid_groups))
    ].copy()

    if not intra_sim.empty:
        intra_sim["ComparisonType"] = intra_sim[col_1] + " (intra)"
        comparisons_to_concat.append(intra_sim)

    # Process Inter-group similarities (Group A vs Group B)
    inter_sim = full_df[
        (full_df[col_1] != full_df[col_2])
        & (full_df[col_1].isin(valid_groups))
        & (full_df[col_2].isin(valid_groups))
    ].copy()

    if not inter_sim.empty:
        # Use min() and max() to sort the names alphabetically. This ensures 'B' and 'A' always becomes 'A vs B', preventing duplicates.
        inter_sim["ComparisonType"] = inter_sim.apply(
            lambda row: f"{min(row[col_1], row[col_2])} vs {max(row[col_1], row[col_2])}",
            axis=1
        )
        comparisons_to_concat.append(inter_sim)

    # Combine and return
    if not comparisons_to_concat:
        print(f"No comparisons generated for {level_name}.")
        return None

    boxplot_data = pd.concat(comparisons_to_concat)

    # Optional: Save intermediate data
    if args.save_intermediate:
        filename = f"all_{level_name}_comparisons.csv"
        boxplot_data.to_csv(filename, index=False)
        print(f"Saved intermediate comparisons: {filename}")

    return boxplot_data


def plot_taxonomic_level(boxplot_data, level_name, args):
    """
    Helper function to plot the boxplots for a specific level.
    """
    if boxplot_data is None or boxplot_data.empty:
        return

    # Count unique comparisons to determine plot height
    num_comparisons = boxplot_data["ComparisonType"].nunique()

    # Base height for margins/titles + per-row height
    base_height = 2.0
    row_height = 0.4  # Adjust this value to increase/decrease spacing per bar

    # Calculate total height, enforcing a minimum of 8 inches
    calculated_height = base_height + (num_comparisons * row_height)
    fig_height = max(8.0, calculated_height)

    plt.figure(figsize=(20, fig_height))

    sns.boxplot(
        data=boxplot_data,
        x="Similarity",
        y="ComparisonType",
        hue="ComparisonType",
        palette="vlag",
        legend=False,
    )

    plt.title(f"Intra- and Inter-{level_name} Similarity (% Shared PCs)")
    plt.xlabel("Percentage of Shared PCs")
    plt.ylabel("Comparison")
    plt.grid(True)
    plt.tight_layout()

    outfile = f"{args.prefix}{level_name}_similarity.svg"
    plt.savefig(outfile, format="svg")
    print(f"Saved plot: {outfile}")

def main():
    # --- Argument Parser ---
    parser = argparse.ArgumentParser(
        description="Calculate and plot intra/inter-group protein cluster similarities (Family/Genus/Subfamily)."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to input CSV file containing the binary protein cluster matrix output from VORTEX.",
    )
    parser.add_argument(
        "-m",
        "--metadata",
        required=True,
        help="Path to metadata CSV file. Looks for 'genome', 'family', 'genus', 'subfamily' columns.",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        default="Taxon_cluster_analysis_",
        help="Prefix for output SVG files (default: Taxon_cluster_analysis_).",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="percentage_shared_proteins.csv",
        help="Optional path to save the resulting percentage similarity matrix.",
    )
    parser.add_argument(
        "--save-intermediate",
        action="store_true",
        help="If set, save intermediate comparison CSVs.",
    )

    args = parser.parse_args()

    # Load and process input matrix
    try:
        data = pd.read_csv(args.input, index_col=0)
    except FileNotFoundError:
        sys.exit(f"Error: Input file '{args.input}' not found.")

    # Filter out non-genome metadata columns in the VORTEX output
    columns_to_ignore = [
        col for col in data.columns 
        if col.startswith(('description','num_unique_genomes','max_paralogs_per_genome','pident_min','pident_max','pident_mean','pident_std','length_min','length_max','length_mean','length_std'))
    ]
    if columns_to_ignore:
        print(f"Ignoring extra metric columns: {columns_to_ignore}")
        data = data.drop(columns=columns_to_ignore)
    
    genome_names = data.columns
    binary_data = (data > 0).astype(int).values

    print("Calculating shared proteins...")
    shared_proteins = calculate_shared_proteins(binary_data)
    percentage_shared_proteins = calculate_percentage_shared_proteins(
        shared_proteins, binary_data
    )

    percentage_shared_proteins_df = pd.DataFrame(
        percentage_shared_proteins, index=genome_names, columns=genome_names
    )

    print("\n=== Percentage of Shared Proteins Matrix (Head) ===")
    print(percentage_shared_proteins_df.head())

    if args.output:
        percentage_shared_proteins_df.to_csv(args.output, index=True)
        print(f"\nSaved percentage shared proteins matrix to: {args.output}")

    # Metadata processing...
    try:
        metadata = pd.read_csv(args.metadata)
        metadata.columns = metadata.columns.str.lower()
    except FileNotFoundError:
        sys.exit(f"Error: Metadata file '{args.metadata}' not found.")

    if "genome" not in metadata.columns:
        sys.exit("Error: Metadata must contain a 'genome' column.")

    metadata.set_index("genome", inplace=True)

    # transform similarity matrix to long format
    percentage_similarity_df = percentage_shared_proteins_df.stack().reset_index()
    percentage_similarity_df.columns = ["Genome1", "Genome2", "Similarity"]
    # Filter for unique pairs (triangle)
    percentage_similarity_df = percentage_similarity_df[
        percentage_similarity_df["Genome1"] < percentage_similarity_df["Genome2"]
    ]

    potential_levels = ["family", "genus", "subfamily"]
    available_levels = [col for col in potential_levels if col in metadata.columns]

    if not available_levels:
        print(
            "\nWARNING: None of the expected columns (family, genus, subfamily) were found in metadata."
        )
        print(f"Available columns: {list(metadata.columns)}")
        sys.exit("Cannot proceed with grouping analysis.")
    else:
        print(f"\nFound taxonomic levels to analyze: {available_levels}")

    # Merge available metadata
    # Merge once for Genome1 and once for Genome2
    percentage_similarity_df = percentage_similarity_df.merge(
        metadata[available_levels], left_on="Genome1", right_index=True
    )
    percentage_similarity_df = percentage_similarity_df.merge(
        metadata[available_levels],
        left_on="Genome2",
        right_index=True,
        suffixes=("_1", "_2"),
    )

    # Iterative Processing
    # Loop through all found levels (family, genus, etc.) and process them
    sns.set_theme(style="whitegrid")

    for level in available_levels:
        boxplot_data = process_taxonomic_level(
            level, percentage_similarity_df, metadata, args
        )
        if boxplot_data is not None:
            plot_taxonomic_level(boxplot_data, level, args)

    print("\nDone.")


if __name__ == "__main__":
    main()