#!/usr/bin/env python3
import os
import glob
import argparse
import scipy.sparse
import pandas as pd

def convert_sparse_to_dense(target_dir):
    # Search for any file ending in '_matrix.npz' in the target directory
    search_pattern = os.path.join(target_dir, "*_matrix.npz")
    matrix_files = glob.glob(search_pattern)

    if not matrix_files:
        print(f"No sparse matrix files (*_matrix.npz) found in: {target_dir}")
        return

    print(f"Found {len(matrix_files)} sparse matrix file(s). Processing...\n")

    for matrix_path in matrix_files:
        # Dynamically extract the prefix by stripping off the '_matrix.npz' suffix
        prefix = matrix_path.replace('_matrix.npz', '')
        
        rows_path = f"{prefix}_rows.csv"
        cols_path = f"{prefix}_cols.txt"
        out_path = f"{prefix}_combined.csv"

        # Verify the matching metadata files exist before attempting to parse
        if not (os.path.exists(rows_path) and os.path.exists(cols_path)):
            print(f"[SKIPPED] Missing matching _rows.csv or _cols.txt for: {os.path.basename(matrix_path)}")
            continue

        print(f"Converting: {os.path.basename(prefix)}...")

        # 1. Load the sparse matrix
        sparse_mat = scipy.sparse.load_npz(matrix_path)

        # 2. Load the entire row metadata (retaining all stats)
        rows_df = pd.read_csv(rows_path)

        # 3. Load the column identifiers (genomes)
        with open(cols_path, 'r') as f:
            col_names = [line.strip() for line in f.readlines() if line.strip()]

        # 4. Convert the sparse matrix to a dense pandas DataFrame
        matrix_df = pd.DataFrame(sparse_mat.toarray(), columns=col_names)

        # 5. Concatenate the statistics dataframe and the matrix dataframe side-by-side
        combined_df = pd.concat([rows_df, matrix_df], axis=1)

        # 6. Export to a single CSV file
        combined_df.to_csv(out_path, index=False)
        print(f"[SUCCESS] Saved to: {out_path}\n")

def main():
    parser = argparse.ArgumentParser(description="Convert VORTEX sparse matrices back to dense CSVs with all metadata retained.")
    
    # Add the directory argument
    parser.add_argument(
        "-d", "--dir", 
        required=True, 
        help="Path to the directory containing the VORTEX sparse matrix files (e.g., results/)"
    )
    
    args = parser.parse_args()
    
    # Ensure the provided path is a valid directory
    if not os.path.isdir(args.dir):
        print(f"Error: The directory '{args.dir}' does not exist.")
        return

    convert_sparse_to_dense(args.dir)

if __name__ == "__main__":
    main()