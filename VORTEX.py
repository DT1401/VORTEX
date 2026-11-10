#!/usr/bin/env python3

import os
import glob
import subprocess
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
import logging
import shutil
import sys
import collections
import uuid
import tempfile
import re
import scipy.sparse as sp
from itertools import combinations
import hashlib
import warnings

# Version information
__version__ = "0.9.3b"

TQDM_DISABLE = not sys.stdout.isatty()

# DuckDB for Out-Of-Core Parquet Processing
try:
    import duckdb
    HAS_DUCKDB = True
except Exception:
    HAS_DUCKDB = False

# Visualization imports (Heatmaps & Dendrograms)
try:
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import linkage, to_tree, fcluster, leaves_list
    from scipy.spatial.distance import squareform
    HAS_CORE_VIS = True
except Exception:
    HAS_CORE_VIS = False

# Sparse Math, Networks & Scalable Visualization imports
try:
    import scipy.sparse as sp
    from scipy.sparse import csgraph
    from sklearn.decomposition import TruncatedSVD
    from sklearn.metrics import silhouette_score
    import plotly.express as px
    HAS_SCALABLE_VIS = True
except Exception:
    HAS_SCALABLE_VIS = False

try:
    import umap
    HAS_UMAP = True
except Exception:
    HAS_UMAP = False

# UpSet Plot imports
try:
    from upsetplot import UpSet, from_memberships
    HAS_UPSET = True
except Exception:
    HAS_UPSET = False

# Networkit imports
try:
    import networkit as nk
    HAS_NETWORKIT = True
except ImportError:
    HAS_NETWORKIT = False

# Logging
def setup_logger(log_file=None, console_level=logging.INFO):
    logger = logging.getLogger("protein_pipeline")
    logger.setLevel(logging.DEBUG)
    if logger.handlers:
        for h in list(logger.handlers): logger.removeHandler(h)

    ch = logging.StreamHandler()
    ch.setLevel(console_level)
    ch_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", "%H:%M:%S")
    ch.setFormatter(ch_formatter)
    logger.addHandler(ch)

    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setLevel(logging.DEBUG)
        fh_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(name)s - %(message)s")
        fh.setFormatter(fh_formatter)
        logger.addHandler(fh)

    return logger

def filter_duplicate_fasta(input_file, output_file, logger=None):
    """
    Reads a FASTA file and writes it to an output file, skipping any 
    duplicate IDs to ensure downstream dictionary operations and MMseqs2 analysis 
    do not silently drop or merge sequences.
    """
    seen_ids = set()
    skipped_count = 0
    
    with open(output_file, "w") as out_f:
        for record in SeqIO.parse(input_file, "fasta"):
            if record.id in seen_ids:
                skipped_count += 1
                continue # Skip writing this record
            
            seen_ids.add(record.id)
            SeqIO.write(record, out_f, "fasta")
            
    if skipped_count > 0 and logger:
        logger.warning(f"Skipped {skipped_count} duplicate FASTA records to maintain ID uniqueness.")

def extract_proteins(genbank_file, out_handle, protein_counter, logger=None):
    try:
        base_name = os.path.splitext(os.path.basename(genbank_file))[0]
        safe_genome_id = re.sub(r'[^A-Za-z0-9_.]', '_', base_name)
        with open(genbank_file, "r") as gb_handle:
            for record in SeqIO.parse(gb_handle, "genbank"):
                # 1. Define the genome ID using the RECORD ID, not the filename
                # (Strip problematic characters just in case)
                safe_genome_id = str(record.id).replace(".", "_").replace(" ", "_")
                
                # 2. Reset the counter to 1 for every new record
                protein_counter = 1 
                
                for feature in record.features:
                    if feature.type == "CDS" and "translation" in feature.qualifiers:
                        protein_sequence = feature.qualifiers["translation"][0]
                        product = feature.qualifiers.get("product", ["Unknown"])[0]
                        protein_acc = feature.qualifiers.get("protein_id", ["Unknown"])[0]
                        
                        protein_id = f"{safe_genome_id}_{protein_counter}"
                        out_handle.write(f">{protein_id} {product} {protein_acc}\n{protein_sequence}\n")
                        
                        protein_counter += 1

    except Exception as e:
        if logger: logger.error(f"Error extracting proteins from {genbank_file}: {e}")
        raise
    return protein_counter


def run_mmseqs(input_file, cluster_tsv_file, log, threads, config, db_path, cluster_db_path, tmp_dir, logger=None):
    if shutil.which("mmseqs") is None:
        raise FileNotFoundError("MMseqs2 binary not found. Please install MMseqs2 and ensure it is in PATH.")

    params = {
        "db": "proteome", "db_path": db_path, "cluster_db": "clusters", "cluster_db_path": cluster_db_path,
        "min_seq_id": config["mmseqs2"].get("min_seq_id", 0.7), "cov_mode": config["mmseqs2"].get("cov_mode", 0),
        "sensitivity": config["mmseqs2"].get("sensitivity", 7.5), "coverage": config["mmseqs2"].get("coverage", 0.5),
    }

    os.makedirs(tmp_dir, exist_ok=True)
    os.makedirs(params["db_path"], exist_ok=True)
    os.makedirs(params["cluster_db_path"], exist_ok=True)

    createdb_cmd = f"mmseqs createdb {input_file} {params['db_path']}/{params['db']}"
    cluster_cmd = f"mmseqs cluster {params['db_path']}/{params['db']} {params['cluster_db_path']}/{params['cluster_db']} {tmp_dir} --min-seq-id {params['min_seq_id']} --cov-mode {params['cov_mode']} --threads {threads} -s {params['sensitivity']} -c {params['coverage']} --cluster-mode 1"
    align_db = f"{params['cluster_db_path']}/{params['cluster_db']}_aligned"
    align_cmd = f"mmseqs align {params['db_path']}/{params['db']} {params['db_path']}/{params['db']} {params['cluster_db_path']}/{params['cluster_db']} {align_db} --threads {threads} -a -c {params['coverage']} --min-seq-id {params['min_seq_id']}"
    createtsv_cmd = f"mmseqs createtsv {params['db_path']}/{params['db']} {params['db_path']}/{params['db']} {params['cluster_db_path']}/{params['cluster_db']} {cluster_tsv_file} --threads {threads}"

    try:
        if not os.path.exists(f"{params['db_path']}/{params['db']}.dbtype"):
            if logger: logger.info("Running: mmseqs createdb")
            with open(log["createdb"], "w") as fh: subprocess.run(createdb_cmd, shell=True, check=True, stdout=fh, stderr=subprocess.STDOUT)

        if not os.path.exists(f"{params['cluster_db_path']}/{params['cluster_db']}.dbtype"):
            if logger: logger.info("Running: mmseqs cluster")
            with open(log["cluster"], "w") as fh: subprocess.run(cluster_cmd, shell=True, check=True, stdout=fh, stderr=subprocess.STDOUT)

        if not os.path.exists(f"{align_db}.dbtype"):
            if logger: logger.info("Running: mmseqs align")
            with open(log.get("align", os.devnull), "a") as fh: subprocess.run(align_cmd, shell=True, check=True, stdout=fh, stderr=subprocess.STDOUT)

        if not (os.path.exists(cluster_tsv_file) and os.path.getsize(cluster_tsv_file) > 0):
            if logger: logger.info("Running: mmseqs createtsv")
            with open(log.get("createtsv", os.devnull), "w") as fh: subprocess.run(createtsv_cmd, shell=True, check=True, stdout=fh, stderr=subprocess.STDOUT)

    except subprocess.CalledProcessError as e:
        if logger: logger.error(f"MMseqs2 failed: {e}")
        raise

    try:
        if os.path.exists(tmp_dir): shutil.rmtree(tmp_dir)
    except Exception: pass


def generate_clusters(cluster_tsv_file, output_fasta, clusters_dir, logger=None):
    os.makedirs(clusters_dir, exist_ok=True)
    with open(cluster_tsv_file, "r") as tsv_file:
        clusters = {}
        for line in tsv_file:
            parts = line.strip().split("\t")
            if len(parts) >= 2: clusters.setdefault(parts[0], []).append(parts[1])

    records = {record.id: record for record in SeqIO.parse(output_fasta, "fasta")}
    cluster_mapping, cluster_length_stats = {}, {}

    for i, (repr_id, prots) in enumerate(tqdm(list(clusters.items()), desc="Generating cluster FASTA", disable=TQDM_DISABLE), 1):
        cluster_id = f"cluster_{i}"
        cluster_mapping[repr_id] = cluster_id
        cluster_fasta = os.path.join(clusters_dir, f"{cluster_id}.faa")
        
        seqs = [records[prot] for prot in prots if prot in records]
        if seqs:
            SeqIO.write(seqs, cluster_fasta, "fasta")
            lengths = [len(s.seq) for s in seqs]
            cluster_length_stats[cluster_id] = {
                "length_min": np.min(lengths), "length_max": np.max(lengths),
                "length_mean": np.mean(lengths), "length_std": np.std(lengths) if len(lengths) > 1 else 0.0,
            }

    if logger: logger.info(f"Cluster FASTA files saved to {clusters_dir}")
    return cluster_mapping, cluster_length_stats


def process_mmseqs_parquet_partitioned(
    query_db, target_db, align_db, threads, duckdb_mem, cluster_mapping, cluster_length_stats, parquet_dir, logger=None
):
    if not HAS_DUCKDB: raise ImportError("DuckDB is required for large-scale Parquet processing.")

    con = duckdb.connect(database=':memory:')
    con.execute(f"PRAGMA threads={threads}")
    con.execute(f"PRAGMA memory_limit='{duckdb_mem}'")
    
    # Disable insertion order to prevent massive memory buffering during PARTITION_BY
    con.execute("SET preserve_insertion_order=false;")
    
    # et a temporary directory so DuckDB can safely spill to disk instead of crashing
    spill_dir = os.path.join(tempfile.gettempdir(), "duckdb_spill")
    os.makedirs(spill_dir, exist_ok=True)
    con.execute(f"SET temp_directory='{spill_dir}';")
    
    mapping_df = pd.DataFrame(list(cluster_mapping.items()), columns=['query', 'cluster_name'])
    con.register('mapping_df', mapping_df)

    if glob.glob(os.path.join(parquet_dir, "**", "*.parquet"), recursive=True):
        if logger: logger.info("Found existing Parquet partitions. Skipping stream creation.")
    else:
        if logger: logger.info("Writing MMseqs alignments to temporary disk file...")
        os.makedirs(parquet_dir, exist_ok=True)
        temp_tsv_path = tempfile.NamedTemporaryFile(delete=False, suffix=".tsv").name

        convertalis_cmd = f"mmseqs convertalis {query_db} {target_db} {align_db} {temp_tsv_path} --format-output query,target,bits,qcov,tcov,pident --threads {threads} -v 1"

        try:
            subprocess.run(convertalis_cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            if logger: logger.info("Partitioning TSV into Parquet via DuckDB...")
            
            con.execute(f"""
                COPY (
                    SELECT COALESCE(m.cluster_name, t.query) AS cluster_name,
                           COALESCE(NULLIF(regexp_extract(t.target, '^(.*)_[0-9]+$', 1), ''), t.target) AS genome,
                           t.bits, t.pident, (t.bits * least(t.qcov, t.tcov)) AS weighted_bits
                    FROM read_csv_auto('{temp_tsv_path}', delim='\t', header=False, names=['query', 'target', 'bits', 'qcov', 'tcov', 'pident']) AS t
                    LEFT JOIN mapping_df m ON t.query = m.query
                    ORDER BY genome
                ) TO '{parquet_dir}' (FORMAT PARQUET, PARTITION_BY genome);
            """)
            
        except subprocess.CalledProcessError as e:
            if logger: logger.error(f"MMseqs convertalis failed: {e.stderr.decode('utf-8') if e.stderr else 'Unknown'}")
            raise
        finally:
            if os.path.exists(temp_tsv_path): os.remove(temp_tsv_path)

    if logger: logger.info("Calculating cluster metrics from Parquet...")
    robustness_df = con.execute(f"""
        WITH mmseqs_data AS (SELECT * FROM read_parquet('{os.path.join(parquet_dir, "**", "*.parquet").replace(chr(92), "/")}'))
        SELECT s.cluster_name as cluster_id, count(DISTINCT s.genome) as num_unique_genomes,
               max(p.paralogs) as max_paralogs_per_genome,
               ROUND(min(s.pident), 2) AS pident_min, ROUND(max(s.pident), 2) AS pident_max,
               ROUND(avg(s.pident), 2) AS pident_mean, ROUND(COALESCE(stddev_samp(s.pident), 0.0), 2) AS pident_std
        FROM mmseqs_data s
        JOIN (SELECT cluster_name, genome, count(*) as paralogs FROM mmseqs_data GROUP BY cluster_name, genome) p
        ON s.cluster_name = p.cluster_name AND s.genome = p.genome
        GROUP BY s.cluster_name
    """).df()
    con.close()

    length_df = pd.DataFrame.from_dict(cluster_length_stats, orient="index").reset_index(names="cluster_id").round(2)
    final_robustness_df = pd.merge(robustness_df, length_df, on="cluster_id", how="left").fillna(0)

    if logger: logger.info(f"Bitscore matrix retained natively in Parquet format at: {parquet_dir}")
    return final_robustness_df


def build_summary_and_records(fasta_dir, output_csv, output_fasta, logger=None):
    prot_to_desc = {}
    for record in SeqIO.parse(output_fasta, "fasta"):
        if record.id not in prot_to_desc:
            # Strip the leading protein_id
            raw_desc = record.description.split(record.id, maxsplit=1)[1].strip() if record.id in record.description else record.description
            
            clean_desc = raw_desc
            if " " in raw_desc:
                base_desc, last_word = raw_desc.rsplit(" ", 1)
        
                # Drop the last word if it is unknown or 
                if last_word.lower() == "unknown" or any(c.isdigit() for c in last_word) or "_" in last_word:
                    clean_desc = base_desc
            
            # Save it back to the dictionary
            prot_to_desc[record.id] = clean_desc

    all_genomes, cluster_records, desc_dict, summary_csv_rows = set(), [], {}, []

    if logger: logger.info("Scanning clusters to build single-pass presence maps...")
    for cluster_file in tqdm(sorted([f for f in os.listdir(fasta_dir) if f.endswith(".faa")]), desc="Parsing clusters", disable=TQDM_DISABLE):
        cluster_name = cluster_file.replace(".faa", "")
        genome_ids, seen_prots, desc_counter = [], set(), collections.Counter()
        
        for record in SeqIO.parse(os.path.join(fasta_dir, cluster_file), "fasta"):
            if record.id in seen_prots: continue
            seen_prots.add(record.id)
            desc = prot_to_desc.get(record.id, "Unknown")
            desc_counter[desc] += 1
            summary_csv_rows.append([record.id, desc, cluster_name])
            if "_" in record.id:
                genome_id = record.id.rsplit("_", 1)[0]
                genome_ids.append(genome_id)
                all_genomes.add(genome_id)
                
        desc_dict[cluster_name] = " | ".join([f"{d} ({c})" for d, c in desc_counter.items()]) if desc_counter else "Unknown"
        cluster_records.append((cluster_name, genome_ids))

    pd.DataFrame(summary_csv_rows, columns=["protein", "description", "cluster_id"]).to_csv(output_csv, index=False)
    return desc_dict, cluster_records, all_genomes


def create_sparse_matrix(cluster_records, all_genomes, output_prefix, desc_dict, robustness_df=None, logger=None, binary=False):
    sorted_genomes = sorted(list(all_genomes))
    genome_to_idx = {g: i for i, g in enumerate(sorted_genomes)}
    
    if logger: logger.info(f"Building SciPy native sparse matrix: {output_prefix}_matrix.npz")

    rob_cols = [c for c in robustness_df.columns if c != "cluster_id"] if robustness_df is not None else []
    robustness_dict = robustness_df.set_index("cluster_id").to_dict("index") if robustness_df is not None else {}

    row_idx, col_idx, data, row_metadata = [], [], [], []

    for r_id, (cluster_name, genome_ids) in enumerate(tqdm(cluster_records, desc="Encoding sparse rows", disable=TQDM_DISABLE)):
        row_dict = {"cluster_id": cluster_name, "description": desc_dict.get(cluster_name, "Unknown")}
        for col in rob_cols: row_dict[col] = robustness_dict.get(cluster_name, {}).get(col, 0)
        
        # Synchronize all genome counts directly with the parsed .faa files
        import collections
        genome_counts = collections.Counter(genome_ids)
        
        # The number of unique genomes is the number of keys in the counter
        row_dict["num_unique_genomes"] = len(genome_counts)
        
        # The max paralogs is the highest count value found for any single genome
        row_dict["max_paralogs_per_genome"] = max(genome_counts.values()) if genome_counts else 0

        
        row_metadata.append(row_dict)

        if binary:
            for g in set(genome_ids):
                if g in genome_to_idx:
                    row_idx.append(r_id); col_idx.append(genome_to_idx[g]); data.append(1)
        else:
            for g, count in collections.Counter(genome_ids).items():
                if g in genome_to_idx:
                    row_idx.append(r_id); col_idx.append(genome_to_idx[g]); data.append(count)

    # Save Sparse Matrix
    sparse_mat = sp.coo_matrix((data, (row_idx, col_idx)), shape=(len(cluster_records), len(sorted_genomes))).tocsr()
    sp.save_npz(f"{output_prefix}_matrix.npz", sparse_mat)

    # Save Row Index (Clusters + Metadata)
    pd.DataFrame(row_metadata).to_csv(f"{output_prefix}_rows.csv", index=False)

    # Save Column Index (Genomes)
    with open(f"{output_prefix}_cols.txt", "w") as f: f.write("\n".join(sorted_genomes) + "\n")


def extract_core_clusters(binary_prefix, core_output_csv, logger=None):
    mat_path, rows_path = f"{binary_prefix}_matrix.npz", f"{binary_prefix}_rows.csv"
    if not os.path.exists(mat_path) or not os.path.exists(rows_path): return

    if logger: logger.info("Extracting core clusters from native sparse matrix...")
    mat = sp.load_npz(mat_path)
    
    # Row sums. If sum == num_genomes, it is core.
    core_indices = np.where(mat.sum(axis=1).A1 == mat.shape[1])[0]
    
    rows_df = pd.read_csv(rows_path)
    core_df = rows_df.iloc[core_indices]
    core_df[["cluster_id", "description"]].to_csv(core_output_csv, index=False)
    if logger: logger.info(f"Found {len(core_df)} core clusters. Saved to {core_output_csv}")


# Legacy - could use to make rank assignments static between runs
#def extent_hash(extent_array):
#    # Stable hash of sorted genome index array
#    return hashlib.blake2b(extent_array.tobytes(), digest_size=16).hexdigest()

def load_and_calculate_distances(binary_prefix, logger=None):
    mat_path, cols_path = f"{binary_prefix}_matrix.npz", f"{binary_prefix}_cols.txt"
    if not HAS_SCALABLE_VIS or not os.path.exists(mat_path): return None, None, None, None, None
    
    if logger: logger.info("Loading sparse matrix and calculating exact Jaccard and Sørensen-Dice distances...")
    with open(cols_path, "r") as f: genomes = [line.strip() for line in f if line.strip()]
    num_genomes = len(genomes)
    
    # Load and transpose to (Genomes x Clusters) for the dot product
    A = sp.load_npz(mat_path).transpose().tocsr().astype(np.float32)
    row_sums = A.sum(axis=1).A1 
    
    square_dist_jaccard = np.zeros((num_genomes, num_genomes), dtype=np.float32)
    square_dist_dice = np.zeros((num_genomes, num_genomes), dtype=np.float32)
    chunk_size = 2000
    
    for i in tqdm(range(0, num_genomes, chunk_size), desc="Calculating distance batches", disable=TQDM_DISABLE):
        end = min(i + chunk_size, num_genomes)
        int_chunk = A[i:end].dot(A.T).toarray()
        
        # Jaccard denominator: |A| + |B| - |A n B|
        union_chunk = row_sums[i:end, None] + row_sums - int_chunk
        
        # Dice denominator: |A| + |B|
        sum_chunk = row_sums[i:end, None] + row_sums
        
        with np.errstate(divide='ignore', invalid='ignore'):
            jaccard_chunk = int_chunk / union_chunk
            jaccard_chunk[np.isnan(jaccard_chunk)] = 0.0
            
            dice_chunk = (2.0 * int_chunk) / sum_chunk
            dice_chunk[np.isnan(dice_chunk)] = 0.0
            
        square_dist_jaccard[i:end, :] = 1.0 - jaccard_chunk
        square_dist_dice[i:end, :] = 1.0 - dice_chunk
        
    np.fill_diagonal(square_dist_jaccard, 0.0)
    np.fill_diagonal(square_dist_dice, 0.0)
    
    square_dist_jaccard = np.clip(square_dist_jaccard, 0.0, 1.0)
    square_dist_dice = np.clip(square_dist_dice, 0.0, 1.0)
    
    return genomes, A, square_dist_jaccard, squareform(square_dist_jaccard, checks=False), square_dist_dice

def run_phage_hierarchy_analysis(binary_prefix, output_dir, min_shared_clusters=2, metadata_file=None, colour_col="Family", network_cols="all", logger=None, square_dist_jaccard=None, square_dist_dice=None):
    if not HAS_NETWORKIT:
        if logger: logger.error("NetworKit not installed. Skipping hierarchy.")
        return
        
    if square_dist_jaccard is None or square_dist_dice is None:
        if logger: logger.error("Distance matrices are required for Complete Linkage clustering. Skipping hierarchy.")
        return
    
    mat_path = f"{binary_prefix}_matrix.npz"
    cols_path = f"{binary_prefix}_cols.txt"
    rows_path = f"{binary_prefix}_rows.csv"
    
    with open(cols_path) as f: 
        genomes = np.array([l.strip() for l in f if l.strip()])
    rows_df = pd.read_csv(rows_path)
    A = sp.load_npz(mat_path).transpose().tocsr().astype(np.float32)
    num_genomes = A.shape[0]

    # Load the corresponding count matrix to get the total number of proteins per genome
    count_mat_path = mat_path.replace("_binary", "_count")
    if os.path.exists(count_mat_path):
        count_mat = sp.load_npz(count_mat_path).transpose().tocsr()
        protein_counts = count_mat.sum(axis=1).A1
        genome_prot_counts = {genomes[i]: int(protein_counts[i]) for i in range(num_genomes)}
    else:
        genome_prot_counts = collections.defaultdict(int)

    # Loading of external metadata for export of network
    meta_dict = {}
    target_cols = []
    
    if metadata_file and os.path.exists(metadata_file):
        try:
            m_df = pd.read_csv(metadata_file, dtype=str)
            if "Genome" not in m_df.columns:
                m_df.rename(columns={m_df.columns[0]: "Genome"}, inplace=True)
                
            # Determine which columns to keep based on the CLI input
            if network_cols.lower() == "all":
                target_cols = [c for c in m_df.columns if c != "Genome"]
            else:
                target_cols = [c.strip() for c in network_cols.split(",") if c.strip() in m_df.columns]
                
            # Convert to a nested dictionary: {genome_id: {col1: val1, col2: val2}}
            m_df.set_index("Genome", inplace=True)
            meta_dict = m_df.to_dict(orient="index")
            
        except Exception as e:
            if logger: logger.warning(f"Could not load metadata for graph mapping: {e}")

    def build_graph(indices, min_sim, save_edgelist_path=None):
        sub_A = A[indices]
        G = nk.graph.Graph(len(indices), weighted=True)
        r_sums = sub_A.sum(axis=1).A1
        inter = sub_A.dot(sub_A.T).toarray()
        
        edges_for_export = []
        connected_nodes = set()
        
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                intersection = inter[i, j]
                
                if intersection >= min_shared_clusters:
                    union = r_sums[i] + r_sums[j] - intersection
                    sim = intersection / union if union > 0 else 0
                    
                    if sim >= min_sim: 
                        G.addEdge(i, j, sim)
                        
                        if save_edgelist_path:
                            connected_nodes.update([i, j]) 
                            # These two lines define the variables and must come first
                            source_genome = genomes[indices[i]]
                            target_genome = genomes[indices[j]]
                            
                            # Base edge data
                            edge_data = {
                                "Source": source_genome,
                                "Target": target_genome,
                                "Weight": round(sim, 4),
                                "Shared_Clusters": int(intersection)
                            }
                            
                            # Dynamically add all requested metadata columns
                            for col in target_cols:
                                edge_data[f"Source_{col}"] = meta_dict.get(source_genome, {}).get(col, "Unknown")
                                edge_data[f"Target_{col}"] = meta_dict.get(target_genome, {}).get(col, "Unknown")
                                
                            edges_for_export.append(edge_data)
        
        if save_edgelist_path and edges_for_export:
            for i in range(len(indices)):
                if i not in connected_nodes:
                    source_genome = genomes[indices[i]]
                    edge_data = {
                        "Source": source_genome,
                        "Target": "",  # Blank target creates a singleton in Cytoscape
                        "Weight": 0.0,
                        "Shared_Clusters": int(r_sums[i]) # Stores total proteins for the singleton
                    }
                    for col in target_cols:
                        edge_data[f"Source_{col}"] = meta_dict.get(source_genome, {}).get(col, "Unknown")
                        edge_data[f"Target_{col}"] = ""
                    edges_for_export.append(edge_data)

            if edges_for_export:
                pd.DataFrame(edges_for_export).to_csv(save_edgelist_path, index=False)
                if logger: logger.info(f"Exported {len(edges_for_export)} network elements (including singletons) to {save_edgelist_path}")
            
        return G

    def get_core(indices, threshold=0.90):
        group_mat = A[indices]
        return np.where(group_mat.sum(axis=0).A1 >= (threshold * len(indices)))[0]

    def get_partition_map(graph):      
        cc = nk.components.ConnectedComponents(graph)
        cc.run()       
      
        return cc.getPartition().subsetMap()

    lineage_res, sig_res = [], []
    if logger: logger.info("Building Strict Phage Taxonomy via Hybrid Network/Complete Linkage...")

    # Broad Network Sweep (Find disconnected "Islands")
    root_indices = list(range(num_genomes))
    global_network_file = os.path.join(output_dir, "phage_global_network_edges.csv")
    
    # Use a permissive Jaccard (0.05) to group distant relatives into islands
    islands_map = get_partition_map(build_graph(root_indices, 0.05, save_edgelist_path=global_network_file))    
    global_fam_counter = 1
    global_subfam_counter = 1
    global_genus_counter = 1

    # Parse and split families using complete Linkage
    for island_id, island_local_idx in islands_map.items():
        if len(island_local_idx) < 1: continue
        island_idx = [root_indices[i] for i in island_local_idx]
        
    # If it is a single genome, it represents its own lineage
        if len(island_idx) == 1:
            f_label = f"Family_{global_fam_counter}"
            sf_label = f"{f_label}_SF_{global_subfam_counter}"
            g_label = f"{sf_label}_G_{global_genus_counter}"
            
            # Calculate actual core genes for the singleton
            single_idx = [island_idx[0]]
            singleton_core = get_core(single_idx, threshold=0.90)
            core_count = len(singleton_core)
            
            # Record signatures for each unique lineage
            for ci in singleton_core:
                c_id = rows_df.iloc[ci]["cluster_id"]
                desc = rows_df.iloc[ci]["description"]
                sig_res.append({"Rank": "Family", "Lineage": f_label, "Cluster": c_id, "Desc": desc})
                sig_res.append({"Rank": "Subfamily", "Lineage": sf_label, "Cluster": c_id, "Desc": desc})
                sig_res.append({"Rank": "Genus", "Lineage": g_label, "Cluster": c_id, "Desc": desc})
            
            lineage_res.append({
                "Genome": genomes[island_idx[0]], "Total_Proteins": genome_prot_counts[genomes[island_idx[0]]], "Family": f_label, "Subfamily": sf_label, "Genus": g_label,
                "Family_Core": core_count, "Subfamily_Core": core_count, "Genus_Core": core_count
            })
            global_fam_counter += 1; global_subfam_counter += 1; global_genus_counter += 1
            continue

        # Extract the distance matrices just for this Island
        island_dist_jaccard = square_dist_jaccard[np.ix_(island_idx, island_idx)]
        island_dist_dice = square_dist_dice[np.ix_(island_idx, island_idx)]
        
        # Convert to condensed 1D arrays for SciPy linkage
        island_condensed_jaccard = squareform(island_dist_jaccard, checks=False)
        island_condensed_dice = squareform(island_dist_dice, checks=False)
        
        # Run Complete Linkage on both metrics
        Z_jaccard = linkage(island_condensed_jaccard, method='complete')
        Z_dice = linkage(island_condensed_dice, method='complete')
        
        # Define taxonomic cuts based on the hybrid strategy; DICE for f and sf, Jaccard for g
        fam_labels = fcluster(Z_dice, t=0.95, criterion='distance')      # Dice for deep structural splits
        subfam_labels = fcluster(Z_dice, t=0.70, criterion='distance')   # Dice for deep structural splits
        genus_labels = fcluster(Z_jaccard, t=0.40, criterion='distance') # Jaccard for high-resolution tip splits
        
       
        # ----------------- FAMILY LEVEL -----------------
        fam_map = collections.defaultdict(list)
        for i, f_id in enumerate(fam_labels): fam_map[f_id].append(i)

        for f_id, f_loc in fam_map.items():
            f_idx = [island_idx[i] for i in f_loc]
            f_label = f"Family_{global_fam_counter}"
            global_fam_counter += 1
            f_core = get_core(f_idx, threshold=0.90)
            
            for ci in f_core:
                sig_res.append({"Rank": "Family", "Lineage": f_label, "Cluster": rows_df.iloc[ci]["cluster_id"], "Desc": rows_df.iloc[ci]["description"]})

            # ----------------- SUBFAMILY LEVEL -----------------
            sf_map = collections.defaultdict(list)
            for i in f_loc: sf_map[subfam_labels[i]].append(i)
            
            for sf_id, sf_loc in sf_map.items():
                sf_idx = [island_idx[i] for i in sf_loc]
                sf_label = f"{f_label}_SF_{global_subfam_counter}"
                global_subfam_counter += 1
                sf_core = get_core(sf_idx, threshold=0.90)
                
                for ci in sf_core:
                    sig_res.append({"Rank": "Subfamily", "Lineage": sf_label, "Cluster": rows_df.iloc[ci]["cluster_id"], "Desc": rows_df.iloc[ci]["description"]})

                # ----------------- GENUS LEVEL -----------------
                g_map = collections.defaultdict(list)
                for i in sf_loc: g_map[genus_labels[i]].append(i)
                
                for g_id, g_loc in g_map.items():
                    g_idx = [island_idx[i] for i in g_loc]
                    g_label = f"{sf_label}_G_{global_genus_counter}"
                    global_genus_counter += 1
                    g_core = get_core(g_idx, threshold=0.90)
                    
                    for ci in g_core:
                        sig_res.append({"Rank": "Genus", "Lineage": g_label, "Cluster": rows_df.iloc[ci]["cluster_id"], "Desc": rows_df.iloc[ci]["description"]})
                    
                    for i_val in g_idx:
                        lineage_res.append({
                            "Genome": genomes[i_val], "Total_Proteins": genome_prot_counts[genomes[i_val]], "Family": f_label, "Subfamily": sf_label, "Genus": g_label, 
                            "Family_Core": len(f_core), "Subfamily_Core": len(sf_core), "Genus_Core": len(g_core)
                        })

    # Export Final Taxonomy and Signatures
    tax_df = pd.DataFrame(lineage_res)
    
    # Merge external metadata into the taxonomy export if available
    if metadata_file and os.path.exists(metadata_file):
        try:
            m_df = pd.read_csv(metadata_file, dtype=str)
            if "Genome" not in m_df.columns:
                m_df.rename(columns={m_df.columns[0]: "Genome"}, inplace=True)
                
            # Rename overlapping columns e.g. Family, Genus to Ref_Family, Ref_Genus
            rename_dict = {}
            for col in m_df.columns:
                if col != "Genome" and col in tax_df.columns:
                    rename_dict[col] = f"Ref_{col}"
            
            # Only proceed with the merge if there are columns to compare
            if rename_dict:
                m_df.rename(columns=rename_dict, inplace=True)
                
                # Keep only the Genome column and the reference taxonomy columns
                cols_to_keep = ["Genome"] + list(rename_dict.values())
                m_df_subset = m_df[cols_to_keep]
                
                # Left merge to keep all VORTEX genomes, adding only the comparative taxonomy
                tax_df = pd.merge(tax_df, m_df_subset, on="Genome", how="left")
                tax_df = tax_df.fillna("Unknown")
                if logger: logger.info("Merged comparative reference taxonomy into phage_taxonomy.csv")
            
        except Exception as e:
            if logger: logger.warning(f"Could not merge metadata into taxonomy export: {e}")

    tax_df.to_csv(os.path.join(output_dir, "phage_taxonomy.csv"), index=False)
    pd.DataFrame(sig_res).to_csv(os.path.join(output_dir, "phage_signatures.csv"), index=False)
    if logger: logger.info("Hierarchical analysis complete!")

def generate_jaccard_nexus(genomes, square_dist, output_nex, logger=None):
    if genomes is None or square_dist is None: return
    if logger: logger.info("Generating Jaccard distance Nexus file...")

    taxa_names_clean = [re.sub(r'[^A-Za-z0-9_.]', '_', str(n)) for n in genomes]
    nexus_content = f"#NEXUS\n\nBEGIN TAXA;\n\tDIMENSIONS NTAX={len(genomes)};\n\tTAXLABELS {' '.join(taxa_names_clean)};\nEND;\n\nBEGIN DISTANCES;\n\tFORMAT TRIANGLE=BOTH DIAGONAL LABELS;\n\tMATRIX\n"
    for i, taxon in enumerate(taxa_names_clean):
        nexus_content += f"\t{taxon}\t{' '.join([f'{d:.6f}' for d in square_dist[i]])}\n"
    
    with open(output_nex, "w") as f: f.write(nexus_content + "\t;\nEND;\n")

def linkage_to_newick(linkage_matrix, labels):
    tree, _ = to_tree(linkage_matrix, rd=True)
    def build_newick(node):
        if node.is_leaf(): return f"{labels[node.id]}:{node.dist:.6f}"
        return f"({build_newick(node.left)},{build_newick(node.right)}):{node.dist:.6f}"
    return build_newick(tree) + ";"

def run_visualization(genomes, genomes_by_clusters_sparse, square_dist, condensed_dist, base_out, fig_format, dpi, logger=None, no_vis=False):
    if not HAS_CORE_VIS or genomes is None: return

    row_linkage = linkage(condensed_dist, method="complete", optimal_ordering=True)
    with open(f"{base_out}_complete_dendrogram.newick", "w") as f: f.write(linkage_to_newick(row_linkage, genomes) + "\n")

    if not no_vis and genomes_by_clusters_sparse.shape[1] <= 50000:
        if logger: logger.info("Rendering binary clustermap...")
        sns.set_theme(font_scale=1)
        bin_df = pd.DataFrame(genomes_by_clusters_sparse.toarray(), index=genomes)
        ax1 = sns.clustermap(bin_df, cmap="Blues", row_linkage=row_linkage, col_cluster=True, yticklabels=True, xticklabels=False, figsize=(20, 10))
        ax1.ax_col_dendrogram.set_visible(False)
        ax1.ax_heatmap.set_yticklabels(ax1.ax_heatmap.get_yticklabels(), rotation=0, fontsize=max(1.5, min(10, 100/len(genomes))))
        plt.savefig(f"{base_out}_binary.{fig_format}", format=fig_format, dpi=dpi, bbox_inches="tight"); plt.close("all")

    shared_pct = pd.DataFrame((1.0 - square_dist) * 100.0, index=genomes, columns=genomes)
    with open(f"{base_out}_shared_proteins_dendrogram.newick", "w") as f: f.write(linkage_to_newick(row_linkage, genomes) + "\n")

    if not no_vis:
        if logger: logger.info("Rendering shared-percent clustermap...")
        ax2 = sns.clustermap(shared_pct, cmap="Reds", row_linkage=row_linkage, col_linkage=row_linkage, yticklabels=True, xticklabels=False, figsize=(12, 10))
        ax2.ax_col_dendrogram.set_visible(False)
        ax2.ax_heatmap.set_yticklabels(ax2.ax_heatmap.get_yticklabels(), rotation=0, fontsize=max(1.5, min(10, 100/len(genomes))))
        plt.savefig(f"{base_out}_shared_pct.{fig_format}", format=fig_format, dpi=dpi, bbox_inches="tight"); plt.close("all")
        try: shared_pct.iloc[ax2.dendrogram_row.reordered_ind, ax2.dendrogram_col.reordered_ind].to_csv(f"{base_out}_shared_pct.csv")
        except Exception: shared_pct.to_csv(f"{base_out}_shared_pct.csv")
    else: shared_pct.to_csv(f"{base_out}_shared_pct.csv")


def run_scalable_visualization(genomes, genomes_by_clusters_sparse, base_out, metadata_file=None, colour_col="Family", logger=None):

    if not HAS_SCALABLE_VIS or genomes is None: return

    meta_df, color_kwarg = None, None
    if metadata_file and os.path.exists(metadata_file):
        try:
            meta_df = pd.read_csv(metadata_file, dtype=str)
            if "Genome" not in meta_df.columns: meta_df.rename(columns={meta_df.columns[0]: "Genome"}, inplace=True)
            if logger: logger.info(f"Loaded metadata. Coloring plots by '{colour_col}'...")
        except Exception: meta_df = None

    try:
        if logger: logger.info("Running TruncatedSVD (Sparse PCA)...")
        svd = TruncatedSVD(n_components=3, random_state=42)
        components = svd.fit_transform(genomes_by_clusters_sparse)
        var_ratio = svd.explained_variance_ratio_ * 100

        df_svd = pd.DataFrame(components, columns=["PC1", "PC2", "PC3"])
        df_svd["Genome"] = genomes

        if meta_df is not None:
            df_svd = df_svd.merge(meta_df, on="Genome", how="left")
            if colour_col in df_svd.columns:
                df_svd[colour_col] = df_svd[colour_col].fillna("Unknown")
                color_kwarg = colour_col

        fig_svd = px.scatter_3d(df_svd, x="PC1", y="PC2", z="PC3", hover_name="Genome", color=color_kwarg, title="Viral Lineage Projection (Sparse SVD)", template="plotly_white")
        fig_svd.update_traces(marker=dict(size=5, opacity=0.4))
        fig_svd.write_html(f"{base_out}_lineage_svd_3d.html")
    except Exception as e:
        if logger: logger.error(f"SVD projection failed: {e}")

    if HAS_UMAP:
        try:
            if logger: logger.info("Running UMAP (2D)...")
            reducer = umap.UMAP(n_components=2, random_state=42, metric="jaccard", init="random")
            df_umap = pd.DataFrame(reducer.fit_transform(genomes_by_clusters_sparse), columns=["UMAP1", "UMAP2"])
            df_umap["Genome"] = genomes

            if meta_df is not None: df_umap = df_umap.merge(meta_df, on="Genome", how="left")

            fig_umap = px.scatter(df_umap, x="UMAP1", y="UMAP2", hover_name="Genome", color=color_kwarg, title="Viral Topological Clusters (2D UMAP)", template="plotly_white")
            fig_umap.update_traces(marker=dict(size=6, opacity=0.4))
            fig_umap.write_html(f"{base_out}_lineage_umap_2d.html")
        except Exception as e:
            if logger: logger.error(f"UMAP projection failed: {e}")


def plot_upset_clusters(binary_prefix, output_plot, logger=None):
    mat_path, cols_path = f"{binary_prefix}_matrix.npz", f"{binary_prefix}_cols.txt"
    if not HAS_UPSET or not os.path.exists(mat_path): return

    try:
        with open(cols_path, "r") as f: genome_cols = [line.strip() for line in f if line.strip()]
        num_genomes = len(genome_cols)
        if num_genomes < 2 or num_genomes > 20: return
        if logger: logger.info("Generating UpSet plot via sparse bitmask dot product...")

        mat = sp.load_npz(mat_path).tocsr()
        powers_of_two = (1 << np.arange(num_genomes)).astype(np.int64)
        
        # Use a dense dot product to replaces loops
        bitmasks = mat.dot(powers_of_two)
        global_bit_counts = np.bincount(bitmasks)

        active_masks = np.nonzero(global_bit_counts)[0]
        active_masks = active_masks[active_masks > 0]
        if len(active_masks) == 0: return

        cols_array = np.array(genome_cols)
        memberships, counts = [], []
        for mask in active_masks:
            memberships.append(tuple(cols_array[(mask & powers_of_two) > 0]))
            counts.append(global_bit_counts[mask])

        upset_series = from_memberships(memberships, data=counts)
        temp_df = pd.DataFrame({"count": upset_series.values, "degree": upset_series.index.to_frame().sum(axis=1)}, index=upset_series.index).sort_values(by=["degree", "count"], ascending=[False, False]).head(40)

        try:
            pd.options.mode.copy_on_write = False
        except Exception:
            pass
        
        fig = plt.figure(figsize=(16, 10))
        UpSet(temp_df["count"], show_counts=True, sort_by=None).plot(fig=fig)
        plt.savefig(output_plot, dpi=300, bbox_inches="tight"); plt.close("all")

    except Exception as e:
        if logger: logger.error(f"UpSet plot failed: {e}")


def main():
    parser = argparse.ArgumentParser(description="VORTEX: Viral Orthologous Relationships, Taxonomy and Exploration")
    parser.add_argument("-i", "--input_dir", type=str, help="folder containing GenBank files")
    parser.add_argument("-if", "--input_fasta", type=str, help="path to a multi-fasta file")
    parser.add_argument("-o", "--output_dir", type=str, required=True, help="path to store the results")
    parser.add_argument("-t", "--threads", type=int, default=16, help="Number of threads.")
    parser.add_argument("-mem", "--duckdb_mem", type=str, required=True, default="16GB", help="Amount of RAM in GB allocated to database, e.g. 16GB")
    parser.add_argument("-ids", "--min_seq_ids", type=str, default="0.7", help="comma-separated list of minimum sequence identities to assess, e.g. 0.3,0.5,0.7")
    parser.add_argument("-c", "--coverage", type=float, default=0.5, help="minimum coverage threshold to pass to MMSeqs2, e.g. 0.8 = 80%% coverage")
    parser.add_argument("-min", "--min_shared_clusters", type=int, default=2, help="minimum number of core shared protein clusters to define network groups.")
    parser.add_argument("-f", "--fig_format", type=str, default="png", choices=["png", "pdf", "svg"], help="format for static figures; binary presence/absence, Jaccard similarity heatmap and UpSet plot")
    parser.add_argument("--force_jaccard", action="store_true", help="Force the generation of the dense Jaccard distance matrix and Nexus format dendrogram regardless of the number of input genomes")
    parser.add_argument("--dpi", type=int, default=300, help="dpi resolution for png output files")
    parser.add_argument("--metadata", type=str, default=None, help="csv metadata file to annotated UMAP and SVD outputs, must contain genome accessions")
    parser.add_argument("--colour_by", type=str, default="Family", help="column in metadata file used to colour genomes in UMAP and SVD outputs")
    parser.add_argument("--network_cols", type=str, default="all", help="comma-separated list of metadata columns to map to the network export, or 'all' to map every column")
    parser.add_argument("--no-vis", action="store_true", help="skip all visualisation steps. Note that static figures are not generated when >1000 genomes are present in the input")
    parser.add_argument("--log-level", type=str, default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Log level.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}", help="Show program's version number and exit")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    
    # Setup Global Logs Directory
    global_logs_dir = os.path.join(args.output_dir, "logs")
    os.makedirs(global_logs_dir, exist_ok=True)
    logger = setup_logger(os.path.join(global_logs_dir, "run.log"), getattr(logging, args.log_level.upper(), logging.INFO))

    if not args.input_dir and not args.input_fasta:
        logger.error("Must provide either --input_dir or --input_fasta."); sys.exit(1)

    output_fasta = os.path.join(args.output_dir, "all_proteins.fasta")
    if args.input_fasta:
        logger.info("Processing input FASTA and filtering out duplicate IDs...")
        filter_duplicate_fasta(args.input_fasta, output_fasta, logger)
    elif args.input_dir:
        gb_files = sorted([f for f in os.listdir(args.input_dir) if f.endswith(".gb") or f.endswith(".gbk")])
        if not gb_files: logger.error("No GenBank files found."); sys.exit(1)
        with open(output_fasta, "w") as out_handle:
            for f in tqdm(gb_files, disable=TQDM_DISABLE):
                # We no longer set protein_counter = 1 here. 
                # extract_proteins handles it internally per-record.
                extract_proteins(os.path.join(args.input_dir, f), out_handle, logger)

    for min_seq_id in [float(x.strip()) for x in args.min_seq_ids.split(",") if x.strip()]:
        seqid_str = str(min_seq_id).replace(".", "_")
        
        # Set up run-specific subfolders
        run_base_dir = os.path.join(args.output_dir, f"mmseqs_seqid_{seqid_str}")
        logs_dir = os.path.join(run_base_dir, "logs")
        results_dir = os.path.join(run_base_dir, "results")
        
        os.makedirs(run_base_dir, exist_ok=True)
        os.makedirs(logs_dir, exist_ok=True)
        os.makedirs(results_dir, exist_ok=True)

        db_path = os.path.join(run_base_dir, "db")
        cluster_db_path = os.path.join(run_base_dir, "cluster_results")
        tmp_dir = os.path.join(run_base_dir, "tmp")
        cluster_tsv_file = os.path.join(run_base_dir, "cluster_mapping.tsv")
        
        # Route logs to logs_dir
        log_files = {
            "createdb": os.path.join(logs_dir, "createdb.log"), 
            "cluster": os.path.join(logs_dir, "cluster.log"), 
            "align": os.path.join(logs_dir, "align.log"), 
            "createtsv": os.path.join(logs_dir, "createtsv.log")
        }

        run_logger = setup_logger(os.path.join(logs_dir, "run.log"), getattr(logging, args.log_level.upper(), logging.INFO))
        run_logger.info(f"Running MMseqs2 (min_seq_id={min_seq_id})")

        if shutil.which("mmseqs") is None: run_logger.error("mmseqs binary not found. Skipping."); continue

        try:
            run_mmseqs(output_fasta, cluster_tsv_file, log_files, args.threads, {"mmseqs2": {"min_seq_id": min_seq_id, "cov_mode": 0, "sensitivity": 7.5, "coverage": args.coverage, "evalue": 0.001, "cluster_mode": 0}}, db_path, cluster_db_path, tmp_dir, run_logger)
        except Exception: continue

        if not os.path.exists(cluster_tsv_file) or os.path.getsize(cluster_tsv_file) == 0: continue

        clusters_dir = os.path.join(run_base_dir, "clusters")
        cluster_mapping, cluster_length_stats = generate_clusters(cluster_tsv_file, output_fasta, clusters_dir, run_logger)

        parquet_dir = os.path.join(run_base_dir, "mmseqs_parquet_data")
        robustness_df = process_mmseqs_parquet_partitioned(os.path.join(db_path, "proteome"), os.path.join(db_path, "proteome"), f"{cluster_db_path}/clusters_aligned", args.threads, args.duckdb_mem, cluster_mapping, cluster_length_stats, parquet_dir, run_logger)

        # Route final analytical files to results_dir
        summary_csv = os.path.join(results_dir, f"protein_cluster_summary_seqid_{seqid_str}.csv")
        desc_dict, cluster_records, all_genomes = build_summary_and_records(clusters_dir, summary_csv, output_fasta, run_logger)

        count_prefix = os.path.join(results_dir, f"protein_cluster_matrix_count_seqid_{seqid_str}")
        create_sparse_matrix(cluster_records, all_genomes, count_prefix, desc_dict, robustness_df, run_logger, binary=False)

        binary_prefix = os.path.join(results_dir, f"protein_cluster_matrix_binary_seqid_{seqid_str}")
        create_sparse_matrix(cluster_records, all_genomes, binary_prefix, desc_dict, robustness_df, run_logger, binary=True)

        core_clusters_csv = os.path.join(results_dir, f"core_protein_clusters_seqid_{seqid_str}.csv")
        extract_core_clusters(binary_prefix, core_clusters_csv, run_logger)
        
        # 1. Capture the new Dice output
        genomes, genomes_by_clusters_sparse, square_dist_jaccard, condensed_dist, square_dist_dice = load_and_calculate_distances(binary_prefix, run_logger)
        
        if genomes is not None:
            num_genomes = len(genomes) 
            if num_genomes < 1000 or args.force_jaccard:
                jaccard_sim_matrix = pd.DataFrame(1.0 - square_dist_jaccard, index=genomes, columns=genomes)
                global_linkage = linkage(condensed_dist, method="complete", optimal_ordering=True)
                ordered_indices = leaves_list(global_linkage)
                jaccard_sim_matrix = jaccard_sim_matrix.iloc[ordered_indices, ordered_indices]
                jaccard_sim_matrix.to_csv(os.path.join(results_dir, f"protein_cluster_matrix_seqid_{seqid_str}.csv"))
            else:
                if run_logger:
                    run_logger.info(f"Dataset contains {num_genomes} genomes (>=1000). Skipping dense Jaccard CSV generation. Use --force_jaccard to override.")
                    
            # 2. Pass both matrices into the hierarchy analysis
            run_phage_hierarchy_analysis(
                binary_prefix, results_dir, min_shared_clusters=args.min_shared_clusters, 
                metadata_file=args.metadata, colour_col=args.colour_by, network_cols=args.network_cols, 
                logger=run_logger, square_dist_jaccard=square_dist_jaccard, square_dist_dice=square_dist_dice
            )

            base_vis = os.path.join(results_dir, f"cluster_vis_seqid_{seqid_str}")
            if len(genomes) > 1000:
                run_scalable_visualization(
                    genomes, genomes_by_clusters_sparse, base_vis, 
                    metadata_file=os.path.join(results_dir, "phage_taxonomy.csv"), 
                    colour_col="Family", 
                    logger=run_logger
                )
            if num_genomes <= 1000 or args.force_jaccard:
                nexus_out = os.path.join(results_dir, f"jaccard_distances_seqid_{seqid_str}.nex")
                generate_jaccard_nexus(genomes, square_dist_jaccard, nexus_out, run_logger)
            else:
                if run_logger: 
                    run_logger.info(f"Dataset contains {num_genomes} genomes (>1000). Skipping dense NEXUS matrix generation.")
                
            if len(genomes) <= 1000: run_visualization(genomes, genomes_by_clusters_sparse, square_dist_jaccard, condensed_dist, base_vis, args.fig_format, args.dpi, run_logger, args.no_vis)
            else: run_scalable_visualization(genomes, genomes_by_clusters_sparse, base_vis, args.metadata, args.colour_by, run_logger)

        if not args.no_vis: 
            plot_upset_clusters(binary_prefix, os.path.join(results_dir, f"cluster_upset_plot_seqid_{seqid_str}.{args.fig_format}"), run_logger)
        
        run_logger.info(f"Completed seqid {min_seq_id}")

if __name__ == "__main__":
    main()
