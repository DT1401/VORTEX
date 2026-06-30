<p align="center">
  <img src="img/VORTEX.png" alt="vortex logo" width="250">
</p>

# VORTEX - Viral Ortholog Relationships, Taxonomy, and EXploration

## About VORTEX
VORTEX is designed to cluster bacteriophage proteins and enable the examination of relationships between genomes based on shared proteins. To achieve this in a timely manner, VORTEX employs MMSeqs2 cluster and then stores information on the clusters using convertalis. The protein clusters are exported as multi-fasta files and mapped to each genome alongside statistics of protein lengths, minimum and maximum sequence identity. The percentage of protein clusters shared between pairs of genomes is used to calculate the Jaccard distance.

To assign genomes into groups, a large connected component network is created which is then broken down using a minimum Jaccard threshold of 0.05. Heirarchical clustering (complete linkage) of the network is then used to classify the genomes into three levels (loosely corresponding to the taxonomic ranks of family, subfamily and genus). At each level, the number and identity of core protein clusters is determined. 

VORTEX is agnostic to the formal ICTV taxonomy, and does not use the established taxonomy to train assignments. This is deliberate, since ranks up to and including family should demonstrate the presence of core genes that are detectable by standard methods of determining sequence similarity. Therefore, this approach enables an independent assessment of the current taxonomy. 

## Installation
VORTEX can be installed using miniforge. Clone the github repository and use the yaml file to create an environment with the required packages and dependancies
```
conda env create --file environment.yml
conda activate vortex
```

## Usage

### Preparing Input files
#### GenBank files
To provide a standardised input dataset it is recommended that all genomes should be annotated using the same tool with identical parameters (unless you have phages which are known to use alternative translation tables). VORTEX will extract all protein coding regions where each protein has a header of LOCUS_CDS-number and the information contained in the \product qualifier.

#### Protein fasta files
Alternatively, the vConTACT2_proteins.faa file from the excellent [INPHARED](https://github.com/RyanCook94/inphared) can be used as input, with proteins from new genomes appended using `cat`. The header for protein fasta records should follow the format LOCUS_CDS-number and the information contained in the \product qualifier, e.g.
```
>LC190123_001 Large terminase subunit
```

### Running VORTEX

```
python VORTEX.py --input_fasta inphared_proteins/14Apr2025_vConTACT2_proteins.faa --output_dir inphared_03seqid --min_seq_ids 0.3 --coverage 0.8 --min_shared_clusters 4 --metadata inphared_proteins/inphared_test_metadata.csv --color_by Family --duckdb_mem 32GB --threads 12 --fig_format png --dpi 600
```

## Command-Line Parameters
Inputs and the output directory are specified using the flags
```
-i, --input_dir - path to a directory of GenBank files
-if, --input_fasta - path to a multi-fasta file of protein sequences
-o, --output_dir - path to store the results
```

VORTEX is written to handle jobs in parallel, with data from MMSeqs convertalis streamed to a DuckDB database using PyArrow. Specify the number of threads using `-t, --threads` The amount of RAM for the database streaming is specified using ``-mem, --duckdb_mem``, e.g. 
```
--duckdb_mem 16GB
```

Multiple minimum sequence identity thresholds `-ids, --min_seq_ids` can be used by providing a comma-separated list, e.g.
```
--min_seq_ids 0.3,0.5,0.7
```
A separate subdirectory for each sequence identity threshold is created in the output folder. Each minimum sequence identity threshold will be run using the minimum coverage threshold specified using the `-c, --coverage` flag

The `-min, --min_shared_clusters` flag specified the minimum number of shared protein clusters required to consider two phages to be linked at the highest rank assignment. The default value is `4`.

The `-m, --metadata` and `--color_by` flags allow metadata to be mapped to the network. The metadata file should be a csv file with a header row where the first column is labelled `Genome` and contains a list of accessions (omitting the version number). 

Figures can be saved in either `.png`, `.pdf` or `.svg` format, e.g.
```
--fig_format svg
```
The `--dpi` flag allows the user the specify the quality if `.png` outputs are chosen. Alternatively visualisations can be skipped entirely using `--no_vis`.


### Output files
The outputs provided by VORTEX differ based on the number of genomes present within the input files. This has been implemented to avoid excessive calculation times and generation of large (and unreadable) figures.

### Cases where the number of genomes is <1000
`cluster_vis_seqid_binary.{png/tiff/svg}`  
Binary presence-absence heatmap where genomes are rows and columns represent protein clusters.

`cluster_vis_seqid_shared_pct.{png/tiff/svg}`  
Heirarchically clusterd heatmap of jaccard similarity based on the percentage of shared protein clusters between genomes.

`jaccard_distances_seqid.phy`  
A tree of the heirarchically clustered jaccard distances in phylip format. 

`phage_global_network_edges.csv`  
Export of the network for import into tools such as [Cytoscape](http://cytoscape.org/) or [Gephi](https://gephi.org/). Genomes are nodes and two edge weights are provided for each pair of genomes, the jaccard distance and the number of shared protein clusters. If a metadata file and the `color_by` flag have been specified, this information is provided as two additional columns for source and target nodes.

`phage_taxonomy.csv`  
A CSV files that maps each genome to three ranks alongside the number of core proteins determined at each rank.

`phage_signatures.csv`  
A CSV file that details each of the core protein clusters identified in the `phage_taxonomy.csv` output. The product descriptions for protein clusters are also provided. Where these differ, counts are provided for each unique product description, e.g.
```
Family_1    cluster_1234    small terminase subunit (6) | hypothetical protein (4)
```

`protein_cluster_summary_seqid.csv`  
A csv file mapping all proteins to their respective protein cluster.


#### Sparse matrices  
By default, the matrices of shared proteins and jaccard similarities are stored as sparse matrices to avoid the creation of extremely large csv output files.

`protein_cluster_matrix_binary_seqid`  
The binary matrix presents protein clusters as rows and genomes as columns, where the presence or absence of a protein is denoted in binary as 1 or 0, respectively. 

Additional metrics include the number of unique sequences, maximum number of paralogs as well as the minimum, maximum, average and standard deviation of percentage identity and amino acid lengths for proteins within each cluster.

`protein_cluster_matrix_seqid`
An all-vs-all matrix of the Jaccard similarity based on the percentage of shared protein clusters between genomes. 

### Cases where the number of genomes is >1000
Where the number of input genomes exceeds 1000, VORTEX automatically skips the generation of visualisation files and the phylip-format dendrogram. All other outputs are generated.

## Utility Scripts
Where large-scale analysis of >1000 genomes are performed, it can be useful to be able to downsample these results to return information on a specific group of phages.
Rather than re-running VORTEX, two utility scripts are provided to achieve this. 

### shared_protein_clusters.py
The shared protein clusters script takes a list of accession numbers and returns just the protein clusters associated with the specified genomes. This enables the identification of all phages which share one or more protein clusters with the target genomes specified.

Command line arguments:
```
--matrix_prefix", required=True, help="Prefix of the binary matrix files (e.g., results/protein_cluster_matrix_binary_seqid_0_7)
--summary_csv - Path to the protein cluster summary CSV (e.g., results/protein_cluster_summary_seqid_0_7.csv)
--accessions - Text file containing genome accessions to query (one per line, no header)
--out_dir - Path to save the results
```

### subset_cluster_results.py
The results of a large-scale analysis can be downsampled using a list of accession numbers. This script will re-generate the output files to only include the genomes provided. Note however that the metrics provided in the binary matrix output will not be changed.

Command line arguments:
```
--matrix_prefix - Prefix of the binary matrix files (e.g., results/protein_cluster_matrix_binary_seqid_0_7) 
--cluster_dir - Path to the original clusters directory containing .faa files
--out_dir - Path to save the output to
--taxonomy_csv - Path to phage_taxonomy.csv (optional)    
--signatures_csv - Path to phage_signatures.csv (optional)
--fig_format - format to save figures (png|pdf|svg)
--dpi - dots per inch of figure (default = 300)
--accessions - Text file containing genome accessions to keep (one per line, no header)
```

## To-do
> - Simply input flags, allow mixed input types e.g. fasta + GenBank?
> - Map all annotation columns to network output, not just the once specified
> - Add utility script to scrape VMR and generate a metadata file
> - Add utility script to calculate and plot intra- and inter-taxon metrics
> - Add utility script to convert sparse matrix into standard csv
> - Change output names to avoid 'taxonomy'
