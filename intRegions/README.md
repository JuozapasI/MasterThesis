# Intergenic Read Clustering Pipeline

This Nextflow pipeline extracts unassigned reads from a BAM file, filters those that do not overlap with known transcriptomic annotations (GTF),
clusters them, and returns a list of such intergenic regions with a distances to the closest gene.

---

## Requirements

A single conda environment file containing all required packages (exept basic command line tools) is given in (`env.yaml`).
To use it run 

```bash
conda create -f env.yaml
conda activate intRegions
```
---

## Usage

Basic usage:

```bash
nextflow run intergenic.nf \
  --bam path/to/input.bam \
  --gtf path/to/annotations.gtf
```

All parameters:

| Parameter              | Description                                                                 | Default        |
|------------------------|-----------------------------------------------------------------------------|----------------|
| `--bam`                | **(Required)** Input BAM file                                               | —              |
| `--gtf`                | **(Required)** Transcriptomic reference GTF file                            | —              |
| `--cpm_threshold`      | Cluster size threshold in counts-per-million (CPM)                          |`5`             |
|                        | (clusters with number of reads less than this value will be filtered out)   |                |
| `--threshold`          | Absolute cluster size threshold (overrides CPM-based threshold if provided) | `null`         |
| `--read_length`        | Maximum allowed read length (to prevent spliced reads spanning to create    | `1000`         |
|                        | artificial clusters or merge different clusters)                            |                |
| `--unassigned_pattern` | Pattern to detect unassigned reads in bam files (e.g. `GN:Z:-`)             | `"GN:Z:-"`     |

## Output format

The output file (`annotated_clusters.tsv`) is saved in the working directory.
The format is following:

| Column # | Content                                                   |
|----------|-----------------------------------------------------------|
| 1        | Chromosome                                                |
| 2        | Start position                                            |
| 3        | End position                                              |
| 4        | `"."` (placeholder)                                       |
| 5        | Number of reads in the cluster                            |
| 6        | Strand (`+` or `-`)                                       |
| 7        | Closest gene upstream (same strand)                       |
| 8        | Distance to closest gene upstream (same strand)           |
| 9        | Closest gene downstream (same strand)                     |
| 10       | Distance to closest gene downstream (same strand)         |
| 11       | Closest gene upstream (opposite strand)                   |
| 12       | Distance to closest gene upstream (opposite strand)       |
| 13       | Closest gene downstream (opposite strand)                 |
| 14       | Distance to closest gene downstream (opposite strand)     |
