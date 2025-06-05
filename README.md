# hdock-analysis-pipeline
A fully automated pipeline for parsing, scoring, and visualizing HDOCK protein–peptide docking results
This repository contains “hdock_analysis.py,” a comprehensive Python-based pipeline designed to parse HDOCK docking output (compiled_hdock_results.xlsx) together with peptide/receptor FASTA sequences (seq_fasta.xlsx), compute physicochemical properties (GRAVY and net charge at pH 7.4), generate composite scores (with reduced emphasis on ligand RMSD), perform per-complex RMSD‐based clustering, and produce publication-ready outputs: a “long-form” CSV of all poses, a per-complex summary CSV, and three per-complex plots (Docking Score vs. Normalized RMSD with a consistent colorbar, peptide interface-frequency bar charts, and binary contact maps). The tool also supports optional anchor-residue filtering via local PDBs. It was built to streamline large-scale HDOCK analyses and make downstream data inspection, ranking, and visualization effortless.

# HDOCK Analysis Pipeline

[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/pandas)](https://pypi.org/project/pandas/)  

## Overview

**HDOCK Analysis Pipeline** (`hdock_analysis.py`) is a fully automated Python toolkit for parsing, scoring, and visualizing HDOCK protein–peptide docking results. Starting from a single Excel file (`compiled_hdock_results.xlsx`) and an optional FASTA list (`seq_fasta.xlsx`), the script performs:

1. **Parsing**  
   • Extracts per‐pose metrics (Docking Score, Confidence Score, Ligand RMSD, Interface Residues) from the HDOCK “Summary” sheet.  

2. **Quality Control**  
   • Drops empty or invalid rows, imputes missing scores, and filters Z‐score outliers in Docking Score.

3. **Physicochemical Annotation**  
   • For each peptide, computes GRAVY and net charge at pH 7.4 (Kyte–Doolittle scale).

4. **Composite Scoring**  
   • Normalizes Docking Score, Confidence Score, and RMSD (scaled by √(peptide length)), then computes a weighted composite score (default weights: 0.5, 0.4, 0.1).

5. **Optional Anchor‐Residue Filtering**  
   • If you supply local PDB files (**see “Providing PDBs”** below) and a peptide anchor residue number, the pipeline can drop poses whose anchor residues lie > 4 Å from the receptor.

6. **Clustering**  
   • Performs agglomerative clustering (Ward linkage) on pairwise RMSD among poses (default RMSD cutoff 2.0 Å). Assigns each pose a ClusterID and identifies cluster centroids.

7. **Outputs**  
   • **Long‐form CSV** (`hdock_longform_peptide.csv`): one row per model, with all computed scores and annotations.  
   • **Summary CSV** (`hdock_summary_peptide.csv`): one row per complex, listing best models by Docking Score, Confidence, RMSD, and Composite Score, plus averages and cluster‐centroid info.  
   • **Plots (per‐complex)**:  
     - `*_score_vs_rmsd.png`: Docking Score vs. Normalized RMSD (0–1), with a single, global colorbar for Composite Score (exactly five ticks on X and Y axes plus colorbar).  
     - `*_peptide_iface_freq.png`: Horizontal bar plot of the top peptide residues by interface‐frequency (5 X ticks).  
     - `*_contact_map.png`: Binary contact heatmap (models × peptide residues), with up to 5 ticks on X and Y axes.

## Table of Contents

- [Quick Start](#quick-start)  
- [Installation](#installation)  
- [Usage](#usage)  
  - [Basic Run](#basic-run)  
  - [Using Anchor‐Residue Filtering](#using-anchor-residue-filtering)  
  - [Changing Weights, Cutoffs, and Chains](#changing-weights-cutoffs-and-chains)  
- [Output Files](#output-files)  
- [Providing PDB Files (Anchor Filtering)](#providing-pdb-files-anchor-filtering)  
- [Examples](#examples)  
- [Dependencies](#dependencies)  
- [Citation](#citation)

---

## Quick Start

1. **Clone this repository**  
   ```bash
   git clone https://github.com/SidSin0809/hdock-analysis-pipeline.git
   cd hdock-analysis-pipeline
---
2. **Install dependencies**
   pip install -r requirements.txt
---
3. **Run the pipeline**
   python hdock_analysis.py \
  --input compiled_hdock_results.xlsx \
  --seq_excel seq_fasta.xlsx \
  --output_dir ./results
---
Installation
Create (and activate) a Python 3 virtual environment, then install:
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

Usage
Basic Run
python hdock_analysis.py \
  --input path/to/compiled_hdock_results.xlsx \
  --seq_excel path/to/seq_fasta.xlsx \
  --output_dir ./results

--input (-i): Path to the HDOCK results Excel (compiled_hdock_results.xlsx).

--seq_excel: Path to an Excel file containing FASTA entries (one entry per cell, no header). 

--output_dir (-o): Directory where all CSVs and plots will be saved. It will be created if it doesn’t exist.

Weight Tuning and Cutoff
--weights w_dock w_conf w_rmsd: Three floats that sum to 1.0, controlling how much each z‐score contributes to the composite score. Default: 0.5 0.4 0.1 (lowest emphasis on RMSD).

--cluster_cutoff <float>: RMSD cutoff in Å for clustering (default 2.0).

Changing Peptide Chain ID
By default, the script assumes the peptide is chain B when reading PDBs for anchor filtering. If your peptide uses a different chain ID (e.g. P), pass:
--peptide_chain P

Using Anchor‐Residue Filtering
If you want to filter out poses where a specific peptide residue (e.g. residue 5) lies > 4 Å from the receptor, supply:
--anchor_res 5

**This will cause the script to look for local PDB files as described in Providing PDB Files.**

---

4. Output Files
   After a successful run, your --output_dir will contain:
   
4.1. Long‐Form CSV
hdock_longform_peptide.csv

Columns include:
Complex, Model, DockingScore, ConfidenceScore, Ligandrmsd, InterfaceResidues,
PeptideSeq, PeptideLength, GRAVY, NetCharge, Z_Dock, Z_Conf, Z_RMSD, CompositeScore, ClusterID

4.2. Summary CSV
   hdock_summary_peptide.csv

   One row per unique Complex, with fields:
   BestDockModel, BestDockScore, BestConfModel, BestConfScore, BestRmsdModel, BestRmsd, BestCompositeModel, BestCompositeScore, AvgDockingScore, AvgConfidenceScore, AvgLigandrmsd, AvgGRAVY, AvgNetCharge, TopClusterID, ClusterCentroidModel

Per‐Complex PNGs (for each <Complex>, e.g. 6pb0-233):

<Complex>_score_vs_rmsd.png

<Complex>_peptide_iface_freq.png

<Complex>_contact_map.png

Each scatter (*_score_vs_rmsd.png) has exactly five X/Y ticks and a five‐tick colorbar (global CompositeScore range).

---

5. **Providing PDB Files (Anchor Filtering)**
   If you plan to use --anchor_res, you must supply PDB files in a subdirectory of your output folder. Specifically:

Create a folder called model_pdbs/ inside --output_dir.

mkdir -p results/model_pdbs


For every pose and complex, place two PDB files in results/model_pdbs/:

Peptide–receptor docking pose:

<Complex>_<Model>.pdb

Receptor‐only structure:

<Complex>_receptor.pdb

The script will parse <Complex>_<Model>.pdb, locate the Cα atom of the specified anchor residue on the peptide chain (default chain B), and measure its closest distance to any atom in <Complex>_receptor.pdb. Any pose whose anchor‐to‐receptor distance > 4 Å will be dropped.

Run with anchor filtering:
   python hdock_analysis.py \
  --input compiled_hdock_results.xlsx \
  --seq_excel seq_fasta.xlsx \
  --output_dir ./results \
  --anchor_res 5

---

6. **Examples**
Standard pipeline (no anchor filtering):

python hdock_analysis.py \
  --input compiled_hdock_results.xlsx \
  --seq_excel seq_fasta.xlsx \
  --output_dir ./results

With custom weights and RMSD cutoff:

python hdock_analysis.py \
  --input compiled_hdock_results.xlsx \
  --seq_excel seq_fasta.xlsx \
  --output_dir ./results \
  --weights 0.6 0.3 0.1 \
  --cluster_cutoff 1.5

With anchor filtering (peptide chain ‘P’, anchor residue 17):

mkdir -p results/model_pdbs
# (populate results/model_pdbs/ with PDB files as described above)
python hdock_analysis.py \
  --input compiled_hdock_results.xlsx \
  --seq_excel seq_fasta.xlsx \
  --output_dir ./results \
  --peptide_chain P \
  --anchor_res 17

---

7. **requirements**
pandas>=1.2
numpy>=1.19
scipy>=1.5
scikit-learn>=0.24
matplotlib>=3.3
openpyxl>=3.0
biopython>=1.78
---

**If you use this pipeline in your research, please cite**

