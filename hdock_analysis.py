#!/usr/bin/env python3
"""
hdock_analysis5.py

An HDOCK protein–protein/peptide analysis pipeline that:
  • Parses multiple "<PDB>-<LigID>" complexes from compiled_hdock_results.xlsx
  • Reads receptor/peptide sequences from seq_fasta.xlsx
  • Calculates peptide GRAVY & NetCharge at pH 7.4
  • Computes composite scores with reduced emphasis on ligand RMSD
  • Performs per-complex clustering based on RMSD
  • Outputs: long-form CSV, summary CSV, and per-complex plots (Score vs Norm RMSD, Interface‐freq, Contact Map)

Usage:
    python hdock_analysis5.py \
        --input compiled_hdock_results.xlsx \
        --seq_excel seq_fasta.xlsx \
        --output_dir ./results \
        [--weights 0.5 0.4 0.1] \
        [--cluster_cutoff 2.0] \
        [--peptide_chain B] \
        [--anchor_res 5] \
        [--logfile hdock_pep.log]

Dependencies:
    pandas, numpy, scipy, scikit-learn, matplotlib, openpyxl, argparse, logging, biopython
"""

import os
import sys
import re
import argparse
import logging

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, FormatStrFormatter
from Bio.PDB import PDBParser

# ---------------------------------------------------------
# 1. Logging Setup
# ---------------------------------------------------------
def setup_logging(logfile=None):
    log_format = "%(asctime)s [%(levelname)s] %(message)s"
    if logfile:
        logging.basicConfig(level=logging.INFO, format=log_format, filename=logfile)
    else:
        logging.basicConfig(level=logging.INFO, format=log_format)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------
# 2. Parse HDOCK Results
# ---------------------------------------------------------
def parse_hdock_results(xlsx_path: str, sheet_name: str = "Summary") -> pd.DataFrame:
    """
    Parses HDOCK Excel output into a DataFrame with columns:
      Complex, Model, DockingScore, ConfidenceScore, Ligandrmsd, InterfaceResidues
    """
    logger.info(f"Parsing Excel: {xlsx_path} (sheet: {sheet_name})")
    try:
        df_raw = pd.read_excel(xlsx_path, sheet_name=sheet_name, header=None, dtype=str)
    except Exception as e:
        logger.error(f"Error reading Excel: {e}")
        sys.exit(1)

    nrows, ncols = df_raw.shape
    records = []
    i = 0
    docking_pat = re.compile(r"^Docking\s*Score$", re.IGNORECASE)

    while i < nrows - 4:
        cell0 = df_raw.iat[i, 0]
        if isinstance(cell0, str) and cell0.strip():
            nxt = df_raw.iat[i + 1, 0]
            if isinstance(nxt, str) and docking_pat.match(nxt.strip()):
                complex_name = cell0.strip()
                docking_row = df_raw.iloc[i + 1]
                conf_row    = df_raw.iloc[i + 2]
                rmsd_row    = df_raw.iloc[i + 3]
                int_row     = df_raw.iloc[i + 4]

                # collect model IDs from the “InterfaceResidues” row
                model_ids = []
                for col in range(1, ncols):
                    val = int_row.iat[col]
                    if pd.isna(val):
                        break
                    model_ids.append(str(val).strip())
                if not model_ids:
                    logger.warning(f"No models found for '{complex_name}' at row {i}")
                    i += 1
                    continue

                docking_vals, conf_vals, rmsd_vals = [], [], []
                for col in range(1, 1 + len(model_ids)):
                    # DockingScore
                    try:
                        docking_vals.append(float(docking_row.iat[col]))
                    except:
                        docking_vals.append(np.nan)
                    # ConfidenceScore
                    try:
                        conf_vals.append(float(conf_row.iat[col]))
                    except:
                        conf_vals.append(np.nan)
                    # Ligandrmsd
                    try:
                        rmsd_vals.append(float(rmsd_row.iat[col]))
                    except:
                        rmsd_vals.append(np.nan)

                # InterfaceResidues (e.g. “B:15,A:234,…”) per model
                interface_lists = []
                for col in range(1, 1 + len(model_ids)):
                    raw = int_row.iat[col]
                    if isinstance(raw, str) and raw.strip():
                        residues = [r.strip() for r in raw.split(',') if r.strip()]
                        interface_lists.append(residues)
                    else:
                        interface_lists.append([])

                # build records
                for idx, model in enumerate(model_ids):
                    records.append({
                        'Complex': complex_name,
                        'Model': model,
                        'DockingScore': docking_vals[idx],
                        'ConfidenceScore': conf_vals[idx],
                        'Ligandrmsd': rmsd_vals[idx],
                        'InterfaceResidues': interface_lists[idx]
                    })
                i += 6
                continue
        i += 1

    df_long = pd.DataFrame.from_records(records)
    if df_long.empty:
        logger.error("No data parsed; check Excel format.")
        sys.exit(1)
    logger.info(f"Parsed {len(df_long)} pose entries.")
    return df_long

# ---------------------------------------------------------
# 3. Quality Control
# ---------------------------------------------------------
def quality_control(df: pd.DataFrame, z_thresh: float = 3.0) -> pd.DataFrame:
    df_qc = df.copy()

    # Drop rows where DockingScore, ConfidenceScore, Ligandrmsd are all NaN
    mask_all_na = df_qc[['DockingScore','ConfidenceScore','Ligandrmsd']].isna().all(axis=1)
    if mask_all_na.sum():
        logger.info(f"Dropping {mask_all_na.sum()} fully-missing rows.")
        df_qc = df_qc.loc[~mask_all_na].reset_index(drop=True)

    # Impute missing ConfidenceScore → 0.0
    na_conf = df_qc['ConfidenceScore'].isna().sum()
    if na_conf:
        logger.info(f"Imputing {na_conf} missing ConfidenceScore(s) with 0.0.")
        df_qc.loc[df_qc['ConfidenceScore'].isna(), 'ConfidenceScore'] = 0.0

    # Impute missing Ligandrmsd → (max + 1.0)
    max_r = df_qc['Ligandrmsd'].max(skipna=True)
    na_rmsd = df_qc['Ligandrmsd'].isna().sum()
    if na_rmsd:
        fill_val = (max_r if pd.notna(max_r) else 0.0) + 1.0
        logger.info(f"Imputing {na_rmsd} missing Ligandrmsd(s) with {fill_val:.2f}.")
        df_qc.loc[df_qc['Ligandrmsd'].isna(), 'Ligandrmsd'] = fill_val

    # Drop outliers by Z-score of DockingScore within each Complex
    drop_idx = []
    for cplx, grp in df_qc.groupby('Complex'):
        scores = grp['DockingScore'].astype(float)
        if len(scores) < 2:
            continue
        z = np.abs(stats.zscore(scores, nan_policy='omit'))
        outliers = grp.index[z > z_thresh]
        if len(outliers):
            logger.info(f"Complex '{cplx}': dropping {len(outliers)} outlier(s) by Z>{z_thresh}.")
            drop_idx.extend(outliers.tolist())

    if drop_idx:
        df_qc = df_qc.drop(index=drop_idx).reset_index(drop=True)
    return df_qc

# ---------------------------------------------------------
# 4. Compute GRAVY & NetCharge at pH 7.4
# ---------------------------------------------------------
def compute_peptide_physicochem(df: pd.DataFrame, seq_dict: dict, pH: float = 7.4) -> pd.DataFrame:
    """
    For each Complex in df, look up peptide ID (after the dash) in seq_dict.
    Compute GRAVY (Kyte-Doolittle) and NetCharge at pH.
    """
    # Kyte-Doolittle scale
    kd_scale = {
        'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,'H':-3.2,'I':4.5,
        'K':-3.9,'L':3.8,'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,'R':-4.5,'S':-0.8,
        'T':-0.7,'V':4.2,'W':-0.9,'Y':-1.3
    }
    pka = {'D':3.9,'E':4.1,'C':8.3,'Y':10.1,'H':6.0,'K':10.5,'R':12.5,'N_term':9.0,'C_term':2.0}

    def compute_gravy(seq):
        if not seq:
            return 0.0
        return sum(kd_scale.get(res.upper(), 0.0) for res in seq) / max(len(seq), 1)

    def compute_net_charge(seq, pH):
        pos, neg = 0.0, 0.0
        # N-terminus
        pos += 1.0 / (1.0 + 10 ** (pH - pka['N_term']))
        # C-terminus
        neg += 1.0 / (1.0 + 10 ** (pka['C_term'] - pH))
        for res in seq.upper():
            if res == 'K':
                pos += 1.0 / (1.0 + 10 ** (pH - pka['K']))
            elif res == 'R':
                pos += 1.0 / (1.0 + 10 ** (pH - pka['R']))
            elif res == 'H':
                pos += 1.0 / (1.0 + 10 ** (pH - pka['H']))
            elif res == 'D':
                neg += 1.0 / (1.0 + 10 ** (pka['D'] - pH))
            elif res == 'E':
                neg += 1.0 / (1.0 + 10 ** (pka['E'] - pH))
            elif res == 'C':
                neg += 1.0 / (1.0 + 10 ** (pka['C'] - pH))
            elif res == 'Y':
                neg += 1.0 / (1.0 + 10 ** (pka['Y'] - pH))
        return pos - neg

    # Map each Complex to its peptide sequence
    peptides = {}
    for cplx in df['Complex'].unique():
        lig_id = cplx.split('-')[-1] if '-' in cplx else None
        seq = seq_dict.get(lig_id, '')
        peptides[cplx] = seq

    df_phys = df.copy()
    df_phys['PeptideSeq']    = df_phys['Complex'].map(peptides)
    df_phys['PeptideLength'] = df_phys['PeptideSeq'].map(lambda s: len(s) if s else 1)
    df_phys['GRAVY']         = df_phys['PeptideSeq'].map(lambda s: compute_gravy(s))
    df_phys['NetCharge']     = df_phys['PeptideSeq'].map(lambda s: compute_net_charge(s, pH))

    return df_phys

# ---------------------------------------------------------
# 5. Compute Composite Scores (reduced RMSD weight)
# ---------------------------------------------------------
def compute_composite_scores(df: pd.DataFrame, weights=(0.5, 0.4, 0.1)) -> pd.DataFrame:
    """
    For each Complex, compute z‐scores of:
      - negative DockingScore  (lower docking → higher z)
      - positive ConfidenceScore
      - negative (Ligandrmsd / sqrt(PeptideLength))
    Then CompositeScore = w_dock * Z_dock + w_conf * Z_conf + w_rmsd * Z_rmsd.
    """
    w_dock, w_conf, w_rmsd = weights
    df_scores = df.copy()

    all_z_d, all_z_c, all_z_r = [], [], []
    for cplx, grp in df_scores.groupby('Complex'):
        docks  = grp['DockingScore'].astype(float).values
        confs  = grp['ConfidenceScore'].astype(float).values
        rmsds  = grp['Ligandrmsd'].astype(float).values
        pep_len = grp['PeptideLength'].iloc[0]

        # Normalize RMSD by sqrt(peptide length)
        norm_rmsds = rmsds / np.sqrt(max(pep_len, 1))

        # z‐scores:
        z_d = np.nan_to_num(stats.zscore(-docks, nan_policy='omit'), nan=0.0)
        z_c = np.nan_to_num(stats.zscore(confs, nan_policy='omit'), nan=0.0)
        z_r = np.nan_to_num(stats.zscore(-norm_rmsds, nan_policy='omit'), nan=0.0)

        all_z_d.extend(z_d.tolist())
        all_z_c.extend(z_c.tolist())
        all_z_r.extend(z_r.tolist())

    df_scores['Z_Dock'] = all_z_d
    df_scores['Z_Conf'] = all_z_c
    df_scores['Z_RMSD'] = all_z_r
    df_scores['CompositeScore'] = (
        w_dock * df_scores['Z_Dock'] +
        w_conf * df_scores['Z_Conf'] +
        w_rmsd * df_scores['Z_RMSD']
    )
    return df_scores

# ---------------------------------------------------------
# 6. Clustering by RMSD
# ---------------------------------------------------------
def cluster_models_by_rmsd(df: pd.DataFrame, cutoff: float = 2.0) -> pd.DataFrame:
    """
    Agglomerative clustering (Ward) on pairwise Euclidean distance of Ligandrmsd.
    Assign each pose a ClusterID per Complex.
    """
    df_clust = df.copy()
    df_clust['ClusterID'] = -1

    for cplx, grp in df_clust.groupby('Complex'):
        idxs = grp.index.tolist()
        if len(idxs) < 2:
            df_clust.loc[idxs, 'ClusterID'] = 0
            continue

        rmsd_vals = grp['Ligandrmsd'].astype(float).values.reshape(-1, 1)
        dist_mat = squareform(pdist(rmsd_vals, metric='euclidean'))
        condensed = squareform(dist_mat)
        link = hierarchy.linkage(condensed, method='ward')
        labels = hierarchy.fcluster(link, t=cutoff, criterion='distance')

        for idx, lbl in zip(idxs, labels):
            df_clust.at[idx, 'ClusterID'] = lbl

    return df_clust

# ---------------------------------------------------------
# 7. Interface Frequency Computation (Peptide‐centric)
# ---------------------------------------------------------
def compute_interface_frequencies(df: pd.DataFrame, peptide_chain: str = 'B') -> dict:
    """
    Count how often each peptide residue (e.g. “B:15”) appears in InterfaceResidues
    across all models of a given Complex.
    """
    freq_dict = {}
    for cplx, grp in df.groupby('Complex'):
        counts = {}
        for residues in grp['InterfaceResidues']:
            for r in residues:
                parts = r.split(':')
                if len(parts) == 2 and parts[0] == peptide_chain:
                    resnum = parts[1]
                    counts[resnum] = counts.get(resnum, 0) + 1
        freq_dict[cplx] = counts
    return freq_dict

# ---------------------------------------------------------
# 8. Plotting Functions (with uniform ticks & colourbar range)
# ---------------------------------------------------------
def plot_score_vs_rmsd(df: pd.DataFrame, complex_name: str, output_dir: str,
                       vmin: float, vmax: float):
    """
    Scatter: DockingScore vs Normalized RMSD (0–1), coloured by CompositeScore using fixed [vmin, vmax].
    X‐axis: exactly 5 ticks. Y‐axis: exactly 5 ticks. Colourbar: exactly 5 ticks.
    """
    grp = df[df['Complex'] == complex_name]
    if grp.empty:
        return

    # Normalize RMSD (per‐complex) to [0, 1]
    rmsds = grp['Ligandrmsd'].astype(float)
    min_r, max_r = rmsds.min(), rmsds.max()
    if max_r > min_r:
        norm_rmsd = (rmsds - min_r) / (max_r - min_r)
    else:
        norm_rmsd = np.zeros_like(rmsds)

    grp = grp.copy()
    grp['NormRMSD'] = norm_rmsd

    fig, ax = plt.subplots(figsize=(6, 5))
    sc = ax.scatter(
        grp['DockingScore'],
        grp['NormRMSD'],
        c=grp['CompositeScore'],
        cmap='viridis',
        vmin=vmin,
        vmax=vmax,
        s=40,
        alpha=0.8
    )
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('Composite Score')

    # Force colourbar to have exactly 5 ticks (linearly spaced between vmin/vmax)
    cbar.locator = MaxNLocator(nbins=5)
    cbar.update_ticks()

    ax.set_xlabel('Docking Score')
    ax.set_ylabel('Normalized RMSD')
    ax.set_title(f'{complex_name}: Score vs Normalized RMSD')

    # Uniform tick formatting: exactly 5 ticks on both axes
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    plt.tight_layout()
    out_png = os.path.join(output_dir, f"{complex_name}_score_vs_rmsd.png")
    plt.savefig(out_png, dpi=300)
    plt.close()
    logger.info(f"Saved plot: {out_png}")

def plot_peptide_contact_freq(freq_dict: dict, complex_name: str, output_dir: str):
    """
    Horizontal bar plot of the top 10 peptide residues by frequency of appearance
    in InterfaceResidues. X‐axis forced to 5 ticks.
    """
    counts = freq_dict.get(complex_name, {})
    if not counts:
        return

    sorted_pairs = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    top_pairs = sorted_pairs[:min(10, len(sorted_pairs))]
    residues, vals = zip(*top_pairs)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.barh(residues, vals, color='teal')
    ax.set_xlabel('Frequency (models)')
    ax.set_ylabel('Peptide Residue Number')
    ax.set_title(f'{complex_name}: Top Peptide Interface Residues')
    ax.invert_yaxis()

    # Uniform tick formatting: exactly 5 ticks on X axis
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

    plt.tight_layout()
    out_png = os.path.join(output_dir, f"{complex_name}_peptide_iface_freq.png")
    plt.savefig(out_png, dpi=300)
    plt.close()
    logger.info(f"Saved peptide interface-frequency plot: {out_png}")

def plot_contact_matrix(df: pd.DataFrame, complex_name: str, output_dir: str, peptide_chain: str='B'):
    """
    Binary contact map: rows = models (in their original ranking order),
    columns = peptide residues (sorted). Cells = 1 if that residue appears
    in InterfaceResidues for that model, else 0. Both axes forced to 5 ticks max.
    """
    grp = df[df['Complex'] == complex_name]
    if grp.empty:
        return

    # Gather all unique peptide residue numbers in this Complex
    peptide_residues = sorted(
        {r.split(':')[1] for sub in grp['InterfaceResidues'] for r in sub if r.startswith(peptide_chain + ':')},
        key=lambda x: int(x)
    )
    if not peptide_residues:
        return

    # Build binary matrix: shape = (n_models, n_residues)
    mat = np.zeros((len(grp), len(peptide_residues)), dtype=int)
    for i, residues in enumerate(grp['InterfaceResidues']):
        for r in residues:
            parts = r.split(':')
            if len(parts) == 2 and parts[0] == peptide_chain:
                resnum = parts[1]
                j = peptide_residues.index(resnum)
                mat[i, j] = 1

    fig, ax = plt.subplots(figsize=(6, 4))
    im = ax.imshow(mat, aspect='auto', cmap='Blues', interpolation='nearest')
    ax.set_xlabel('Peptide Residue Number')
    ax.set_ylabel('Models (ranked)')
    ax.set_title(f'{complex_name}: Peptide Contact Map')

    # X‐axis: force at most 5 ticks
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.set_xticklabels([str(r) for r in peptide_residues], rotation=90, fontsize=6)

    # Y‐axis: show at most 5 ticks
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.set_yticklabels(grp['Model'].tolist(), fontsize=6)

    plt.tight_layout()
    out_png = os.path.join(output_dir, f"{complex_name}_contact_map.png")
    plt.savefig(out_png, dpi=300)
    plt.close()
    logger.info(f"Saved contact map: {out_png}")

# ---------------------------------------------------------
# 9. Summary Generation
# ---------------------------------------------------------
def summarize_by_complex(df: pd.DataFrame) -> pd.DataFrame:
    recs = []
    for cplx, grp in df.groupby('Complex'):
        idx_d = grp['DockingScore'].idxmin()
        bm_dock    = grp.at[idx_d, 'Model']
        bd_score   = grp.at[idx_d, 'DockingScore']

        idx_c = grp['ConfidenceScore'].idxmax()
        bm_conf  = grp.at[idx_c, 'Model']
        bc_score = grp.at[idx_c, 'ConfidenceScore']

        idx_r = grp['Ligandrmsd'].idxmin()
        bm_rmsd = grp.at[idx_r, 'Model']
        br_val  = grp.at[idx_r, 'Ligandrmsd']

        idx_comp = grp['CompositeScore'].idxmax()
        bm_comp   = grp.at[idx_comp, 'Model']
        comp_score = grp.at[idx_comp, 'CompositeScore']

        avg_d = grp['DockingScore'].mean()
        avg_c = grp['ConfidenceScore'].mean()
        avg_r = grp['Ligandrmsd'].mean()
        avg_g = grp['GRAVY'].mean()       if 'GRAVY' in grp else np.nan
        avg_n = grp['NetCharge'].mean()   if 'NetCharge' in grp else np.nan

        top_cl = None
        cent_model = None
        if 'ClusterID' in grp.columns:
            top_cl = int(grp['ClusterID'].value_counts().idxmax())
            members = grp[grp['ClusterID'] == top_cl]
            idx_cent = members['Ligandrmsd'].idxmin()
            cent_model = members.at[idx_cent, 'Model']

        recs.append({
            'Complex': cplx,
            'BestDockModel': bm_dock, 'BestDockScore': bd_score,
            'BestConfModel': bm_conf, 'BestConfScore': bc_score,
            'BestRmsdModel': bm_rmsd, 'BestRmsd': br_val,
            'BestCompositeModel': bm_comp, 'BestCompositeScore': comp_score,
            'AvgDockingScore': avg_d, 'AvgConfidenceScore': avg_c, 'AvgLigandrmsd': avg_r,
            'AvgGRAVY': avg_g, 'AvgNetCharge': avg_n,
            'TopClusterID': top_cl, 'ClusterCentroidModel': cent_model
        })

    df_sum = pd.DataFrame.from_records(recs).set_index('Complex')
    return df_sum

# ---------------------------------------------------------
# 10. Main Pipeline
# ---------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="HDOCK Protein–Peptide Pipeline (Enhanced)")
    parser.add_argument('--input', '-i', required=True, help="Path to compiled_hdock_results.xlsx")
    parser.add_argument('--seq_excel', type=str, default=None, help="Path to Excel with FASTA entries (no header)")
    parser.add_argument('--output_dir', '-o', required=True, help="Directory for outputs (CSV, PNG)")
    parser.add_argument('--weights', nargs=3, type=float, default=[0.5, 0.4, 0.1], 
                        help="Weights: w_dock w_conf w_rmsd (sum to 1.0, reduce RMSD emphasis)")
    parser.add_argument('--cluster_cutoff', type=float, default=2.0, help="RMSD cutoff (Å) for clustering")
    parser.add_argument('--peptide_chain', type=str, default='B', help="Chain ID of peptide in PDBs")
    parser.add_argument('--anchor_res', type=int, default=None, help="Peptide residue number for anchor filtering")
    parser.add_argument('--logfile', type=str, default=None, help="Optional logfile path")
    args = parser.parse_args()

    setup_logging(args.logfile)
    if not os.path.isfile(args.input):
        logger.error(f"Input not found: {args.input}")
        sys.exit(1)
    os.makedirs(args.output_dir, exist_ok=True)

    # 1) Read FASTA sequences (if provided) into seq_dict
    seq_dict = {}
    if args.seq_excel:
        if not os.path.isfile(args.seq_excel):
            logger.error(f"Sequence Excel not found: {args.seq_excel}")
            sys.exit(1)
        try:
            df_seq_raw = pd.read_excel(args.seq_excel, header=None, dtype=str)
        except Exception as e:
            logger.error(f"Failed to read sequence Excel: {e}")
            sys.exit(1)

        for idx, row in df_seq_raw.iterrows():
            cell = row[0]
            if not isinstance(cell, str) or not cell.startswith('>'):
                continue
            parts = cell.split('\n', 1)
            if len(parts) != 2:
                continue
            lig_id = parts[0].lstrip('>')
            seq    = parts[1].strip()
            seq_dict[lig_id] = seq
        logger.info(f"Parsed {len(seq_dict)} FASTA entries from {args.seq_excel}")

    # 2) Parse the HDOCK results
    df_long = parse_hdock_results(args.input)
    if seq_dict:
        missing = []
        for cplx in df_long['Complex'].unique():
            lig_id = cplx.split('-')[-1] if '-' in cplx else None
            if not lig_id or lig_id not in seq_dict:
                missing.append(cplx)
        if missing:
            logger.warning(f"Complexes missing FASTA sequences: {', '.join(missing)}; will fallback to interface length.")

    # 3) Quality control
    df_qc = quality_control(df_long)

    # 4) Compute GRAVY & NetCharge
    df_phys = compute_peptide_physicochem(df_qc, seq_dict, pH=7.4)

    # 5) Compute CompositeScore with weights
    df_scored = compute_composite_scores(df_phys, weights=tuple(args.weights))

    # 6) Anchor‐residue filtering (optional)
    if args.anchor_res is not None:
        pdb_dir = os.path.join(args.output_dir, "model_pdbs")
        if not os.path.isdir(pdb_dir):
            logger.warning("--anchor_res provided but no 'model_pdbs' directory found; skipping anchor filter.")
        else:
            parser_pdb = PDBParser(QUIET=True)
            to_drop = []
            for idx, row in df_scored.iterrows():
                cplx  = row['Complex']
                model = row['Model']
                pep_file = os.path.join(pdb_dir, f"{cplx}_{model}.pdb")
                rec_file = os.path.join(pdb_dir, f"{cplx}_receptor.pdb")
                if not (os.path.isfile(pep_file) and os.path.isfile(rec_file)):
                    continue
                try:
                    pep_struc = parser_pdb.get_structure("pep", pep_file)[0]
                    rec_struc = parser_pdb.get_structure("rec", rec_file)[0]
                    pep_chain = list(pep_struc.get_chains())[0]
                    pep_res   = pep_chain[(' ', int(args.anchor_res), ' ')]
                    pep_coord = pep_res['CA'].get_coord() if 'CA' in pep_res else pep_res.child_list[0].get_coord()

                    min_d = float('inf')
                    for rc in rec_struc.get_chains():
                        for r in rc:
                            atom = r['CA'] if 'CA' in r else r.child_list[0]
                            d = np.linalg.norm(atom.get_coord() - pep_coord)
                            min_d = min(min_d, d)
                    if min_d > 4.0:
                        to_drop.append(idx)
                except Exception:
                    continue
            if to_drop:
                logger.info(f"Dropping {len(to_drop)} poses failing anchor contact (res {args.anchor_res}).")
                df_scored = df_scored.drop(index=to_drop).reset_index(drop=True)

    # 7) Clustering by RMSD
    df_clustered = cluster_models_by_rmsd(df_scored, cutoff=args.cluster_cutoff)

    # 8) Compute interface frequencies (peptide‐centric)
    freq_dict = compute_interface_frequencies(df_clustered, peptide_chain=args.peptide_chain)

    # 9) Save long‐form CSV
    long_csv = os.path.join(args.output_dir, 'hdock_longform_peptide.csv')
    df_clustered.to_csv(long_csv, index=False)
    logger.info(f"Long-form CSV saved: {long_csv}")

    # 10) Summarize by Complex
    df_summary = summarize_by_complex(df_clustered)
    summary_csv = os.path.join(args.output_dir, 'hdock_summary_peptide.csv')
    df_summary.to_csv(summary_csv)
    logger.info(f"Summary CSV saved: {summary_csv}")

    # ────────────────────────────────────────────────────────────────────
    # NEW: Compute global vmin, vmax for CompositeScore so that 0 maps
    #      to the same colour in every plot. These will be passed to
    #      plot_score_vs_rmsd() below.
    global_min = df_clustered['CompositeScore'].min()
    global_max = df_clustered['CompositeScore'].max()
    logger.info(f"Global CompositeScore range: [{global_min:.3f}, {global_max:.3f}]")
    # ────────────────────────────────────────────────────────────────────

    # 11) Generate per‐Complex plots
    for cplx in df_clustered['Complex'].unique():
        plot_score_vs_rmsd(df_clustered, cplx, args.output_dir,
                           vmin=global_min, vmax=global_max)
        plot_peptide_contact_freq(freq_dict, cplx, args.output_dir)
        plot_contact_matrix(df_clustered, cplx, args.output_dir, peptide_chain=args.peptide_chain)

    logger.info("Pipeline completed successfully.")

if __name__ == '__main__':
    main()
