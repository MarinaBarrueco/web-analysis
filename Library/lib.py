
#!/usr/bin/env python3
"""
Peptide Analysis Core Library

Core functions for peptide library analysis including:
- Pattern-based peptide filtering and grouping
- CPM (Counts Per Million) normalization and filtering  
- DESeq2-based differential expression analysis
- Statistical testing and significance classification
- Data visualization (MA plots, volcano plots, histograms)
- File I/O and clustering result processing

Author: Peptide Analysis Pipeline
Version: 2.0
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import subprocess
import re
import numpy as np
from matplotlib.patches import Patch
# from rpy2.robjects import pandas2ri
# import rpy2.robjects as ro
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from typing import Any


# get the rute from the current file
CURRENT_DIR = os.path.join("/".join(os.path.dirname(__file__).split("/")[:-1]))

GIBS_CLUSTER_BIN = os.path.join(CURRENT_DIR, "gibbscluster-2.0/gibbscluster")



def group_by_peptide(df: pd.DataFrame, conditions: dict, pattern: str):
    """
    Groups peptides based on a given pattern(the library peptide form) and aggregates their counts.

    This function processes a DataFrame containing peptide sequences and their respective counts
    across multiple samples. It calculates the total peptide count, filters peptides matching a
    specific pattern, and aggregates them based on a cleaned version of the peptide sequence.

    Parameters:
    -----------
    df : pandas.DataFrame
        A DataFrame where each row represents a peptide with its associated counts.
    conditions : dict
        Dictionary mapping sample column names to condition labels (Control/Experiment).
    pattern : str
        A regular expression pattern used to filter and extract peptides.

    Returns:
    --------
    tuple: (pandas.DataFrame, dict)
        A DataFrame grouped by the cleaned peptide sequences, with counts aggregated and sorted,
        and a summary dictionary with statistics.
    """
    
    # Input validation
    if df.empty:
        raise ValueError("Input DataFrame is empty")
    
    if "peptide" not in df.columns:
        raise ValueError("DataFrame must contain a 'peptide' column")
    
    # Validate that all condition columns exist
    missing_cols = [col for col in conditions.keys() if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Condition columns not found in DataFrame: {missing_cols}")
    
    # Create working copy and handle null values
    df = df[["peptide"] + list(conditions.keys())].copy()
    
    # Remove rows with null peptides and replace NaN values with 0
    df = df.dropna(subset=["peptide"])
    df[list(conditions.keys())] = df[list(conditions.keys())].fillna(0)
    
    # Convert peptide column to string and remove obvious invalid entries
    df["peptide"] = df["peptide"].astype(str)
    df = df[df["peptide"] != "nan"]
    df = df[df["peptide"] != ""]
    df = df[df["peptide"].str.len() > 0]
    
    # Calculate total count with validation for zero division
    df["Total Count"] = df[list(conditions.keys())].sum(axis=1)
    
    # Calculate total peptides and valid peptides
    total_peptides = df["Total Count"].sum()
    
    if total_peptides == 0:
        print("Warning: No peptides with non-zero counts found")
        return pd.DataFrame(), {"Total Peptides": 0, "Valid Peptides": 0, "Success Rate": "0.00%"}
    
    # Pattern matching with error handling
    try:
        df["Consistent"] = df["peptide"].str.contains(pattern, regex=True, na=False)
    except re.error as e:
        raise ValueError(f"Invalid regular expression pattern: {pattern}. Error: {e}")
    
    valid_peptides = df.loc[df["Consistent"], "Total Count"].sum()
    
    summary = {
        "Total Peptides": int(total_peptides),
        "Valid Peptides": int(valid_peptides),
        "Success Rate": f"{100 * valid_peptides / total_peptides:.2f}%" if total_peptides > 0 else "0.00%"
    }

    # Extract the clean peptide matching the pattern with error handling
    try:
        df["Clean Peptide"] = df['peptide'].str.extract(f"({pattern})", expand=False)
    except re.error as e:
        raise ValueError(f"Error extracting pattern: {e}")
    
    # From here..
    # Handle amino acid replacements safely
    df["Clean Peptide"] = df['Clean Peptide'].fillna("").astype(str).str.replace("*", "Q", regex=False)

    # Drop rows where Clean Peptide is empty or NaN
    df = df[df["Clean Peptide"].str.len() > 0]
    df = df[df["Clean Peptide"] != "nan"]
    # to here.. makes the difference with the previous version
    if df.empty:
        print("Warning: No peptides matched the specified pattern")
        return pd.DataFrame(), summary

    # Aggregate by clean peptide and sort
    df_grouped = df.groupby("Clean Peptide", as_index=False).sum(numeric_only=True)
    df_grouped = df_grouped.sort_values("Total Count", ascending=False)
    df_grouped.reset_index(inplace=True, drop=True)

    # Drop unnecessary columns
    df_grouped = df_grouped.drop(columns=["Consistent"], errors='ignore')

    return df_grouped, summary

def plot_histograms(df_before, df_after) -> plt.Figure:
    """
    Returns a matplotlib figure with histograms comparing CPM and counts before/after filtering.
    Use st.pyplot() to display this figure in Streamlit.
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    sns.set_style("whitegrid")

    # Filter out zero/negative values for log scale plots
    def filter_positive(data):
        return data[data > 0]

    # Controls
    ctrl_before = filter_positive(df_before["Control Count"])
    ctrl_after = filter_positive(df_after["Control Count"])
    if len(ctrl_before) > 0 and len(ctrl_after) > 0:
        sns.histplot(ctrl_before, bins=50, color="blue", alpha=0.5, label="Before", ax=axes[0, 0], log_scale=True)
        sns.histplot(ctrl_after, bins=50, color="red", alpha=0.5, label="After", ax=axes[0, 0], log_scale=True)
    axes[0, 0].set_title("Control Count Distribution")
    axes[0, 0].legend()

    # Experiments
    exp_before = filter_positive(df_before["Experiment Count"])
    exp_after = filter_positive(df_after["Experiment Count"])
    if len(exp_before) > 0 and len(exp_after) > 0:
        sns.histplot(exp_before, bins=50, color="blue", alpha=0.5, label="Before", ax=axes[0, 1], log_scale=True)
        sns.histplot(exp_after, bins=50, color="red", alpha=0.5, label="After", ax=axes[0, 1], log_scale=True)
    axes[0, 1].set_title("Experiment Count Distribution")
    axes[0, 1].legend()

    # CPM Control
    cpm_ctrl_before = filter_positive(df_before["Control CPM"])
    cpm_ctrl_after = filter_positive(df_after["Control CPM"])
    if len(cpm_ctrl_before) > 0 and len(cpm_ctrl_after) > 0:
        sns.histplot(cpm_ctrl_before, bins=50, color="blue", alpha=0.5, label="Before", ax=axes[1, 0], log_scale=True)
        sns.histplot(cpm_ctrl_after, bins=50, color="red", alpha=0.5, label="After", ax=axes[1, 0], log_scale=True)
    axes[1, 0].set_title("CPM Control Distribution")
    axes[1, 0].legend()

    # CPM Experiment
    cpm_exp_before = filter_positive(df_before["Experiment CPM"])
    cpm_exp_after = filter_positive(df_after["Experiment CPM"])
    if len(cpm_exp_before) > 0 and len(cpm_exp_after) > 0:
        sns.histplot(cpm_exp_before, bins=50, color="blue", alpha=0.5, label="Before", ax=axes[1, 1], log_scale=True)
        sns.histplot(cpm_exp_after, bins=50, color="red", alpha=0.5, label="After", ax=axes[1, 1], log_scale=True)
    axes[1, 1].set_title("CPM Experiment Distribution")
    axes[1, 1].legend()

    plt.tight_layout()
    return fig

# def filter_by_CPM(df: pd.DataFrame,conditions: dict, threshold_1=None, threshold_2=None, plot=False) -> tuple[pd.DataFrame, plt.Figure | None]:
#     """
#     Filters peptides based on CPM thresholds and returns filtered dataframe + histogram figure.
#     """

#     control_cols = [col for col, value in conditions.items() if value == "Control" and col in df.columns]
#     experiment_cols = [col for col, value in conditions.items() if value == "Experiment" and col in df.columns]

#     df["Control Count"] = df[control_cols].sum(axis=1)
#     df["Experiment Count"] = df[experiment_cols].sum(axis=1)

#     total_ctrl, total_exp = df["Control Count"].sum(), df["Experiment Count"].sum()
#     df["Control CPM"] = (df["Control Count"] / total_ctrl * 1e6) if total_ctrl > 0 else 0
#     df["Experiment CPM"] = (df["Experiment Count"] / total_exp * 1e6) if total_exp > 0 else 0

#     df_before = df.copy()
#     if threshold_1 is not None:
#         df = df[df["Experiment Count"] > threshold_1]
#     if threshold_2 is not None:
#         df = df[df["Experiment CPM"] > threshold_2]

#     df = df.reset_index(drop=True)
    

#     fig = plot_histograms(df_before, df) if plot else None
#     return df, fig

def filter_by_CPM(df: pd.DataFrame,
                  conditions: dict,
                  cpm_thresh: float,
                  min_samples: int,
                  plot: bool = False):
    """
    Filters peptides based on CPM (counts per million) thresholds.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input DataFrame with peptide counts
    conditions : dict
        Dictionary mapping column names to conditions (Control/Experiment)  
    cpm_thresh : float
        CPM threshold for filtering
    min_samples : int
        Minimum number of samples that must meet the threshold
    plot : bool
        Whether to generate comparison plots
        
    Returns:
    --------
    tuple: (filtered_df, figure or None)
    """
    
    # Input validation
    if df.empty:
        print("Warning: Input DataFrame is empty")
        return df.copy(), None
    
    if cpm_thresh < 0:
        raise ValueError("CPM threshold must be non-negative")
    
    if min_samples < 0:
        raise ValueError("Minimum samples must be non-negative")
    
    # Separate control and experiment columns
    ctrl_cols = [c for c, v in conditions.items() if v == "Control" and c in df.columns]
    exp_cols = [c for c, v in conditions.items() if v == "Experiment" and c in df.columns]
    
    if not exp_cols:
        print("Warning: No experiment columns found for CPM filtering")
        return df.copy(), None
    
    # Calculate library sizes (total counts per sample) with validation
    lib_sizes = df[exp_cols].sum()
    
    # Check for zero library sizes
    zero_lib_samples = lib_sizes[lib_sizes == 0].index.tolist()
    if zero_lib_samples:
        print(f"Warning: Zero library sizes found in samples: {zero_lib_samples}")
        # Remove zero library size samples from consideration
        exp_cols = [c for c in exp_cols if c not in zero_lib_samples]
        if not exp_cols:
            print("Warning: All experiment samples have zero counts")
            return df.copy(), None
        lib_sizes = lib_sizes[exp_cols]
    
    # Calculate CPM with safe division
    cpm = df[exp_cols].div(lib_sizes, axis=1) * 1e6
    
    # Apply filtering criteria
    keep = (cpm >= cpm_thresh).sum(axis=1) >= min_samples
    
    df_before = df.copy()
    df_after = df[keep].copy()
    
    # Calculate summary statistics for plotting
    if ctrl_cols:
        df_before["Control Count"] = df_before[ctrl_cols].sum(axis=1)
        ctrl_total_before = df_before["Control Count"].sum()
        df_before["Control CPM"] = (df_before["Control Count"] / ctrl_total_before * 1e6) if ctrl_total_before > 0 else 0
        
        df_after["Control Count"] = df_after[ctrl_cols].sum(axis=1) if not df_after.empty else pd.Series(dtype=float)
        ctrl_total_after = df_after["Control Count"].sum() if not df_after.empty else 0
        df_after["Control CPM"] = (df_after["Control Count"] / ctrl_total_after * 1e6) if ctrl_total_after > 0 else 0
    else:
        df_before["Control Count"] = 0
        df_before["Control CPM"] = 0
        df_after["Control Count"] = 0
        df_after["Control CPM"] = 0
    
    # Experiment counts and CPM
    df_before["Experiment Count"] = df_before[exp_cols].sum(axis=1)
    exp_total_before = df_before["Experiment Count"].sum()
    df_before["Experiment CPM"] = (df_before["Experiment Count"] / exp_total_before * 1e6) if exp_total_before > 0 else 0
    
    df_after["Experiment Count"] = df_after[exp_cols].sum(axis=1) if not df_after.empty else pd.Series(dtype=float)
    exp_total_after = df_after["Experiment Count"].sum() if not df_after.empty else 0  
    df_after["Experiment CPM"] = (df_after["Experiment Count"] / exp_total_after * 1e6) if exp_total_after > 0 else 0

    df_after.reset_index(drop=True, inplace=True)
    
    # Generate plot if requested and we have data
    fig = None
    if plot and not df_before.empty and not df_after.empty:
        try:
            fig = plot_histograms(df_before, df_after)
        except Exception as e:
            print(f"Warning: Could not generate comparison plot: {e}")
            fig = None

    print(f"CPM filtering: {len(df_before)} -> {len(df_after)} peptides retained")
    
    return df_after, fig

def filter_by_CPM_v2_style(
        df: pd.DataFrame,
        conditions: dict,
        threshold_count: int | None = None,
        threshold_cpm: float | None = None,
        output_path: str | None = None,
        plot: bool = False
) -> tuple[pd.DataFrame, Any | None]:
    """
    Filters peptides using version-2 logic (global Experiment Count + CPM_Exp),
    but follows the structure and defensive programming style of version 1.

    Parameters
    ----------
    df : pd.DataFrame
        Peptide-level count matrix.
    conditions : dict
        {column_name: "Control" | "Experiment"} mapping.
    threshold_count : int, optional
        Minimum aggregate Experiment Count required (version-2 `threshold_1`).
    threshold_cpm : float, optional
        Minimum aggregate CPM_Exp required (version-2 `threshold_2`).
    output_path : str, optional
        CSV destination for the filtered table.  If None, nothing is saved.
    plot : bool, default False
        Whether to return a before/after histogram figure.

    Returns
    -------
    (filtered_df, figure_or_None)
    """

    # ─────────────────────────── 1. Input validation ──────────────────────────
    if df.empty:
        print("Warning: input DataFrame is empty.")
        return df.copy(), None

    if threshold_count is not None and threshold_count < 0:
        raise ValueError("`threshold_count` must be non-negative.")

    if threshold_cpm is not None and threshold_cpm < 0:
        raise ValueError("`threshold_cpm` must be non-negative.")

    # ─────────────────────── 2. Column set-up via `conditions` ───────────────
    ctrl_cols = [c for c, v in conditions.items() if v == "Control"   and c in df.columns]
    exp_cols  = [c for c, v in conditions.items() if v == "Experiment" and c in df.columns]

    if not ctrl_cols or not exp_cols:
        raise ValueError("No valid control or experiment columns identified.")

    # ─────────────── 3. Aggregate counts & CPM (version-2 logic) ─────────────
    df_work = df.copy()                       # do **not** mutate caller’s df
    df_work["Control Count"]    = df_work[ctrl_cols].sum(axis=1)
    df_work["Experiment Count"] = df_work[exp_cols].sum(axis=1)

    total_ctrl = df_work["Control Count"].sum()
    total_exp  = df_work["Experiment Count"].sum()

    df_work["CPM_Control"] = (df_work["Control Count"]    / total_ctrl * 1e6) if total_ctrl else 0
    df_work["CPM_Exp"]     = (df_work["Experiment Count"] / total_exp  * 1e6) if total_exp  else 0

    # ───────────────────── 4. Store “before” snapshot for plots ──────────────
    df_before = df_work.copy()

    # ────────────────────────── 5. Apply dual thresholds ─────────────────────
    mask = pd.Series(True, index=df_work.index)
    if threshold_count is not None:
        mask &= df_work["Experiment Count"] > threshold_count
    if threshold_cpm   is not None:
        mask &= df_work["CPM_Exp"]         > threshold_cpm

    df_after = df_work.loc[mask].copy()
    df_after.reset_index(drop=True, inplace=True)

    # ──────────────────────── 6. Optional peptide cleaning ───────────────────
    if "Clean Peptide" in df_after.columns:
        # substitute stop codon (*) with Q, preserving other chars
        df_after["Clean Peptide"] = df_after["Clean Peptide"].str.replace(r"\*", "Q", regex=True)

    # ──────────────────────── 7. Optional disk output ────────────────────────
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        df_after.to_csv(output_path, index=False)

    # ───────────────────────── 8. Optional plotting ──────────────────────────
    fig = None
    if plot and not df_before.empty and not df_after.empty:
        try:
            fig = plot_histograms(df_before, df_after)   # ← your helper
        except Exception as e:
            print(f"Warning: could not generate histogram: {e}")

    # ───────────────────────────── 9. Reporting ──────────────────────────────
    print(f"Global-CPM filtering: {len(df_before)} → {len(df_after)} peptides retained")

    return df_after, fig

def draw_correlation_matrix(df: pd.DataFrame) -> plt.Figure:
    """
    Returns a matplotlib figure with pairwise correlations of CPM data.
    """
    cpm_df = df.filter(regex=r'^(Clean Peptide|CPM_.*)$')
    cpm_df.columns = ["peptide"] + [col.replace("CPM_", "") for col in cpm_df.columns[1:]]

    fig = sns.pairplot(cpm_df, diag_kind="hist", plot_kws={"alpha": 0.5, "s": 10}).fig
    fig.suptitle("CPM Pairwise Correlations", y=1.02)

    for ax in fig.axes:
        ax.set_xscale("log")
        ax.set_yscale("log")

    fig.tight_layout()
    return fig




def prepare_data(df: pd.DataFrame, columns_condition: dict):
    """
    Prepares count and metadata dataframes for DESeq2 analysis.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe containing peptide sequences and expression values (output of filter_by_CPM)
    columns_condition : dict
        Dictionary mapping column names to condition labels (Control/Experiment)
    
    Returns:
    --------
    tuple: (count_data, meta_data) 
        count_data: DataFrame with peptides as rows, samples as columns
        meta_data: DataFrame with sample metadata
    """
    
    # Input validation
    if df.empty:
        raise ValueError("Input DataFrame is empty")
    
    if "Clean Peptide" not in df.columns:
        raise ValueError("DataFrame must contain 'Clean Peptide' column")
    
    if not columns_condition:
        raise ValueError("columns_condition dictionary cannot be empty")
    
    # Validate that requested columns exist in the dataframe
    available_condition_cols = [col for col in df.columns if columns_condition.get(col) in ["Control", "Experiment"]]
    
    if not available_condition_cols:
        raise ValueError("No valid condition columns found in DataFrame")
    
    # Create working dataframe with only needed columns
    DE_df = df[["Clean Peptide"] + available_condition_cols].copy()
    
    # Check for duplicate peptides
    if DE_df["Clean Peptide"].duplicated().any():
        print("Warning: Duplicate peptides found. Aggregating by sum.")
        DE_df = DE_df.groupby("Clean Peptide", as_index=False).sum(numeric_only=True)
    
    # Set 'Clean Peptide' as index for proper alignment
    DE_df.set_index("Clean Peptide", inplace=True)
    
    # Remove peptides with all zero counts
    row_sums = DE_df.sum(axis=1)
    zero_rows = row_sums == 0
    if zero_rows.any():
        n_zero = zero_rows.sum()
        print(f"Warning: Removing {n_zero} peptides with zero counts across all samples")
        DE_df = DE_df[~zero_rows]
    
    if DE_df.empty:
        raise ValueError("No peptides with non-zero counts remain after filtering")
    
    # Create count data - only include columns that exist in both dataframe and conditions
    valid_columns = [col for col in available_condition_cols if col in DE_df.columns]
    count_data = DE_df[valid_columns].copy()
    
    # Ensure all values are non-negative integers
    if (count_data < 0).any().any():
        print("Warning: Negative values found in count data. Converting to absolute values.")
        count_data = count_data.abs()
    
    # Convert to integer type for DESeq2
    count_data = count_data.round().astype(int)
    
    # Create metadata dataframe
    meta_data = pd.DataFrame({
        "condition": [columns_condition[col] for col in valid_columns]
    }, index=valid_columns)
    
    # Convert condition to categorical for proper DESeq2 handling
    meta_data["condition"] = pd.Categorical(
        meta_data["condition"], 
        categories=["Control", "Experiment"], 
        ordered=False
    )
    
    # Final validation
    if count_data.columns.tolist() != meta_data.index.tolist():
        raise ValueError("Count data columns and metadata indices do not match!")
    
    # Check for sufficient samples in each condition
    condition_counts = meta_data["condition"].value_counts()
    for condition, count in condition_counts.items():
        if count < 2:
            print(f"Warning: Only {count} sample(s) in {condition} condition. "
                  "DESeq2 may have issues with statistical testing.")
    
    print(f"Prepared data: {len(count_data)} peptides × {len(count_data.columns)} samples")
    print(f"Conditions: {dict(condition_counts)}")
    
    return count_data, meta_data

def run_deseq2(
    count_data: pd.DataFrame,
    meta_data: pd.DataFrame,
    factor: str = "condition",
    control_label: str = "Control",
    treatment_label: str = "Experiment",
) -> pd.DataFrame:
    """
    Differential expression via PyDESeq2.

    Parameters
    ----------
    count_data : rows=peptides, cols=samples (integers)
    meta_data  : index=samples, one column `factor` with categories including control_label & treatment_label
    out_dir    : if non-empty, writes MA-plot there
    factor     : name of the column in meta_data to test
    control_label   : baseline level in factor
    treatment_label : level to compare against baseline

    Returns
    -------
    res_df : pd.DataFrame
        Indexed by peptide, with columns ['baseMean','log2FoldChange','pvalue','padj',...]
    """
# WARNInin: aqui creo que no esta comparando bien las cosas (Tiene que comparar por columnas)
    # 1) align metadata to samples
    meta = meta_data.reindex(count_data.columns)
    # 2) transpose for PyDESeq2
    counts = count_data.T

    # 3) build and run
    dds = DeseqDataSet(
        counts=counts,
        metadata=meta,
        design=f"~{factor}",
        n_cpus=1
    )
    dds.deseq2()  # size factors, dispersion, fit GLM

    # 4) stats with explicit contrast
    stats = DeseqStats(
        dds,
        contrast=[factor, treatment_label, control_label]
    )
    stats.summary()  # Wald test + shrinkage

    res = stats.results_df.copy()
    res.index = count_data.index  # restore peptides as rows

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(
        res["baseMean"],
        res["log2FoldChange"],
        c=(res["padj"] < 0.05).map({True: "red", False: "grey"}),
        alpha=0.6, s=20, edgecolors="none"
    )
    ax.set_xscale("log")
    ax.set_xlabel("Mean of normalized counts")
    ax.set_ylabel("Log2 fold change")
    ax.set_title(f"MA plot: {treatment_label} vs {control_label}")
    fig.tight_layout()


    return res, fig



def run_vst(
    count_data: pd.DataFrame,
    meta_data: pd.DataFrame,
    factor: str = "condition"
) -> pd.DataFrame:
    """
    Variance‐stabilizing transform via PyDESeq2, using size factors from .obsm.

    1) Fit DESeq2 (to compute size factors).
    2) Pull numeric size factors from dds.obsm["size_factors"].
    3) Normalize and log2-transform counts.

    Returns
    -------
    vst_df : pd.DataFrame
        Log2(normalized_counts + 1), same index/columns as count_data.
    """
    # 1) Align metadata & transpose for PyDESeq2
    meta = meta_data.copy().reindex(count_data.columns)
    counts = count_data.T

    # 2) Run the pipeline (this computes size factors)
    dds = DeseqDataSet(
        counts=counts,
        metadata=meta,
        design=f"~{factor}",
        size_factors_fit_type="ratio",
        n_cpus=1,
        quiet=True
    )
    dds.fit_size_factors()

    # 3) Extract numeric size factors from .obsm
    #    (an array of length n_samples)
    if "size_factors" not in dds.obs:
        raise KeyError("Expected 'size_factors' in dds.obs but not found.")
    sf_array = np.asarray(dds.obs["size_factors"])
    sf = pd.Series(sf_array, index=counts.index, name="size_factor")

    # 4) Divide and log2-transform
    normed = count_data.div(sf, axis=1)
    vst_df = np.log2(normed + 1)

    return vst_df



def significant_DE_peptides(df:pd.DataFrame, significance_threshold:float = 0.05, log2fc_threshold_experiment:float = 1.5):
    """
    Identifies differentially expressed (DE) peptides based on FDR (padj) and log2 fold change thresholds.
    Classifies peptides as upregulated, downregulated, or not significant.
    
    Parameters:
    df (pd.DataFrame): DataFrame containing differential expression results with 'padj' and 'log2FoldChange' columns.
    output_dir (str): Directory where output files will be saved.
    significance_threshold (float, optional): FDR threshold for significance. Default is 0.05.
    log2fc_threshold_experiment (float, optional): Absolute log2 fold change threshold for significance. Default is 1.5.
    
    Returns:
    pd.DataFrame: DataFrame with significant DE peptides, including an 'updown' column.
    
    Outputs:
    - Saves a CSV file of upregulated peptides in the specified output directory.
    - Prints summary statistics of DE peptide counts.
    """


    df["padj"] = df["padj"].astype(float)
    df["log2FoldChange"] = df["log2FoldChange"].astype(float)
    df["-log10(padj)"] = -np.log10(df["padj"] + 1e-10)
    df["-log10(pvalue)"] = -np.log10(df["pvalue"] + 1e-10)
    
    df.loc[:, "updown"] = df.apply(
        lambda row: "up" if (row["log2FoldChange"] > log2fc_threshold_experiment and row["padj"] < significance_threshold)
        else ("down" if (row["log2FoldChange"] < -log2fc_threshold_experiment and row["padj"] < significance_threshold) else "not_significant"),
        axis=1)
    

    # Count up/downregulated genes
    up_count = df[df["updown"] == "up"].shape[0]
    down_count = df[df["updown"] == "down"].shape[0]
    not_significant_count = df[df["updown"] == "not_significant"].shape[0]
    summary = {
        "Upregulated Peptides": up_count,
        "Downregulated Peptides": down_count,
        "Not Significant Peptides": not_significant_count,
        "Total Peptides": df.shape[0],
        "Significance Threshold": significance_threshold,
        "Log2FC Threshold": log2fc_threshold_experiment
    }
    return df, summary



def plot_heatmap(vst_data_df_filtered: pd.DataFrame) -> plt.Figure:
    """
    Returns a matplotlib figure showing a heatmap of VST-normalized counts.
    """
    fig, ax = plt.subplots(figsize=(12, 10))
    sns.heatmap(vst_data_df_filtered, cmap="viridis", annot=False, cbar_kws={'label': 'Expression (VST)'}, ax=ax)
    ax.set_title("Heatmap of VST Normalized Counts for Significant Peptides")
    ax.set_xlabel("Samples")
    ax.set_ylabel("Peptide Sequences")
    plt.xticks(rotation=45)
    plt.yticks(rotation=0)
    plt.tight_layout()
    return fig


def volcano_plot(
    df: pd.DataFrame,
    output_path: str = None,
    significance_threshold=0.05,
    log2fc_threshold=1.5,
    annotate_threshold=0.01,
    title="Volcano Plot"
) -> plt.Figure:
    """
    Creates a volcano plot and returns the matplotlib Figure.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain at least these columns: 'pvalue', 'padj', 'log2FoldChange'.
    output_path : str, optional
        If provided, saves the figure to this PDF path.
    significance_threshold : float or str
        FDR cutoff (e.g. 0.05). Will be cast to float.
    log2fc_threshold : float or str
        Absolute log₂ fold-change cutoff (e.g. 1.5). Will be cast to float.
    annotate_threshold : float or str
        Fraction of FDR used to pick which points to label. Will be cast to float.
    title : str
        Plot title.

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    # 1) Coerce thresholds to float
    significance_threshold = float(significance_threshold)
    log2fc_threshold    = float(log2fc_threshold)
    annotate_threshold  = float(annotate_threshold)

    # 2) Make a safe copy and coerce columns to numeric
    df = df.copy()
    for col in ("padj", "pvalue", "log2FoldChange"):
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # 3) Compute –log10(padj)
    df["neg_log10_padj"] = -np.log10(df["padj"] + 1e-10)

    # 4) Assign colors
    df["color"] = "grey"
    df.loc[df["log2FoldChange"] >  log2fc_threshold , "color"] = "red"
    df.loc[df["log2FoldChange"] < -log2fc_threshold, "color"] = "blue"
    df.loc[df["padj"] > significance_threshold, "color"] = "grey"

    # 5) Start plotting
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.scatterplot(
        x="log2FoldChange",
        y="neg_log10_padj",
        hue="color",
        palette={"grey": "grey", "red": "red", "blue": "blue"},
        data=df,
        alpha=0.6,
        s=50,
        edgecolor="none",
        ax=ax,
        legend=False
    )

    # 6) Reference lines
    ax.axvline( log2fc_threshold, color="black", linestyle="--", linewidth=1)
    ax.axvline(-log2fc_threshold, color="black", linestyle="--", linewidth=1)
    ax.axhline(-np.log10(significance_threshold), color="black", linestyle="--", linewidth=1)

    # 7) Annotate top up-regulated points
    mask_up = (
        (df["log2FoldChange"] >  log2fc_threshold) &
        (df["neg_log10_padj"]   > -np.log10(significance_threshold * annotate_threshold))
    )
    for idx, row in df.loc[mask_up].iterrows():
        ax.text(
            row["log2FoldChange"],
            row["neg_log10_padj"] + 0.1,
            str(idx),
            fontsize=8,
            ha="center"
        )
    
    # 8) Annotate top down-regulated points
    mask_down = (
        (df["log2FoldChange"] < -log2fc_threshold) &
        (df["neg_log10_padj"]   > -np.log10(significance_threshold * annotate_threshold))
    )
    for idx, row in df.loc[mask_down].iterrows():
        ax.text(
            row["log2FoldChange"],
            row["neg_log10_padj"] + 0.1,
            str(idx),
            fontsize=8,
            ha="center"
        )

    # 9) Legend and labels
    legend_elems = [
        Patch(facecolor="grey", label="Not significant"),
        Patch(facecolor="red",  label=f"Up (|FC|>{log2fc_threshold}, FDR>{significance_threshold})"),
        Patch(facecolor="blue", label=f"Down (|FC|<{log2fc_threshold}, FDR>{significance_threshold})"),
    ]
    ax.legend(handles=legend_elems, loc="upper left", frameon=True)

    ax.set_title(title, fontsize=14)
    ax.set_xlabel("Log2 Fold Change", fontsize=12)
    ax.set_ylabel("-Log10 Adjusted p-value", fontsize=12)
    ax.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()

    # 10) Save if requested
    if output_path:
        fig.savefig(output_path, dpi=300)

    return fig

# Function to generate WebLogo plots
def generate_weblogo(input_fasta, output_pdf, title):
    command = f"weblogo -f {input_fasta} -o {output_pdf} --title '{title}' --format pdf --errorbars NO --fineprint ''"
    subprocess.run(command, shell=True)

def to_fasta(peptides: pd.DataFrame, output_file: str = "peptides.fasta") -> str:
    """
    Converts a DataFrame of peptides into a FASTA file format.

    Parameters:
    peptides (pd.DataFrame): DataFrame containing peptide sequences.
    output_file (str): Path to the output FASTA file.

    Returns:
    str: Path to the generated FASTA file.
    """
    with open(output_file, "w") as f:
        for i, row in peptides.iterrows():
            f.write(f">peptide_{i}\n{row['variable_pep']}\n")
    return output_file


def parse_gibbscluster_output(dir,num_clusters):
    """
    Parses the GibbsCluster output file to extract peptide sequences for each cluster.

    Parameters:
    file_path (str): Path to the GibbsCluster output file.

    Returns:
    dict: Dictionary with clusters as keys and lists of peptides as values.
    """
    # Read the file, skipping the header row


    df = pd.read_csv(CURRENT_DIR+"/"+ dir +f"/res/gibbs.{num_clusters}g.ds.out", delim_whitespace=True, header=None, skiprows=1, 
                     names=["G", "Gn", "Num", "Sequence", "Core", "o", "of", "ip", "IP", "il", "IL", "dp", "DP", "dl", "DL", 
                            "Annotation", "sS", "Self", "bgG", "bgG_Val", "bgS", "bgS_Val", "cS", "cScore"])
    df["Gn"]= df["Gn"] + 1
    return df

def merge_data(cluster_df: pd.DataFrame, value_df: pd.DataFrame, value_col) -> pd.DataFrame:
    """
    Merges a dataframe containing peptide sequences and clusters with a specified value column.
    
    Parameters:
        cluster_df (pd.DataFrame): DataFrame containing peptide sequences and cluster assignments.
        value_df (pd.DataFrame): DataFrame containing peptide sequences and corresponding values.
        value_col (str or list): The name of the column(s) in value_df to merge on.
    
    Returns:
        pd.DataFrame: Merged DataFrame including the specified value column.
    """
    if isinstance(value_col, str):
        value_col = [value_col]
    
    # Validate that all requested columns exist
    missing_cols = [col for col in value_col if col not in value_df.columns]
    if missing_cols:
        raise ValueError(f"Columns not found in value_df: {missing_cols}")
    
    # Create subset with proper column selection
    value_subset = value_df[["Clean Peptide"] + value_col].rename(columns={"Clean Peptide": "Sequence"})
    
    # Perform merge with validation
    result = pd.merge(cluster_df, value_subset, on="Sequence", how="left")
    
    # Check for unmatched sequences
    unmatched = result["Sequence"][result[value_col[0]].isna()]
    if len(unmatched) > 0:
        print(f"Warning: {len(unmatched)} sequences from cluster_df not found in value_df")
    
    return result

def box_plot(df: pd.DataFrame, x_col: str, y_col: str, title: str, output_dir: str = "") -> plt.Figure:
    """
    Creates a boxplot with a strip plot overlay and returns the figure.
    
    Parameters:
        df (pd.DataFrame): DataFrame containing data to plot.
        x_col (str): Column name for the x-axis.
        y_col (str): Column name for the y-axis.
        title (str): Title of the plot.
        output_dir (str, optional): Directory path to save the plot. Defaults to "" (no saving).
        
    Returns:
        plt.Figure: The matplotlib figure object
    """
    
    # Input validation
    if df.empty:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "No data to plot", ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title)
        return fig
    
    if x_col not in df.columns:
        raise ValueError(f"Column '{x_col}' not found in DataFrame")
    if y_col not in df.columns:
        raise ValueError(f"Column '{y_col}' not found in DataFrame")
    
    # Remove rows with NaN values in plotting columns
    plot_df = df[[x_col, y_col]].dropna()
    
    if plot_df.empty:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, f"No valid data for columns {x_col}, {y_col}", 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title)
        return fig
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    try:
        # Create boxplot
        sns.boxplot(x=x_col, y=y_col, data=plot_df, showfliers=False, palette="Set2",
                    boxprops=dict(facecolor="none", edgecolor="black"),
                    whiskerprops=dict(color="black"), medianprops=dict(color="red"), ax=ax)
        
        # Add strip plot overlay
        sns.stripplot(x=x_col, y=y_col, data=plot_df, color="blue", size=2, alpha=0.5, jitter=True, ax=ax)
        
    except Exception as e:
        ax.text(0.5, 0.5, f"Error creating plot: {str(e)}", 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title)
        return fig
    
    ax.set_xlabel("Cluster")
    ax.set_ylabel(y_col)
    ax.set_title(title)
    plt.setp(ax.get_xticklabels(), rotation=45)
    ax.grid(axis="y", linestyle="--", alpha=0.7)
    
    plt.tight_layout()
    
    if output_dir:
        try:
            fig.savefig(output_dir, dpi=300, bbox_inches='tight')
        except Exception as e:
            print(f"Warning: Could not save plot to {output_dir}: {e}")
    
    return fig


# Function to remove constant positions (keep only variable positions)
def remove_constant(peptide, constant_pos):
    return "".join(aa for index, aa in enumerate(peptide) if index not in constant_pos)

# Function to save sequences as FASTA format
def save_as_fasta(peptides, filename):
    with open(filename, "w") as f:
        for i, seq in enumerate(peptides):
            f.write(f">seq{i}\n{seq}\n")

def top_n_per_cluster(df:pd.DataFrame,n, col="log2FoldChange"):
    top_pepts = {}
    for i in df["Gn"].unique():
        ith_cluster_df = df.loc[df["Gn"] == i].sort_values([col], ascending=False).reset_index(drop=True)
        n_peptides = round(len(ith_cluster_df.index) * n/100)
        top_pepts[i] = ith_cluster_df.head(n_peptides)["Sequence"].values
    return top_pepts

def constant_positions(pattern: str):
    """
    Given a regex like "A.C.{7}C", returns the 0-based positions of
    the constant letters once the pattern is fully expanded.
    E.g. "A.C.{7}C" -> [0, 2, 10].

    Supported syntax:
      - Literals: A–Z, a–z, 0–9, or escaped characters (e.g. '\\.')
      - Wildcard: .
      - Fixed quantifier: {n} immediately after a literal or .

    Anything else (groups, +, *, ?, alternation, character classes, etc.)
    will raise a ValueError.
    """
    token_re = re.compile(r"""
        (\\.|\.|[A-Za-z0-9])   # either escaped char (\x), dot, or alphanumeric
        (?:\{(\d+)\})?         # optional fixed-quantifier {n}
    """, re.VERBOSE)

    idx = 0
    const_pos = []
    pos = 0
    while pos < len(pattern):
        m = token_re.match(pattern, pos)
        if not m:
            raise ValueError(f"Unsupported regex token at position {pos}: '{pattern[pos:]}'")
        token, quant = m.group(1), m.group(2)
        count = int(quant) if quant is not None else 1

        # Is this a literal?
        if token.startswith("\\"):
            # escaped literal
            for i in range(count):
                const_pos.append(idx + i)
        elif token == ".":
            # wildcard: skip
            pass
        else:
            # plain literal A–Z or 0–9
            for i in range(count):
                const_pos.append(idx + i)

        idx += count
        pos = m.end()

    return const_pos
