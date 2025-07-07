
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import subprocess
import re
from matplotlib.patches import Patch
# from rpy2.robjects import pandas2ri
# import rpy2.robjects as ro
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


# get the rute from the current file
CURRENT_DIR = os.path.join("/".join(os.path.dirname(__file__).split("/")[:-1]))

GIBS_CLUSTER_BIN = os.path.join(CURRENT_DIR, "gibbscluster-2.0/gibbscluster")



def group_by_peptide(df:pd.DataFrame, conditions, pattern):
    """
    Groups peptides based on a given pattern(the library peptide form) and aggregates their counts.

    This function processes a DataFrame containing peptide sequences and their respective counts
    across multiple samples. It calculates the total peptide count, filters peptides matching a
    specific pattern, and aggregates them based on a cleaned version of the peptide sequence (without the GSG linker).

    Parameters:
    -----------
    df : pandas.DataFrame
        A DataFrame where each row represents a peptide with its associated counts. Merged dataframe from the Data_Merge python script.
    n_samples : int
        The number of sample columns to consider for aggregation, normally 6, 3 replicates of the protein panning and three for the control panning.
    pattern : str
        A regular expression pattern used to filter and extract peptides: a regular expression matching the peptide library.
    out_dir : str
        The output directory where the grouped peptide data will be saved as a CSV file.

    Returns:
    --------
    pandas.DataFrame
        A DataFrame grouped by the cleaned peptide sequences, with counts aggregated and sorted.
        The valid peptides that have the form of the library. 
    """

    # Add a column for the total count across specific columns

    # WARNNIG: THE PEPTIDE COLUMN MUST BE NAMED 'peptide' IN THE INPUT DATAFRAME
    df = df[["peptide"] + list(conditions.keys())].copy()
    df["Total Count"] = df[conditions.keys()].sum(axis=1)

    # Calculate total peptides and valid peptides
    total_peptides = df["Total Count"].sum()
    df["Consistent"] = df["peptide"].astype(str).str.contains(pattern)
    valid_peptides = df.loc[df["Consistent"], "Total Count"].sum()

    # Print results

    
    summary = {
        "Total Peptides": total_peptides,
        "Valid Peptides": valid_peptides,
        "Success Rate": f"{100 * valid_peptides / total_peptides:.2f}%"
    }


    # Extract the clean peptide matching the pattern
    df["Clean Peptide"] = df['peptide'].astype(str).str.extract(f"({pattern})", expand=False)
    df["Clean Peptide"] = df['Clean Peptide'].astype(str).str.replace("*", "Q")

    # Aggregate by clean peptide and sort
    df_grouped = df.groupby("Clean Peptide", as_index=False).sum(numeric_only=True)
    df_grouped = df_grouped.sort_values("Total Count", ascending=False)

    df_grouped.reset_index(inplace=True, drop=True)

    # Drop unnecessary column
    df_grouped = df_grouped.drop(columns=["Consistent"], errors='ignore')

    return df_grouped, summary

def plot_histograms(df_before, df_after) -> plt.Figure:
    """
    Returns a matplotlib figure with histograms comparing CPM and counts before/after filtering.
    Use st.pyplot() to display this figure in Streamlit.
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    sns.set_style("whitegrid")

    # Controls
    sns.histplot(df_before["Control Count"], bins=50, color="blue", alpha=0.5, label="Before", ax=axes[0, 0], log_scale=True)
    sns.histplot(df_after["Control Count"], bins=50, color="red", alpha=0.5, label="After", ax=axes[0, 0], log_scale=True)
    axes[0, 0].set_title("Control Count Distribution")
    axes[0, 0].legend()

    # Experiments
    sns.histplot(df_before["Experiment Count"], bins=50, color="blue", alpha=0.5, label="Before", ax=axes[0, 1], log_scale=True)
    sns.histplot(df_after["Experiment Count"], bins=50, color="red", alpha=0.5, label="After", ax=axes[0, 1], log_scale=True)
    axes[0, 1].set_title("Experiment Count Distribution")
    axes[0, 1].legend()

    # CPM Control
    sns.histplot(df_before["Control CPM"], bins=50, color="blue", alpha=0.5, label="Before", ax=axes[1, 0], log_scale=True)
    sns.histplot(df_after["Control CPM"], bins=50, color="red", alpha=0.5, label="After", ax=axes[1, 0], log_scale=True)
    axes[1, 0].set_title("CPM Control Distribution")
    axes[1, 0].legend()

    # CPM Experiment
    sns.histplot(df_before["Experiment CPM"], bins=50, color="blue", alpha=0.5, label="Before", ax=axes[1, 1], log_scale=True)
    sns.histplot(df_after["Experiment CPM"], bins=50, color="red", alpha=0.5, label="After", ax=axes[1, 1], log_scale=True)
    axes[1, 1].set_title("CPM Experiment Distribution")
    axes[1, 1].legend()

    plt.tight_layout()
    return fig

def filter_by_CPM(df: pd.DataFrame,conditions: dict, threshold_1=None, threshold_2=None, plot=False) -> tuple[pd.DataFrame, plt.Figure | None]:
    """
    Filters peptides based on CPM thresholds and returns filtered dataframe + histogram figure.
    """

    control_cols = [col for col, value in conditions.items() if value == "Control" and col in df.columns]
    experiment_cols = [col for col, value in conditions.items() if value == "Experiment" and col in df.columns]

    df["Control Count"] = df[control_cols].sum(axis=1)
    df["Experiment Count"] = df[experiment_cols].sum(axis=1)

    total_ctrl, total_exp = df["Control Count"].sum(), df["Experiment Count"].sum()
    df["Control CPM"] = (df["Control Count"] / total_ctrl * 1e6) if total_ctrl > 0 else 0
    df["Experiment CPM"] = (df["Experiment Count"] / total_exp * 1e6) if total_exp > 0 else 0

    df_before = df.copy()
    if threshold_1 is not None:
        df = df[df["Experiment Count"] > threshold_1]
    if threshold_2 is not None:
        df = df[df["Experiment CPM"] > threshold_2]

    df = df.reset_index(drop=True)
    

    fig = plot_histograms(df_before, df) if plot else None
    return df, fig



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




def prepare_data(df:pd.DataFrame, columns_condition:dict):
    """
    Prepares count and metadata dataframes for DESeq2 analysis. Use as the input dataframe the output of the filter_by_CPM function.
    
    Parameters:
    df (pd.DataFrame): The input dataframe containing peptide sequences and expression values.
    
    Returns:
    tuple: (count_data, meta_data) where count_data contains expression values and meta_data contains sample conditions.
    """
    DE_df = df[["Clean Peptide"] + 
                          [col for col in df.columns if columns_condition.get(col) in ["Control", "Experiment"]]].copy()
    # Set 'Clean Peptide' as index for proper alignment
    DE_df.set_index("Clean Peptide", inplace=True)

    count_data = DE_df[columns_condition.keys()].copy()


    meta_data = pd.DataFrame({"condition": columns_condition.values()},
                             index=count_data.columns)
    meta_data["condition"] = pd.Categorical(meta_data["condition"], categories=["Control", "Experiment"], ordered=False)
    assert count_data.columns.tolist() == meta_data.index.tolist(), "Indices of count_data and meta_data do not match!"
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
    dds.deseq2()

    # 3) Extract numeric size factors from .obsm
    #    (an array of length n_samples)
    if "size_factors" not in dds.obsm:
        raise KeyError("Expected 'size_factors' in dds.obsm but not found.")
    sf_array = np.asarray(dds.obsm["size_factors"])
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
    mask = (
        (df["log2FoldChange"] >  log2fc_threshold) &
        (df["neg_log10_padj"]   > -np.log10(significance_threshold * annotate_threshold))
    )
    for idx, row in df.loc[mask].iterrows():
        ax.text(
            row["log2FoldChange"],
            row["neg_log10_padj"] + 0.1,
            str(idx),
            fontsize=8,
            ha="center"
        )

    # 8) Legend and labels
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

    # 9) Save if requested
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

def gibbs_cluster(input_file,motif_length,num_clusters, output_dir, run_name, logos=True):

    """
    Runs GibbsClustering on a given peptide sequence input file.

    Parameters:
    input_file (str): Path to the input FASTA file containing peptide sequences.
    output_dir (str): Directory where GibbsCluster output will be saved.
    num_clusters (int): Number of clusters to generate.
    run_name (str): Name for the GibbsCluster run.

    Returns:
    str: The path of the clustering files.
    """
    def get_session_id(log_string):

        """
        Extracts the session ID from the GibbsCluster log output.

        Parameters:
        log_string (str): The standard output from GibbsCluster execution.

        Returns:
        str or None: Extracted session ID if found, otherwise None.
        """
        
        match = re.search(r"#Session ID:\s*(\d+)", log_string)
        return match.group(1) if match else None
    
    # save the input file as a fasta file
    if isinstance(input_file, pd.DataFrame):
        input_file = to_fasta(input_file, os.path.join(output_dir, "input_peptides.fasta"))
            
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    input_file_dir = os.path.abspath(input_file)
    print(f"Input file for GibbsCluster: {input_file_dir}")

    
    command = [
        GIBS_CLUSTER_BIN,
        "-f", input_file_dir,
        "-R", output_dir,
        "-P", run_name, 
        "-g", str(num_clusters),
        "-l", str(motif_length),
        "-s8", 
        "-b0.8",
        "-q0",
        "-T",
        "-j2",
        "-S10",
        "-c0",
        "-z0",

    ]
    
    try:
        res = subprocess.run(command, check=True, text=True, capture_output=True)
        print(f"GibbsCluster run '{run_name}' completed successfully.")
        out_path = output_dir+"/" + run_name + "_" + get_session_id(res.stdout)
        print(f"Find the results in :{ out_path}")
        if logos:
            os.makedirs(os.path.join(out_path,"logos"), exist_ok=True)
            for i in range(1, int(num_clusters) + 1):  
                cluster_file = os.path.join(out_path,"cores", f"gibbs.{i}of{num_clusters}.core")  # Adjust filename based on GibbsCluster output
                output_image = os.path.join(out_path,"logos", f"cluster_{i}.pdf")
                print(cluster_file)
                if os.path.exists(cluster_file):  # Ensure the file exists
                    generate_weblogo(cluster_file, output_image, f"cluster {i}")
                    print(f"Generated WebLogo for Cluster {i}: {output_image}")
        else:
            print(f"Cluster file {cluster_file} not found.")
        return out_path
    except subprocess.CalledProcessError as e:
        print(f"Error running GibbsCluster: {e.stderr}")
        return None

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
        value_col (str): The name of the column in value_df to merge on.
    
    Returns:
        pd.DataFrame: Merged DataFrame including the specified value column.
    """
    if type(value_col) == str:
        value_col = [value_col]

    value_df = value_df[["Clean Peptide"] + [value_col]].rename(columns={"Clean Peptide": "Sequence"})
    return pd.merge(cluster_df, value_df, on="Sequence")

def box_plot(df: pd.DataFrame, x_col: str, y_col: str, title: str, output_dir: str = "") -> None:
    """
    Creates a boxplot with a strip plot overlay.
    
    Parameters:
        df (pd.DataFrame): DataFrame containing data to plot.
        x_col (str): Column name for the x-axis.
        y_col (str): Column name for the y-axis.
        title (str): Title of the plot.
        output_dir (str, optional): Directory path to save the plot. Defaults to "" (no saving).
    """
    plt.figure(figsize=(8, 6))
    sns.boxplot(x=x_col, y=y_col, data=df, showfliers=False, palette="Set2",
                boxprops=dict(facecolor="none", edgecolor="black"),
                whiskerprops=dict(color="black"), medianprops=dict(color="red"))
    sns.stripplot(x=x_col, y=y_col, data=df, color="blue", size=2, alpha=0.5, jitter=True)
    
    plt.xlabel("Cluster")
    plt.ylabel(y_col)
    plt.title(title)
    plt.xticks(rotation=45)
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    
    if output_dir:
        plt.savefig(output_dir)
    plt.show()


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
