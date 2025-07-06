import streamlit as st
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from Library.lib import *  # your own library functions

st.set_page_config(page_title="Peptide Analysis App", layout="wide")
st.title("🧬 Peptide Differential Expression & Clustering")

# Upload CSV file
uploaded_file = st.file_uploader("Upload your peptide CSV file", type=["csv"])
if uploaded_file is None:
    st.info("👈 Upload a CSV file to start.")
    st.stop()

# Load data
data = pd.read_csv(uploaded_file)
st.success(f"Loaded data with shape {data.shape}")

# Rename columns
st.header("🔠 Column Renaming")
rename_dict = {}
for col in data.columns:
    new_name = st.text_input(f"Rename column '{col}'", value=col)
    rename_dict[col] = new_name
data = data.rename(columns=rename_dict)

st.write("Renamed Data Preview:", data.head())

# Parameters for grouping/filtering
st.header("⚙️ Analysis Parameters")

regex = st.text_input("Peptide regex for grouping", value=r'A.C.{7}C')
st.write("Example regex: `A.C.{7}C` matches peptides starting with A, followed by any amino acid, then C, and 7 more amino acids ending with C.")
cpm_threshold = st.number_input("CPM filter threshold", min_value=0, value=5)
min_count = st.number_input("Minimum samples with counts >= CPM threshold", min_value=1, value=1)
plot_filtering = st.checkbox("Plot CPM filtering", value=True)

# Run clustering
num_clusters = st.number_input("Number of Clusters", min_value=2, value=4)

# Define experimental design
st.header("🧪 Experimental Conditions")
columns_conditions = {}
for col in data.columns[1:]:  # skip peptide ID column
    condition = st.selectbox(f"Condition for {col}", options=["Control", "Experiment", "Exclude"], index=0)
    if condition != "Exclude":
        columns_conditions[col] = condition

# Run analysis button
if st.button("🚀 Run Analysis"):
    with st.spinner("Running analysis..."):

        # Group by peptide
        grouped = group_by_peptide(data, len(columns_conditions), regex, "Cleaned_Peptides.csv")

        # Filter by CPM
        filtered, fig = filter_by_CPM(grouped , cpm_threshold, min_count, plot=plot_filtering)

        st.pyplot(fig)

        # Prepare data for DESeq2
        count_data, meta_data = prepare_data(filtered, columns_conditions)

        # Differential expression
        res_df, fig = run_deseq2(count_data, meta_data)
        st.pyplot(fig)
        res_df = significant_DE_peptides(res_df)

        # Adjust p-values
        _, padj, _, _ = multipletests(res_df["pvalue"], method="fdr_bh")
        res_df["padj"] = padj
        res_df["-log10(padj)"] = -np.log10(padj + 1e-10)
        res_df = significant_DE_peptides(res_df)

        st.subheader("DE Results Table")
        st.dataframe(res_df.head(40))

        # Run VST for heatmap
        vst_data = run_vst(count_data, meta_data)

        upregulated = res_df.loc[res_df["updown"] == "up"].index
        vst_up = vst_data.loc[upregulated]
        st.subheader("🔺 Heatmap of Upregulated Peptides")
        st.pyplot(plot_heatmap(vst_up))  # assuming plot_heatmap shows or returns a fig

        # Volcano plot
        st.subheader("📈 Volcano Plot of DE Peptides")
        st.pyplot(volcano_plot(res_df))

        # Sequence processing
        const_p = constant_positions(regex)
        st.text("Debug: Constant positions in peptides: " + str(const_p))
        res_df["variable_pep"] = res_df.index.map(lambda x: remove_constant(x, const_p))
        res_df["additional_cys"] = res_df["variable_pep"].apply(lambda pep: "C" in pep)
        adj_DE = res_df.loc[~res_df["additional_cys"]].drop(columns=["additional_cys"]).reset_index()

        # # Save peptide lists
        # adj_DE[adj_DE["updown"] == "up"][["variable_pep"]]
        # adj_DE[adj_DE["updown"] == "up"][["Clean Peptide"]]

        # Generate WebLogo plots
        weblogos = [
            ("upregulated.fasta", "upregulated.pdf", "Enriched - Original"),
            ("upregulated_variable.fasta", "upregulated_variable.pdf", "Enriched - Trimmed"),
        ]
        for input_fasta, output_pdf, title in weblogos:
            generate_weblogo(input_fasta, output_pdf, title)

        st.success("WebLogos generated and saved!")

        
        # button to run clustering
        st.subheader("🔍 Clustering Analysis")
        st.write("Running Gibbs clustering on upregulated peptides...")


        clustering_dir = gibbs_cluster(
            "Upregulated_peptides.pep", 11, num_clusters,
            "Clustering", "Clustering_Results", logos=True
        )
        clusters = parse_gibbscluster_output(f"{clustering_dir}/res/gibbs.{num_clusters}g.ds.out")
        adj_DE.rename(columns={"Clean Peptide": "Sequence"}, inplace=True)
        clusters = pd.merge(clusters, adj_DE[["Sequence", "log2FoldChange", "padj", "-log10(padj)"]])

        # Display cluster boxplots
        st.subheader("📊 Cluster Analysis")
        st.pyplot(box_plot(clusters, "Gn", "log2FoldChange", "log2FC per Cluster"))
        st.pyplot(box_plot(clusters, "Gn", "-log10(padj)", "-log10(padj) per Cluster"))
        st.pyplot(box_plot(clusters, "Gn", "padj", "padj per Cluster"))

        st.success("Analysis complete ✅")

