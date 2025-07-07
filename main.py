import streamlit as st
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from Library.lib import *  # your own library functions

st.set_page_config(page_title="Peptide Analysis App", layout="wide")
st.title("Peptide Differential Expression & Clustering")

# Upload CSV file
uploaded_file = st.file_uploader("Upload your peptide CSV file", type=["csv"])
if uploaded_file is None:
    st.info("Upload a CSV file to start.")
    st.stop()

# Load data
data = pd.read_csv(uploaded_file)
st.success(f"Loaded data with shape {data.shape}")


st.write("Renamed Data Preview:", data.head())

# Parameters for grouping/filtering
st.header("Analysis Parameters")

regex = st.text_input("Peptide library form (it needs to be in the form of a python regular expression)", value=r'A.C.{7}C')
st.write("Example regex: `A.C.{7}C` matches peptides starting with A (for alanine), followed by any natural amino acid, then C, and 7 more amino acids ending with C.")
cpm_threshold = st.number_input("CPM (counts per million) filter threshold", min_value=0, value=5)
min_count = st.number_input("Minimum samples with counts >= CPM threshold", min_value=1, value=1)
plot_filtering = st.checkbox("Plot CPM filtering", value=True)

# Run clustering
num_clusters = st.number_input("Number of Clusters", min_value=2, value=4)

# Define experimental design
st.header("Experimental Conditions")
st.write("Define the experimental conditions for each peptide column (if it is either the protein sample or the control sample). Select 'Exclude' to skip a column from analysis.")
columns_conditions = {}
for col in data.columns[1:]:  # skip peptide ID column
    condition = st.selectbox(f"Condition for {col}", options=["Control", "Experiment", "Exclude"], index=0)
    if condition != "Exclude":
        columns_conditions[col] = condition

# Run analysis button
if st.button("🚀 Run Analysis"):
    with st.spinner("Running analysis..."):

        # Group by peptide
        grouped, summary = group_by_peptide(data, columns_conditions, regex)

        
        # Display grouped data
        st.subheader("Grouped Data")
        st.dataframe(grouped.head(40))

        st.subheader("Grouped Data Summary")
        for key, value in summary.items():
            st.write(f"{key}: {value}")


        # Filter by CPM
        filtered, fig = filter_by_CPM(grouped, columns_conditions , cpm_threshold, min_count, plot=plot_filtering)

        if plot_filtering:
            st.subheader("CPM Filtering Histogram")
            st.write("Histogram of CPM values before and after filtering")
            st.pyplot(fig)

        st.subheader("Filtered Data")
        st.dataframe(filtered.head(40))

        # display the number of peptides before and after filtering
        st.write(f"Number of peptides before filtering: {grouped.shape[0]}")
        st.write(f"Number of peptides after filtering: {filtered.shape[0]}")

        # Prepare data for DESeq2
        count_data, meta_data = prepare_data(filtered, columns_conditions)

        # Differential expression
        res_df, fig = run_deseq2(count_data, meta_data)
        st.subheader("Differential Expression Results")
        st.write("Results of differential expression analysis using DESeq2")
        st.pyplot(fig)

        res_df, summary = significant_DE_peptides(res_df)

        st.subheader("Significant DE Peptides")
        st.write("Peptides with significant differential expression (FDR < 0.05)")

        for key, value in summary.items():
            st.write(f"{key}: {value}")

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

        st.subheader("📑 Upregulated Peptides")
        st.write("List of upregulated peptides after filtering and processing:")
        st.dataframe(adj_DE[adj_DE["updown"] == "up"][["variable_pep", "log2FoldChange", "padj"]])



        # turn the peptides into a fasta file

        upregulated_variable_peps = adj_DE[adj_DE["updown"] == "up"][["variable_pep"]]
        
        # button to run clustering
        st.subheader("🔍 Clustering Analysis")
        st.write("Running Gibbs clustering on upregulated peptides...")


        clustering_dir = gibbs_cluster(
            upregulated_variable_peps , 11, num_clusters,
            "Clustering", "Clustering_Results", logos=True
        )
        clusters = parse_gibbscluster_output(clustering_dir,num_clusters=num_clusters)
        adj_DE.rename(columns={"Clean Peptide": "Sequence"}, inplace=True)
        clusters = pd.merge(clusters, adj_DE[["Sequence", "log2FoldChange", "padj", "-log10(padj)"]])


        # Generate WebLogo plots
        weblogos = [
            ("upregulated.fasta", "upregulated.pdf", "Enriched - Original"),
            ("upregulated_variable.fasta", "upregulated_variable.pdf", "Enriched - Trimmed"),
        ]
        for input_fasta, output_pdf, title in weblogos:
            generate_weblogo(input_fasta, output_pdf, title)

        st.success("WebLogos generated and saved!")


        # Display cluster boxplots
        st.subheader("📊 Cluster Analysis")
        st.pyplot(box_plot(clusters, "Gn", "log2FoldChange", "log2FC per Cluster"))
        st.pyplot(box_plot(clusters, "Gn", "-log10(padj)", "-log10(padj) per Cluster"))
        st.pyplot(box_plot(clusters, "Gn", "padj", "padj per Cluster"))

        st.success("Analysis complete ✅")

