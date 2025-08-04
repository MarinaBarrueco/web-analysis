import streamlit as st
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from Library.lib import *  # your own library functions
from Library.Clustering import gibbs_cluster
from Library.GibbsClusterAdvanced import gibbs_cluster_advanced, GibbsClusterAdvanced
from Library.validation import create_validation_summary, validate_clustering_parameters
from Library.enhanced_gibbs_clustering import enhanced_gibbs_cluster_with_validation

st.set_page_config(page_title="Peptide Analysis App", layout="wide")
st.title("🧬 Peptide Differential Expression & Clustering")

# Add helpful introduction
st.markdown("""
**Welcome to the Peptide Analysis Pipeline!** 

This tool performs comprehensive differential expression analysis and clustering of peptide sequences.

📋 **Required data format:**
- CSV file with peptide sequences and count data
- First column should contain peptide sequences (named 'peptide')
- Additional columns should contain count data for different samples
""")

with st.expander("📖 Analysis Steps Overview"):
    st.markdown("""
    1. **Data Upload & Validation** - Load your CSV file and validate format
    2. **Parameter Configuration** - Set regex pattern and filtering thresholds  
    3. **Experimental Design** - Define control vs experimental conditions
    4. **Differential Expression** - Statistical analysis using DESeq2
    5. **Clustering Analysis** - Gibbs clustering of significant peptides
    6. **Results Visualization** - Heatmaps, volcano plots, and sequence logos
    """)

# Upload CSV file
uploaded_file = st.file_uploader("Upload your peptide CSV file", type=["csv"])
if uploaded_file is None:
    st.info("Upload a CSV file to start.")
    st.stop()

# Load data
try:
    data = pd.read_csv(uploaded_file)
    st.success(f"Loaded data with shape {data.shape}")
    st.write("Data Preview:", data.head())
    
    # Parameters for grouping/filtering
    st.header("Analysis Parameters")

    regex = st.text_input("Peptide library form (it needs to be in the form of a python regular expression)", value=r'A.C.{7}C')
    st.write("Example regex: `A.C.{7}C` matches peptides starting with A (for alanine), followed by any natural amino acid, then C, and 7 more amino acids ending with C.")
    cpm_threshold = st.number_input("CPM (counts per million) filter threshold", min_value=0, value=5)
    min_count = st.number_input("Minimum samples with counts >= CPM threshold", min_value=1, value=1)
    plot_filtering = st.checkbox("Plot CPM filtering", value=True)

    # Run clustering
    st.header("Clustering Parameters")
    use_advanced_clustering = st.checkbox("Use Advanced Gibbs Clustering (more features)", value=True)
    
    # Consensus validation options
    st.subheader("🔬 Consensus Validation (Recommended)")
    use_consensus_validation = st.checkbox("Enable Consensus Clustering Validation", value=True, 
                                         help="Provides robust cluster number selection and quality assessment")
    
    if use_consensus_validation:
        col_val1, col_val2 = st.columns(2)
        with col_val1:
            n_consensus_iterations = st.number_input("Consensus Iterations", min_value=20, max_value=200, value=50,
                                                   help="Number of subsampling iterations for consensus matrix")
            enable_bootstrap = st.checkbox("Enable Bootstrap Validation", value=True,
                                         help="Additional stability assessment (slower but more thorough)")
        with col_val2:
            n_bootstrap_iterations = st.number_input("Bootstrap Iterations", min_value=10, max_value=100, value=30,
                                                   help="Number of bootstrap iterations for stability assessment")
            auto_k_selection = st.checkbox("Automatic K Selection", value=True,
                                         help="Automatically select optimal number of clusters based on validation metrics")
    else:
        n_consensus_iterations = 0
        n_bootstrap_iterations = 0
        enable_bootstrap = False
        auto_k_selection = False

    if use_advanced_clustering:
        col1, col2 = st.columns(2)
        with col1:
            cluster_range = st.selectbox("Cluster Configuration", 
                                       ["Single number", "Range"], 
                                       index=0)
            if cluster_range == "Single number":
                num_clusters = st.number_input("Number of Clusters", min_value=2, value=4)
            else:
                min_clusters = st.number_input("Minimum Clusters", min_value=1, value=2)
                max_clusters = st.number_input("Maximum Clusters", min_value=2, value=5)
                num_clusters = (min_clusters, max_clusters)
            
            num_seeds = st.number_input("Number of Seeds", min_value=1, value=3, 
                                      help="Multiple seeds for better convergence")
            motif_length = st.number_input("Motif Length", min_value=6, value=8)
            
        with col2:
            iterations = st.number_input("Iterations per Temperature Step", min_value=5, value=10)
            temperature_start = st.number_input("Starting Temperature", min_value=0.1, value=1.5)
            temperature_steps = st.number_input("Temperature Steps", min_value=5, value=20)
            lambda_penalty = st.number_input("Lambda Penalty", min_value=0.0, value=0.8,
                                           help="Penalty for inter-cluster similarity")
            
        # Advanced options in expander
        with st.expander("Advanced Options"):
            sigma_weight = st.number_input("Sigma Weight (small clusters)", min_value=0.1, value=5.0)
            use_trash_cluster = st.checkbox("Use Trash Cluster for Outliers", value=False)
            if use_trash_cluster:
                trash_threshold = st.number_input("Trash Cluster Threshold", min_value=0.0, value=0.0)
            else:
                trash_threshold = 0.0
            
            sequence_weighting = st.selectbox("Sequence Weighting", 
                                            options=[0, 1, 2],
                                            format_func=lambda x: {0: "1/ns", 1: "clustering", 2: "none"}[x],
                                            index=0)
            background_model = st.selectbox("Background Model",
                                          options=[0, 1, 2],
                                          format_func=lambda x: {0: "flat", 1: "pre-calculated", 2: "from data"}[x],
                                          index=1)
    else:
        num_clusters = st.number_input("Number of Clusters", min_value=2, value=4)
        motif_length = 8

    # Define experimental design (only after data is loaded)
    st.header("Experimental Conditions")
    st.write("Define the experimental conditions for each peptide column (if it is either the protein sample or the control sample). Select 'Exclude' to skip a column from analysis.")
    columns_conditions = {}
    for col in data.columns[1:]:  # skip peptide ID column
        condition = st.selectbox(f"Condition for {col}", options=["Control", "Experiment", "Exclude"], index=0)
        if condition != "Exclude":
            columns_conditions[col] = condition

    # Pre-analysis validation (only after data is loaded and conditions are set)
    st.subheader("🔍 Input Validation")
    validation_summary = create_validation_summary(data, columns_conditions, regex)

    # Show validation results
    if validation_summary["overall_status"] == "FAILED":
        st.error("❌ Critical validation errors found:")
        for error in validation_summary["critical_errors"]:
            st.error(f"• {error}")
        st.error("Please fix these errors before running analysis.")
        analysis_allowed = False
    elif validation_summary["overall_status"] == "WARNING":
        st.warning("⚠️ Validation warnings found:")
        for warning in validation_summary["warnings"]:
            st.warning(f"• {warning}")
        
        # Allow user to decide whether to continue with warnings
        continue_with_warnings = st.checkbox("Continue analysis despite warnings", key="continue_warnings")
        analysis_allowed = continue_with_warnings
        
        if continue_with_warnings:
            st.info("Analysis will proceed with warnings noted above.")
    else:
        st.success("✅ Input validation passed")
        analysis_allowed = True

    # Run analysis button (only enabled if validation allows)
    if not analysis_allowed:
        st.button("🚀 Run Analysis", disabled=True, help="Fix validation issues above to enable analysis")
    elif st.button("🚀 Run Analysis"):
        # Create progress tracking
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        status_text.text("Starting analysis pipeline...")
        progress_bar.progress(5)
        
        try:
            # Group by peptide
            status_text.text("Grouping peptides by pattern...")
            progress_bar.progress(10)
            grouped, summary = group_by_peptide(data, columns_conditions, regex)
            
            if grouped.empty:
                st.error("No valid peptides found after grouping. Please check your data and regex pattern.")
                st.info("💡 Try adjusting the regex pattern or check that your peptide column contains valid sequences.")
                st.stop()

            # Display grouped data
            progress_bar.progress(20)
            st.subheader("✅ Step 1: Grouped Data")
            st.dataframe(grouped.head(40))

            st.subheader("Grouping Summary")
            for key, value in summary.items():
                st.write(f"{key}: {value}")

            # Filter by CPM
            status_text.text("Applying CPM filtering...")
            progress_bar.progress(30)
            filtered, fig = filter_by_CPM(grouped, columns_conditions, cpm_threshold, min_count, plot=plot_filtering)
            
            if filtered.empty:
                st.error("No peptides remain after CPM filtering.")
                st.info("💡 Try lowering the CPM threshold or minimum sample count, then run the analysis again.")
                st.info(f"Current settings: CPM threshold = {cpm_threshold}, Min samples = {min_count}")
                st.stop()

            # Display filtering results
            progress_bar.progress(40)
            st.subheader("✅ Step 2: CPM Filtering")
            
            if plot_filtering:
                st.write("Histogram of CPM values before and after filtering")
                if fig is not None:
                    st.pyplot(fig)

            st.dataframe(filtered.head(40))
            st.write(f"📊 Peptides: {grouped.shape[0]} → {filtered.shape[0]} after filtering")

            # Prepare data for DESeq2
            status_text.text("Preparing data for differential expression analysis...")
            progress_bar.progress(50)
            count_data, meta_data = prepare_data(filtered, columns_conditions)

            # Differential expression
            status_text.text("Running DESeq2 differential expression analysis...")
            progress_bar.progress(60)
            res_df, fig = run_deseq2(count_data, meta_data)
            
            st.subheader("✅ Step 3: Differential Expression Results")
            st.write("Results of differential expression analysis using DESeq2")
            st.pyplot(fig)

            # Process significant peptides
            status_text.text("Identifying significant peptides...")
            progress_bar.progress(70)
            res_df, summary = significant_DE_peptides(res_df)

            st.subheader("Significant DE Peptides")
            st.write("Peptides with significant differential expression (FDR < 0.05)")

            for key, value in summary.items():
                st.write(f"{key}: {value}")

            st.subheader("DE Results Table")
            st.dataframe(res_df.head(40))

            # Run VST for heatmap
            status_text.text("Generating visualizations...")
            progress_bar.progress(75)
            vst_data = run_vst(count_data, meta_data)

            upregulated = res_df.loc[res_df["updown"] == "up"].index
            if len(upregulated) > 0:
                vst_up = vst_data.loc[upregulated]
                st.subheader("🔺 Heatmap of Upregulated Peptides")
                st.pyplot(plot_heatmap(vst_up))
            else:
                st.warning("No upregulated peptides found for heatmap")

            # Volcano plot
            st.subheader("📈 Volcano Plot of DE Peptides")
            st.pyplot(volcano_plot(res_df))

            # Sequence processing
            status_text.text("Processing peptide sequences...")
            progress_bar.progress(80)
            const_p = constant_positions(regex)
            st.text("Constant positions in peptides: " + str(const_p))
            res_df["variable_pep"] = res_df.index.map(lambda x: remove_constant(x, const_p))
            res_df["additional_cys"] = res_df["variable_pep"].apply(lambda pep: "C" in pep)
            adj_DE = res_df.loc[~res_df["additional_cys"]].drop(columns=["additional_cys"]).reset_index()

            st.subheader("📑 Upregulated Peptides")
            st.write("List of upregulated peptides after filtering and processing:")
            upregulated_df = adj_DE[adj_DE["updown"] == "up"][["variable_pep", "log2FoldChange", "padj"]]
            if not upregulated_df.empty:
                st.dataframe(upregulated_df)
            else:
                st.warning("No upregulated peptides found after processing")

            # Prepare peptides for clustering
            status_text.text("Preparing clustering analysis...")
            progress_bar.progress(85)
            upregulated_variable_peps = adj_DE[adj_DE["updown"] == "up"][["variable_pep"]]
            
            if upregulated_variable_peps.empty:
                st.warning("No upregulated peptides found for clustering analysis.")
                st.info("💡 This could indicate that no peptides were significantly upregulated, or the fold change/p-value thresholds are too strict.")
                progress_bar.progress(100)
                status_text.text("Analysis completed!")
                st.success("Analysis completed up to differential expression stage ✅")
                st.stop()
            
            # Validate peptide lengths for clustering
            peptide_lengths = upregulated_variable_peps["variable_pep"].str.len()
            min_pep_length = peptide_lengths.min()
            max_pep_length = peptide_lengths.max()
            
            if min_pep_length < motif_length:
                st.error(f"Some peptides are shorter ({min_pep_length}) than the specified motif length ({motif_length}).")
                st.info(f"💡 Please reduce the motif length to {min_pep_length} or less and run the analysis again.")
                progress_bar.progress(100)
                status_text.text("Analysis completed!")
                st.success("Analysis completed up to clustering parameter validation ✅")
                st.stop()
            
            # Validate clustering parameters
            is_valid, errors = validate_clustering_parameters(num_clusters, motif_length, len(upregulated_variable_peps))
            if not is_valid:
                st.error("❌ Clustering parameter validation failed:")
                for error in errors:
                    st.error(f"• {error}")
                st.info("💡 Please adjust the clustering parameters and run the analysis again.")
                progress_bar.progress(100)
                status_text.text("Analysis completed!")
                st.success("Analysis completed up to clustering stage ✅")
                st.stop()
            
            st.info(f"Found {len(upregulated_variable_peps)} upregulated peptides for clustering (length range: {min_pep_length}-{max_pep_length})")
            
            # Clustering analysis section
            st.subheader("✅ Step 4: Clustering Analysis")
            status_text.text("Running clustering analysis...")
            progress_bar.progress(90)
            
            if use_advanced_clustering:
                st.write("Running Advanced Gibbs clustering on upregulated peptides...")
                st.write(f"Configuration: {num_clusters} clusters, {num_seeds} seeds, motif length {motif_length}")
                
                try:
                    # Run advanced clustering
                    clusters_df, pwms, logos = gibbs_cluster_advanced(
                        upregulated_variable_peps["variable_pep"],
                        motif_length=motif_length,
                        num_clusters=num_clusters,
                        num_seeds=num_seeds,
                        iterations=iterations,
                        temperature_start=temperature_start,
                        temperature_steps=temperature_steps,
                        lambda_penalty=lambda_penalty,
                        sigma_weight=sigma_weight,
                        use_trash_cluster=use_trash_cluster,
                        trash_threshold=trash_threshold,
                        sequence_weighting=sequence_weighting,
                        background_model=background_model,
                        seed=42  # For reproducibility
                    )
                    
                    st.success("Advanced clustering completed successfully!")
                    
                    # Show clustering summary if using range
                    if isinstance(num_clusters, tuple):
                        st.subheader("📊 Clustering Quality Analysis")
                        
                        # Create a temporary advanced clusterer to get full results
                        temp_clusterer = GibbsClusterAdvanced(
                            motif_length=motif_length,
                            num_clusters=num_clusters,
                            num_seeds=num_seeds,
                            iterations=iterations,
                            temperature_start=temperature_start,
                            temperature_steps=temperature_steps,
                            lambda_penalty=lambda_penalty,
                            sigma_weight=sigma_weight,
                            use_trash_cluster=use_trash_cluster,
                            trash_threshold=trash_threshold,
                            sequence_weighting=sequence_weighting,
                            background_model=background_model,
                            seed=42
                        )
                        temp_clusterer.fit(upregulated_variable_peps["variable_pep"])
                        summary = temp_clusterer.get_summary()
                        
                        # Create KLD comparison plot
                        cluster_nums = list(summary['cluster_results'].keys())
                        kld_values = [summary['cluster_results'][k]['avg_kld'] for k in cluster_nums]
                        
                        fig_kld, ax_kld = plt.subplots(figsize=(8, 4))
                        ax_kld.bar(cluster_nums, kld_values, color='steelblue', alpha=0.7)
                        ax_kld.set_xlabel('Number of Clusters')
                        ax_kld.set_ylabel('Average KLD')
                        ax_kld.set_title('Clustering Quality vs Number of Clusters')
                        ax_kld.grid(True, alpha=0.3)
                        st.pyplot(fig_kld)
                        
                        # Show detailed results table
                        results_data = []
                        for k in cluster_nums:
                            result = summary['cluster_results'][k]
                            results_data.append({
                                'Clusters': k,
                                'Avg KLD': f"{result['avg_kld']:.4f}",
                                'Non-empty Clusters': result['non_empty_clusters'],
                                'Cluster Sizes': ', '.join(map(str, result['cluster_sizes']))
                            })
                        
                        results_df = pd.DataFrame(results_data)
                        st.dataframe(results_df)
                    
                except Exception as e:
                    st.error(f"Advanced clustering failed: {str(e)}")
                    st.write("Falling back to basic clustering...")
                    
                    clusters_df, pwms, logos = gibbs_cluster(
                        upregulated_variable_peps["variable_pep"],
                        motif_length=motif_length,
                        num_clusters=num_clusters if isinstance(num_clusters, int) else num_clusters[1]
                    )
            else:
                st.write("Running basic Gibbs clustering on upregulated peptides...")
                
                clusters_df, pwms, logos = gibbs_cluster(
                    upregulated_variable_peps["variable_pep"],
                    motif_length=motif_length,
                    num_clusters=num_clusters
                )

            # merge stats
            clusters_df = clusters_df.merge(
                adj_DE.reset_index()[["Clean Peptide", "log2FoldChange", "padj", "-log10(padj)"]],
                left_on="Sequence", right_on="Clean Peptide", how="left"
            )

            # show logos
            st.subheader("🎨 Sequence logos")
            for k in sorted(logos):
                st.pyplot(logos[k])

            st.success("WebLogos generated and saved!")

            # Display cluster boxplots
            st.subheader("📊 Cluster Analysis")
            st.pyplot(box_plot(clusters_df, "Gn", "log2FoldChange", "log2FC per Cluster"))
            st.pyplot(box_plot(clusters_df, "Gn", "-log10(padj)", "-log10(padj) per Cluster"))
            st.pyplot(box_plot(clusters_df, "Gn", "padj", "padj per Cluster"))

            # Final completion
            progress_bar.progress(100)
            status_text.text("Analysis completed!")
            st.success("🎉 Complete Analysis Finished Successfully!")
            
        except Exception as e:
            st.error(f"❌ Analysis failed with error: {str(e)}")
            st.exception(e)
            st.error("Please check your input data and parameters, then try again.")
        
except Exception as e:
    st.error(f"Error loading CSV file: {str(e)}")
    st.error("Please ensure your file is a valid CSV with proper formatting.")
    st.info("💡 Common issues:")
    st.info("• File encoding should be UTF-8")
    st.info("• First row should contain column headers")
    st.info("• At least one column should be named 'peptide' or contain peptide sequences")
    st.stop()
