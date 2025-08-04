import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from Library.lib import *
from Library.Clustering import gibbs_cluster
from Library.GibbsClusterAdvanced import gibbs_cluster_advanced, GibbsClusterAdvanced
from Library.validation import create_validation_summary, validate_clustering_parameters
from Library.enhanced_gibbs_clustering import enhanced_gibbs_cluster_with_validation

# Configure page layout
st.set_page_config(
    page_title="Peptide Analysis Dashboard", 
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for dashboard styling
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: 600;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 2rem;
    }
    .metric-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 1rem;
        border-radius: 10px;
        color: white;
        text-align: center;
        margin: 0.5rem 0;
    }
    .status-success {
        background: linear-gradient(135deg, #4CAF50 0%, #45a049 100%);
        padding: 0.5rem 1rem;
        border-radius: 5px;
        color: white;
        text-align: center;
    }
    .status-warning {
        background: linear-gradient(135deg, #ff9800 0%, #f57c00 100%);
        padding: 0.5rem 1rem;
        border-radius: 5px;
        color: white;
        text-align: center;
    }
    .status-error {
        background: linear-gradient(135deg, #f44336 0%, #d32f2f 100%);
        padding: 0.5rem 1rem;
        border-radius: 5px;
        color: white;
        text-align: center;
    }
    .expandable-card {
        border: 1px solid #ddd;
        border-radius: 10px;
        padding: 1rem;
        margin: 1rem 0;
        background: #fafafa;
    }
    .plot-container {
        border: 2px solid #e0e0e0;
        border-radius: 10px;
        padding: 1rem;
        margin: 1rem 0;
        background: white;
    }
</style>
""", unsafe_allow_html=True)

# Main dashboard header
st.markdown('<h1 class="main-header">🧬 Peptide Analysis Dashboard</h1>', unsafe_allow_html=True)

# Initialize session state for data persistence
if 'analysis_data' not in st.session_state:
    st.session_state.analysis_data = {}
if 'analysis_complete' not in st.session_state:
    st.session_state.analysis_complete = False

# Sidebar for navigation and controls
with st.sidebar:
    st.header("📊 Dashboard Control")
    
    # Analysis progress indicator
    progress_steps = [
        "📁 Data Upload",
        "⚙️ Parameters",
        "🔬 Analysis",
        "📈 Results"
    ]
    
    current_step = 0
    if 'data' in st.session_state:
        current_step = 1
    if st.session_state.get('parameters_set', False):
        current_step = 2
    if st.session_state.analysis_complete:
        current_step = 3
    
    st.subheader("Progress")
    for i, step in enumerate(progress_steps):
        if i <= current_step:
            st.markdown(f"✅ {step}")
        else:
            st.markdown(f"⏳ {step}")
    
    st.divider()
    
    # Quick actions
    st.subheader("Quick Actions")
    if st.button("🔄 Reset Analysis", type="secondary"):
        for key in list(st.session_state.keys()):
            if key != 'analysis_data':
                del st.session_state[key]
        st.rerun()
    
    if st.button("💾 Download Results", disabled=not st.session_state.analysis_complete):
        st.info("Results download would be implemented here")

# Main dashboard layout
col1, col2 = st.columns([2, 1])

with col1:
    # Data upload section
    st.subheader("📁 Data Upload & Preview")
    
    uploaded_file = st.file_uploader(
        "Upload your peptide CSV file", 
        type=["csv"],
        help="CSV file with peptide sequences and count data"
    )
    
    if uploaded_file is not None:
        try:
            data = pd.read_csv(uploaded_file)
            st.session_state['data'] = data
            
            # Data preview in expandable container
            with st.expander("📋 Data Preview", expanded=True):
                col_prev1, col_prev2 = st.columns(2)
                
                with col_prev1:
                    st.markdown(f"""
                    <div class="metric-card">
                        <h3>{data.shape[0]:,}</h3>
                        <p>Total Sequences</p>
                    </div>
                    """, unsafe_allow_html=True)
                
                with col_prev2:
                    st.markdown(f"""
                    <div class="metric-card">
                        <h3>{data.shape[1]-1}</h3>
                        <p>Sample Columns</p>
                    </div>
                    """, unsafe_allow_html=True)
                
                st.dataframe(data.head(10), use_container_width=True)
                
        except Exception as e:
            st.error(f"Error loading file: {e}")

with col2:
    # Status panel
    st.subheader("📊 Analysis Status")
    
    if uploaded_file is None:
        st.markdown('<div class="status-warning">⏳ Waiting for data upload</div>', unsafe_allow_html=True)
    elif 'data' in st.session_state:
        st.markdown('<div class="status-success">✅ Data loaded successfully</div>', unsafe_allow_html=True)
        
        # Quick stats
        data = st.session_state['data']
        st.metric("Sequences", f"{data.shape[0]:,}")
        st.metric("Samples", f"{data.shape[1]-1}")
        st.metric("Data Size", f"{data.memory_usage(deep=True).sum() / 1024**2:.1f} MB")

# Parameters section (only show if data is loaded)
if 'data' in st.session_state:
    data = st.session_state['data']
    
    st.divider()
    st.subheader("⚙️ Analysis Parameters")
    
    # Parameters in collapsible sections
    param_col1, param_col2 = st.columns(2)
    
    with param_col1:
        with st.expander("🔍 Pattern & Filtering", expanded=True):
            regex = st.text_input(
                "Peptide Pattern", 
                value=r'A.C.{7}C',
                help="Regular expression for peptide library pattern"
            )
            
            cpm_threshold = st.number_input(
                "CPM Threshold", 
                min_value=0, 
                value=5,
                help="Minimum counts per million for filtering"
            )
            
            min_count = st.number_input(
                "Min Samples", 
                min_value=1, 
                value=1,
                help="Minimum samples exceeding CPM threshold"
            )
            
            plot_filtering = st.checkbox("Show filtering plots", value=True)
    
    with param_col2:
        with st.expander("🧬 Clustering Options", expanded=True):
            use_advanced_clustering = st.checkbox("Advanced Clustering", value=True)
            use_consensus_validation = st.checkbox("Consensus Validation", value=True)
            
            if use_consensus_validation:
                auto_k_selection = st.checkbox("Auto K Selection", value=True)
                if auto_k_selection:
                    k_range_option = st.selectbox(
                        "K Search Range",
                        ["Narrow (2-5)", "Medium (2-7)", "Wide (2-10)"],
                        index=1
                    )
                    range_mapping = {
                        "Narrow (2-5)": (2, 5),
                        "Medium (2-7)": (2, 7), 
                        "Wide (2-10)": (2, 10)
                    }
                    k_range = range_mapping[k_range_option]
                else:
                    num_clusters = st.number_input("Number of Clusters", min_value=2, value=4)
                    k_range = num_clusters
            else:
                num_clusters = st.number_input("Number of Clusters", min_value=2, value=4)
                k_range = num_clusters
            
            motif_length = st.number_input("Motif Length", min_value=6, value=8)
            
            # Advanced clustering hyperparameters
            show_advanced_params = st.checkbox("Show Advanced Parameters", value=False)
            
            if show_advanced_params:
                st.markdown("**🔧 Advanced Clustering Parameters**")
                
                # Sampling parameters
                col_adv1, col_adv2 = st.columns(2)
                with col_adv1:
                    st.markdown("*Sampling:*")
                    num_seeds = st.number_input("Number of Seeds", min_value=1, max_value=10, value=3, help="Independent sampling runs for convergence assessment")
                    iterations = st.number_input("Iterations per Temperature", min_value=5, max_value=50, value=10, help="Gibbs sampling iterations at each temperature")
                    n_iter_basic = st.number_input("Basic Clustering Iterations", min_value=100, max_value=5000, value=1000, help="Total iterations for basic clustering")
                
                with col_adv2:
                    st.markdown("*Temperature Annealing:*")
                    temperature_start = st.number_input("Start Temperature", min_value=0.5, max_value=3.0, value=1.5, step=0.1, help="Initial sampling temperature for exploration")
                    temperature_steps = st.number_input("Temperature Steps", min_value=10, max_value=50, value=20, help="Number of temperature reduction steps")
                
                # Regularization parameters
                col_reg1, col_reg2 = st.columns(2)
                with col_reg1:
                    st.markdown("*Regularization:*")
                    lambda_penalty = st.number_input("Lambda Penalty", min_value=0.0, max_value=2.0, value=0.8, step=0.1, help="Penalty for similarity between clusters")
                    sigma_weight = st.number_input("Sigma Weight", min_value=1.0, max_value=10.0, value=5.0, step=0.5, help="Weight penalty for small clusters")
                
                with col_reg2:
                    st.markdown("*Model Parameters:*")
                    sequence_weighting = st.selectbox("Sequence Weighting", 
                                                    options=[0, 1, 2], 
                                                    index=0,
                                                    format_func=lambda x: {0: "1/ns (Recommended)", 1: "Clustering", 2: "None"}[x],
                                                    help="Weighting scheme for sequences")
                    background_model = st.selectbox("Background Model",
                                                  options=[0, 1, 2],
                                                  index=1,
                                                  format_func=lambda x: {0: "Flat", 1: "BLOSUM (Recommended)", 2: "From Data"}[x],
                                                  help="Background amino acid frequency model")
                
                # Validation parameters
                if use_consensus_validation:
                    st.markdown("*Validation Parameters:*")
                    col_val1, col_val2 = st.columns(2)
                    with col_val1:
                        n_consensus_iterations = st.number_input("Consensus Iterations", min_value=10, max_value=200, value=50, help="Subsampling iterations for stability")
                        enable_bootstrap = st.checkbox("Enable Bootstrap", value=True, help="Additional bootstrap validation")
                    
                    with col_val2:
                        if enable_bootstrap:
                            n_bootstrap_iterations = st.number_input("Bootstrap Iterations", min_value=10, max_value=100, value=30, help="Bootstrap resampling iterations")
                        else:
                            n_bootstrap_iterations = 0
                        
                        sample_fraction = st.number_input("Sample Fraction", min_value=0.5, max_value=0.9, value=0.8, step=0.05, help="Fraction of sequences for subsampling")
                
                # Outlier handling
                st.markdown("*Outlier Handling:*")
                col_out1, col_out2 = st.columns(2)
                with col_out1:
                    use_trash_cluster = st.checkbox("Use Trash Cluster", value=False, help="Create cluster for outlier sequences")
                with col_out2:
                    if use_trash_cluster:
                        trash_threshold = st.number_input("Trash Threshold", min_value=0.0, max_value=0.5, value=0.1, step=0.05, help="Probability threshold for trash assignment")
                    else:
                        trash_threshold = 0.0
            else:
                # Default values when advanced parameters are hidden
                num_seeds = 3
                iterations = 10
                n_iter_basic = 1000
                temperature_start = 1.5
                temperature_steps = 20
                lambda_penalty = 0.8
                sigma_weight = 5.0
                sequence_weighting = 0
                background_model = 1
                n_consensus_iterations = 50
                enable_bootstrap = True
                n_bootstrap_iterations = 30
                sample_fraction = 0.8
                use_trash_cluster = False
                trash_threshold = 0.0
    
    # Experimental conditions
    st.subheader("🔬 Experimental Design")
    
    conditions_col1, conditions_col2 = st.columns(2)
    
    with conditions_col1:
        st.write("**Condition Assignment**")
        columns_conditions = {}
        for col in data.columns[1:]:
            condition = st.selectbox(
                f"{col}", 
                options=["Control", "Experiment", "Exclude"], 
                index=0,
                key=f"condition_{col}"
            )
            if condition != "Exclude":
                columns_conditions[col] = condition
    
    with conditions_col2:
        # Validation preview
        if columns_conditions:
            validation_summary = create_validation_summary(data, columns_conditions, regex)
            
            st.write("**Validation Status**")
            if validation_summary["overall_status"] == "PASSED":
                st.markdown('<div class="status-success">✅ All validations passed</div>', unsafe_allow_html=True)
            elif validation_summary["overall_status"] == "WARNING":
                st.markdown('<div class="status-warning">⚠️ Warnings detected</div>', unsafe_allow_html=True)
                for warning in validation_summary["warnings"][:3]:  # Show first 3 warnings
                    st.warning(f"• {warning}")
            else:
                st.markdown('<div class="status-error">❌ Critical errors found</div>', unsafe_allow_html=True)
                for error in validation_summary["critical_errors"][:3]:  # Show first 3 errors
                    st.error(f"• {error}")
    
    # Analysis button
    st.divider()
    
    analysis_allowed = validation_summary["overall_status"] != "FAILED"
    if validation_summary["overall_status"] == "WARNING":
        continue_with_warnings = st.checkbox("Continue with warnings", key="continue_warnings")
        analysis_allowed = continue_with_warnings
    
    col_btn1, col_btn2, col_btn3 = st.columns([1, 2, 1])
    with col_btn2:
        if st.button(
            "🚀 Run Analysis", 
            disabled=not analysis_allowed,
            use_container_width=True,
            type="primary"
        ):
            st.session_state['parameters_set'] = True
            
            # Store parameters
            st.session_state['analysis_params'] = {
                'regex': regex,
                'cpm_threshold': cpm_threshold,
                'min_count': min_count,
                'columns_conditions': columns_conditions,
                'use_advanced_clustering': use_advanced_clustering,
                'use_consensus_validation': use_consensus_validation,
                'k_range': k_range,
                'motif_length': motif_length,
                'plot_filtering': plot_filtering,
                # Advanced clustering parameters
                'num_seeds': num_seeds,
                'iterations': iterations,
                'n_iter_basic': n_iter_basic,
                'temperature_start': temperature_start,
                'temperature_steps': temperature_steps,
                'lambda_penalty': lambda_penalty,
                'sigma_weight': sigma_weight,
                'sequence_weighting': sequence_weighting,
                'background_model': background_model,
                'n_consensus_iterations': n_consensus_iterations,
                'enable_bootstrap': enable_bootstrap,
                'n_bootstrap_iterations': n_bootstrap_iterations,
                'sample_fraction': sample_fraction,
                'use_trash_cluster': use_trash_cluster,
                'trash_threshold': trash_threshold
            }
            
            # Run analysis
            with st.spinner("Running analysis pipeline..."):
                try:
                    # Create progress container
                    progress_container = st.container()
                    with progress_container:
                        progress_bar = st.progress(0)
                        status_text = st.empty()
                    
                    # Step 1: Group peptides
                    status_text.text("Grouping peptides by pattern...")
                    progress_bar.progress(10)
                    grouped, summary = group_by_peptide(data, columns_conditions, regex)
                    
                    if grouped.empty:
                        st.error("No valid peptides found after grouping.")
                        st.stop()
                    
                    # Step 2: CPM filtering
                    status_text.text("Applying CPM filtering...")
                    progress_bar.progress(30)
                    filtered, fig_filter = filter_by_CPM(grouped, columns_conditions, cpm_threshold, min_count, plot=plot_filtering)
                    
                    if filtered.empty:
                        st.error("No peptides remain after CPM filtering.")
                        st.stop()
                    
                    # Step 3: DESeq2
                    status_text.text("Running differential expression analysis...")
                    progress_bar.progress(50)
                    count_data, meta_data = prepare_data(filtered, columns_conditions)
                    res_df, fig_ma = run_deseq2(count_data, meta_data)
                    
                    # Step 4: Significant peptides
                    status_text.text("Identifying significant peptides...")
                    progress_bar.progress(70)
                    res_df, de_summary = significant_DE_peptides(res_df)
                    
                    # Step 5: VST and visualizations
                    status_text.text("Generating visualizations...")
                    progress_bar.progress(80)
                    vst_data = run_vst(count_data, meta_data)
                    fig_volcano = volcano_plot(res_df)
                    
                    # Step 6: Clustering
                    status_text.text("Running clustering analysis...")
                    progress_bar.progress(90)
                    
                    # Process sequences for clustering
                    const_p = constant_positions(regex)
                    res_df["variable_pep"] = res_df.index.map(lambda x: remove_constant(x, const_p))
                    res_df["additional_cys"] = res_df["variable_pep"].apply(lambda pep: "C" in pep)
                    adj_DE = res_df.loc[~res_df["additional_cys"]].drop(columns=["additional_cys"]).reset_index()
                    upregulated_variable_peps = adj_DE[adj_DE["updown"] == "up"][["variable_pep"]]
                    
                    if not upregulated_variable_peps.empty:
                        if use_consensus_validation:
                            clusters_df, pwms, logos, validation_results = enhanced_gibbs_cluster_with_validation(
                                upregulated_variable_peps["variable_pep"].tolist(),
                                motif_length=motif_length,
                                k_range=k_range,
                                use_advanced=use_advanced_clustering,
                                validation_enabled=True,
                                n_consensus_iterations=n_consensus_iterations,
                                n_bootstrap_iterations=n_bootstrap_iterations,
                                enable_bootstrap=enable_bootstrap,
                                sample_fraction=sample_fraction,
                                random_seed=42,
                                # Advanced clustering parameters
                                num_seeds=num_seeds,
                                iterations=iterations,
                                temperature_start=temperature_start,
                                temperature_steps=temperature_steps,
                                lambda_penalty=lambda_penalty,
                                sigma_weight=sigma_weight,
                                sequence_weighting=sequence_weighting,
                                background_model=background_model,
                                use_trash_cluster=use_trash_cluster,
                                trash_threshold=trash_threshold
                            )
                        else:
                            if use_advanced_clustering:
                                clusters_df, pwms, logos = gibbs_cluster_advanced(
                                    upregulated_variable_peps["variable_pep"],
                                    motif_length=motif_length,
                                    num_clusters=k_range if isinstance(k_range, int) else 4,
                                    seed=42,
                                    # Advanced parameters
                                    num_seeds=num_seeds,
                                    iterations=iterations,
                                    temperature_start=temperature_start,
                                    temperature_steps=temperature_steps,
                                    lambda_penalty=lambda_penalty,
                                    sigma_weight=sigma_weight,
                                    sequence_weighting=sequence_weighting,
                                    background_model=background_model,
                                    use_trash_cluster=use_trash_cluster,
                                    trash_threshold=trash_threshold
                                )
                                validation_results = None
                            else:
                                clusters_df, pwms, logos = gibbs_cluster(
                                    upregulated_variable_peps["variable_pep"],
                                    motif_length=motif_length,
                                    num_clusters=k_range if isinstance(k_range, int) else 4,
                                    n_iter=n_iter_basic,
                                    seed=42
                                )
                                validation_results = None
                        
                        # Merge clustering results with DE data
                        clusters_df = clusters_df.merge(
                            adj_DE.reset_index()[["Clean Peptide", "log2FoldChange", "padj", "-log10(padj)"]],
                            left_on="Sequence", right_on="Clean Peptide", how="left"
                        )
                    else:
                        clusters_df = pd.DataFrame()
                        logos = {}
                        validation_results = None
                    
                    # Store results
                    st.session_state['analysis_results'] = {
                        'grouped': grouped,
                        'filtered': filtered,
                        'res_df': res_df,
                        'de_summary': de_summary,
                        'clusters_df': clusters_df,
                        'logos': logos,
                        'validation_results': validation_results,
                        'vst_data': vst_data,
                        'figures': {
                            'filter': fig_filter,
                            'ma_plot': fig_ma,
                            'volcano': fig_volcano
                        }
                    }
                    
                    progress_bar.progress(100)
                    status_text.text("Analysis completed!")
                    st.session_state.analysis_complete = True
                    
                    st.success("🎉 Analysis completed successfully!")
                    st.rerun()
                    
                except Exception as e:
                    st.error(f"❌ Analysis failed: {str(e)}")
                    st.exception(e)

# Results section (only show if analysis is complete)
if st.session_state.analysis_complete and 'analysis_results' in st.session_state:
    results = st.session_state['analysis_results']
    
    st.divider()
    st.subheader("📈 Analysis Results Dashboard")
    
    # Summary metrics
    col_m1, col_m2, col_m3, col_m4 = st.columns(4)
    
    with col_m1:
        st.markdown(f"""
        <div class="metric-card">
            <h3>{len(results['grouped']):,}</h3>
            <p>Grouped Peptides</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col_m2:
        st.markdown(f"""
        <div class="metric-card">
            <h3>{len(results['filtered']):,}</h3>
            <p>After CPM Filter</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col_m3:
        up_count = results['de_summary']['Upregulated Peptides']
        st.markdown(f"""
        <div class="metric-card">
            <h3>{up_count:,}</h3>
            <p>Upregulated</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col_m4:
        cluster_count = len(results['clusters_df']) if not results['clusters_df'].empty else 0
        st.markdown(f"""
        <div class="metric-card">
            <h3>{cluster_count:,}</h3>
            <p>Clustered</p>
        </div>
        """, unsafe_allow_html=True)
    
    # Results tabs
    tab1, tab2, tab3, tab4 = st.tabs(["📊 Filtering", "🧮 Differential Expression", "🔬 Clustering", "📋 Summary"])
    
    with tab1:
        
        st.subheader("Filtering Summary")
        st.write(f"**Success Rate**: {results['grouped'].shape[0]} / {st.session_state['data'].shape[0]} sequences matched pattern")
        st.write(f"**After CPM Filter**: {results['filtered'].shape[0]} peptides retained")
        
        if not results['filtered'].empty:
            st.dataframe(results['filtered'].head(10))

    with tab2:
        col_de1, col_de2 = st.columns(2)
        
        with col_de1:
            st.subheader("MA Plot")
            if results['figures']['ma_plot']:
                st.pyplot(results['figures']['ma_plot'])
        
        with col_de2:
            st.subheader("Volcano Plot")
            if results['figures']['volcano']:
                st.pyplot(results['figures']['volcano'])
        
        # DE summary
        st.subheader("Differential Expression Summary")
        summary_data = pd.DataFrame([results['de_summary']]).T
        summary_data.columns = ['Count']
        st.dataframe(summary_data, use_container_width=True)
    
    with tab3:
        if not results['clusters_df'].empty:
            # Validation results if available
            if results['validation_results']:
                val_res = results['validation_results']
                st.subheader("🔬 Consensus Validation Results")
                
                col_v1, col_v2, col_v3 = st.columns(3)
                with col_v1:
                    st.metric("Optimal K", val_res['optimal_k'])
                with col_v2:
                    st.metric("Stability Score", f"{val_res['stability_score']:.3f}")
                with col_v3:
                    st.metric("Silhouette Score", f"{val_res['silhouette_score']:.3f}")
                
                # Quality interpretation
                silhouette_score = val_res['silhouette_score']
                if silhouette_score > 0.7:
                    st.success("🎯 Excellent clustering quality")
                elif silhouette_score > 0.5:
                    st.info("✅ Good clustering quality")
                elif silhouette_score > 0.25:
                    st.warning("⚠️ Moderate clustering quality")
                else:
                    st.error("❌ Poor clustering quality")
            
            # Sequence logos
            st.subheader("🎨 Sequence Logos")
            if results['logos']:
                logo_cols = st.columns(min(len(results['logos']), 3))
                for i, (k, fig) in enumerate(sorted(results['logos'].items())):
                    with logo_cols[i % 3]:
                        st.write(f"**Cluster {k}**")
                        st.pyplot(fig)
            
            # Cluster data table
            st.subheader("📊 Cluster Assignments")
            st.dataframe(results['clusters_df'].head(20), use_container_width=True)
        else:
            st.warning("No upregulated peptides found for clustering analysis.")
    
    with tab4:
        st.subheader("📋 Complete Analysis Summary")
        
        # Analysis parameters used
        st.write("**Parameters Used:**")
        params = st.session_state['analysis_params']
        param_df = pd.DataFrame([
            ["Peptide Pattern", params['regex']],
            ["CPM Threshold", params['cpm_threshold']],
            ["Min Samples", params['min_count']],
            ["Clustering Method", "Advanced" if params['use_advanced_clustering'] else "Basic"],
            ["Validation", "Enabled" if params['use_consensus_validation'] else "Disabled"],
            ["Motif Length", params['motif_length']]
        ], columns=['Parameter', 'Value'])
        st.dataframe(param_df, use_container_width=True)
        
        # Download options
        st.subheader("💾 Download Results")
        col_dl1, col_dl2, col_dl3 = st.columns(3)
        
        with col_dl1:
            if st.button("📊 Download DE Results"):
                csv = results['res_df'].to_csv()
                st.download_button(
                    label="Download CSV",
                    data=csv,
                    file_name="differential_expression_results.csv",
                    mime="text/csv"
                )
        
        with col_dl2:
            if st.button("🔬 Download Clusters") and not results['clusters_df'].empty:
                csv = results['clusters_df'].to_csv()
                st.download_button(
                    label="Download CSV",
                    data=csv,
                    file_name="clustering_results.csv",
                    mime="text/csv"
                )
        
        with col_dl3:
            if st.button("📈 Download All Plots"):
                st.info("Plot download functionality would be implemented here")

# Footer
st.divider()
st.markdown("""
<div style="text-align: center; color: #666; padding: 2rem;">
    🧬 Enhanced Peptide Analysis Dashboard | Built with Streamlit
</div>
""", unsafe_allow_html=True)