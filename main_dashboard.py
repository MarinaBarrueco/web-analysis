#!/usr/bin/env python3
"""
Peptide Analysis Dashboard

Interactive Streamlit dashboard for comprehensive peptide library analysis including:
- Pattern-based filtering and CPM normalization
- DESeq2-based differential expression analysis  
- Advanced Gibbs clustering with consensus validation
- Publication-ready visualizations and quality metrics

Usage:
    streamlit run main_dashboard.py

Author: Peptide Analysis Pipeline
Version: 2.0
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests

# Import analysis modules
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

# ===================================================================
# CONSTANTS AND CONFIGURATION
# ===================================================================

# Default parameter values
DEFAULT_PARAMS = {
    'regex': r'A.C.{7}C',
    'cpm_threshold': 5,
    'min_count': 1,
    'motif_length': 8,
    'num_clusters': 4,
    'num_seeds': 3,
    'iterations': 10,
    'n_iter_basic': 1000,
    'temperature_start': 1.5,
    'temperature_steps': 20,
    'lambda_penalty': 0.8,
    'sigma_weight': 5.0,
    'sequence_weighting': 0,
    'background_model': 1,
    'n_consensus_iterations': 50,
    'n_bootstrap_iterations': 30,
    'sample_fraction': 0.8,
    'use_trash_cluster': False,
    'trash_threshold': 0.0
}

# Analysis progress steps
ANALYSIS_STEPS = [
    "📁 Data Upload",
    "⚙️ Parameters", 
    "🔬 Analysis",
    "📈 Results"
]

# K range mappings for auto-selection
K_RANGE_OPTIONS = {
    "Narrow (2-5)": (2, 5),
    "Medium (2-7)": (2, 7),
    "Wide (2-10)": (2, 10)
}

# ===================================================================
# HELPER FUNCTIONS
# ===================================================================

def get_analysis_progress_step():
    """Determine current analysis progress step"""
    if 'data' not in st.session_state:
        return 0
    if not st.session_state.get('parameters_set', False):
        return 1
    if not st.session_state.analysis_complete:
        return 2
    return 3

def validate_data_format(data: pd.DataFrame) -> tuple[bool, list[str]]:
    """Validate uploaded data format and return status with messages"""
    errors = []
    
    if data.empty:
        errors.append("Data file is empty")
    
    if 'peptide' not in data.columns:
        errors.append("Missing 'peptide' column")
    
    if len(data.columns) < 2:
        errors.append("Need at least one sample column besides 'peptide'")
    
    # Check for numeric data in sample columns
    sample_cols = [col for col in data.columns if col != 'peptide']
    for col in sample_cols:
        if not pd.api.types.is_numeric_dtype(data[col]):
            try:
                data[col] = pd.to_numeric(data[col], errors='coerce')
            except:
                errors.append(f"Column '{col}' contains non-numeric data")
    
    return len(errors) == 0, errors

def create_metric_card(title: str, value: str, help_text: str = ""):
    """Create a styled metric card"""
    help_tooltip = f"help='{help_text}'" if help_text else ""
    return f"""
    <div class="metric-card" {help_tooltip}>
        <h3>{value}</h3>
        <p>{title}</p>
    </div>
    """

# ===================================================================
# MAIN DASHBOARD
# ===================================================================

# Main dashboard header
st.markdown('<h1 class="main-header">🧬 Peptide Analysis Dashboard</h1>', unsafe_allow_html=True)

# Help and Instructions Section
with st.expander("📖 How to Use This Dashboard", expanded=False):
    st.markdown("""
    ## 🚀 Quick Start Guide
    
    ### 📁 **Step 1: Prepare Your Data**
    Upload a CSV file with the following structure:
    - **First column**: `peptide` - Contains your peptide sequences
    - **Remaining columns**: Sample data with count values for each condition
    
    **Example CSV format:**
    ```
    peptide,Sample1_Control,Sample2_Control,Sample1_Treatment,Sample2_Treatment
    AYCPFRSWPGCGG,1250,890,2340,1980
    AFCSLWRSGDCGG,890,1100,0,45
    ```
    
    ### ⚙️ **Step 2: Configure Analysis Parameters**
    
    **🔍 Pattern & Filtering:**
    - **Peptide Pattern**: Regular expression to match your library design (e.g., `A.C.{7}C` for A-X-C-XXXXXXX-C)
    - **CPM Threshold**: Minimum counts per million for filtering (default: 5)
    - **Min Samples**: Minimum samples that must exceed CPM threshold (default: 1)
    
    **🧬 Clustering Modes:**
    - **Simple**: Basic clustering with fixed number of clusters (recommended for beginners)
    - **Advanced**: Enhanced clustering features with fixed K
    - **Auto-Selection**: Automatically finds optimal number of clusters using validation metrics
    
    ### 🔬 **Step 3: Set Experimental Conditions**
    Assign each sample column to either:
    - **Control**: Reference condition
    - **Experiment**: Treatment condition
    - **Exclude**: Skip this sample
    
    ### 📊 **Step 4: Run Analysis**
    The pipeline will automatically:
    1. Filter peptides by pattern and CPM thresholds
    2. Perform differential expression analysis (DESeq2)
    3. Identify significant peptides
    4. Cluster upregulated sequences
    5. Generate sequence logos and quality metrics
    
    ## 🎯 **Understanding Results**
    
    **📊 Filtering Tab**: Shows data quality and filtering effectiveness
    
    **🧮 Differential Expression Tab**: 
    - MA Plot: Expression vs. fold change
    - Volcano Plot: Significance vs. fold change
    
    **🔬 Clustering Tab**:
    - Quality metrics assess clustering reliability
    - Sequence logos show motif patterns
    - Box plots display differential expression by cluster
    
    ## 💡 **Tips for Best Results**
    - Ensure adequate sequencing depth (>1000 reads per sample)
    - Use biological replicates for robust statistical analysis  
    - Choose appropriate CPM thresholds based on your data distribution
    - Validate clustering results with known motifs if available
    """)

# Initialize session state for data persistence
if 'analysis_data' not in st.session_state:
    st.session_state.analysis_data = {}
if 'analysis_complete' not in st.session_state:
    st.session_state.analysis_complete = False

# Sidebar for navigation and controls
with st.sidebar:
    st.header("📊 Dashboard Control")
    
    # Analysis progress indicator
    current_step = get_analysis_progress_step()
    
    st.subheader("Progress")
    for i, step in enumerate(ANALYSIS_STEPS):
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
            
            # Validate data format
            is_valid, errors = validate_data_format(data)
            
            if not is_valid:
                st.error("❌ Data validation failed:")
                for error in errors:
                    st.error(f"• {error}")
                st.stop()
            
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
                value=DEFAULT_PARAMS['regex'],
                help="Regular expression for peptide library pattern"
            )
            
            cpm_threshold = st.number_input(
                "CPM Threshold", 
                min_value=0, 
                value=DEFAULT_PARAMS['cpm_threshold'],
                help="Minimum counts per million for filtering"
            )
            
            min_count = st.number_input(
                "Min Samples", 
                min_value=1, 
                value=DEFAULT_PARAMS['min_count'],
                help="Minimum samples exceeding CPM threshold"
            )
            
            plot_filtering = st.checkbox("Show filtering plots", value=True)
    
    with param_col2:
        with st.expander("🧬 Clustering Options", expanded=True):
            # Simplified clustering mode selection
            clustering_mode = st.selectbox(
                "Clustering Mode",
                ["Simple", "Advanced", "Auto-Selection"],
                index=0,  # Default to Simple
                help="Simple: Basic clustering with fixed K | Advanced: Enhanced features | Auto-Selection: Finds optimal K"
            )
            
            # Set parameters based on mode
            if clustering_mode == "Simple":
                use_advanced_clustering = False
                use_consensus_validation = False
                num_clusters = st.number_input("Number of Clusters", min_value=2, max_value=8, value=DEFAULT_PARAMS['num_clusters'])
                k_range = num_clusters
                
            elif clustering_mode == "Advanced":
                use_advanced_clustering = True
                use_consensus_validation = False
                num_clusters = st.number_input("Number of Clusters", min_value=2, max_value=8, value=DEFAULT_PARAMS['num_clusters'])
                k_range = num_clusters
                
            else:  # Auto-Selection
                use_advanced_clustering = True
                use_consensus_validation = True
                k_range_option = st.selectbox(
                    "K Search Range",
                    list(K_RANGE_OPTIONS.keys()),
                    index=1
                )
                k_range = K_RANGE_OPTIONS[k_range_option]
                num_clusters = DEFAULT_PARAMS['num_clusters']  # fallback value
            
            motif_length = st.number_input("Motif Length", min_value=6, max_value=12, value=DEFAULT_PARAMS['motif_length'])
            
            # Advanced clustering hyperparameters (only show for Advanced/Auto modes)
            show_advanced_params = False
            if clustering_mode in ["Advanced", "Auto-Selection"]:
                show_advanced_params = st.checkbox("Show Advanced Parameters", value=False)
            
            if show_advanced_params:
                st.markdown("**🔧 Advanced Clustering Parameters**")
                
                # Sampling parameters
                col_adv1, col_adv2 = st.columns(2)
                with col_adv1:
                    st.markdown("*Sampling:*")
                    num_seeds = st.number_input("Number of Seeds", min_value=1, max_value=10, value=DEFAULT_PARAMS['num_seeds'], help="Independent sampling runs for convergence assessment")
                    iterations = st.number_input("Iterations per Temperature", min_value=5, max_value=50, value=DEFAULT_PARAMS['iterations'], help="Gibbs sampling iterations at each temperature")
                    n_iter_basic = st.number_input("Basic Clustering Iterations", min_value=100, max_value=5000, value=DEFAULT_PARAMS['n_iter_basic'], help="Total iterations for basic clustering")
                
                with col_adv2:
                    st.markdown("*Temperature Annealing:*")
                    temperature_start = st.number_input("Start Temperature", min_value=0.5, max_value=3.0, value=DEFAULT_PARAMS['temperature_start'], step=0.1, help="Initial sampling temperature for exploration")
                    temperature_steps = st.number_input("Temperature Steps", min_value=10, max_value=50, value=DEFAULT_PARAMS['temperature_steps'], help="Number of temperature reduction steps")
                
                # Regularization parameters
                col_reg1, col_reg2 = st.columns(2)
                with col_reg1:
                    st.markdown("*Regularization:*")
                    lambda_penalty = st.number_input("Lambda Penalty", min_value=0.0, max_value=2.0, value=DEFAULT_PARAMS['lambda_penalty'], step=0.1, help="Penalty for similarity between clusters")
                    sigma_weight = st.number_input("Sigma Weight", min_value=1.0, max_value=10.0, value=DEFAULT_PARAMS['sigma_weight'], step=0.5, help="Weight penalty for small clusters")
                
                with col_reg2:
                    st.markdown("*Model Parameters:*")
                    sequence_weighting = st.selectbox("Sequence Weighting", 
                                                    options=[0, 1, 2], 
                                                    index=DEFAULT_PARAMS['sequence_weighting'],
                                                    format_func=lambda x: {0: "1/ns (Recommended)", 1: "Clustering", 2: "None"}[x],
                                                    help="Weighting scheme for sequences")
                    background_model = st.selectbox("Background Model",
                                                  options=[0, 1, 2],
                                                  index=DEFAULT_PARAMS['background_model'],
                                                  format_func=lambda x: {0: "Flat", 1: "BLOSUM (Recommended)", 2: "From Data"}[x],
                                                  help="Background amino acid frequency model")
                
                # Validation parameters
                if use_consensus_validation:
                    st.markdown("*Validation Parameters:*")
                    col_val1, col_val2 = st.columns(2)
                    with col_val1:
                        n_consensus_iterations = st.number_input("Consensus Iterations", min_value=10, max_value=200, value=DEFAULT_PARAMS['n_consensus_iterations'], help="Subsampling iterations for stability")
                        enable_bootstrap = st.checkbox("Enable Bootstrap", value=True, help="Additional bootstrap validation")
                    
                    with col_val2:
                        if enable_bootstrap:
                            n_bootstrap_iterations = st.number_input("Bootstrap Iterations", min_value=10, max_value=100, value=DEFAULT_PARAMS['n_bootstrap_iterations'], help="Bootstrap resampling iterations")
                        else:
                            n_bootstrap_iterations = 0
                        
                        sample_fraction = st.number_input("Sample Fraction", min_value=0.5, max_value=0.9, value=DEFAULT_PARAMS['sample_fraction'], step=0.05, help="Fraction of sequences for subsampling")
                
                # Outlier handling
                st.markdown("*Outlier Handling:*")
                col_out1, col_out2 = st.columns(2)
                with col_out1:
                    use_trash_cluster = st.checkbox("Use Trash Cluster", value=DEFAULT_PARAMS['use_trash_cluster'], help="Create cluster for outlier sequences")
                with col_out2:
                    if use_trash_cluster:
                        trash_threshold = st.number_input("Trash Threshold", min_value=0.0, max_value=0.5, value=DEFAULT_PARAMS['trash_threshold'], step=0.05, help="Probability threshold for trash assignment")
                    else:
                        trash_threshold = DEFAULT_PARAMS['trash_threshold']
            else:
                # Default values when advanced parameters are hidden
                num_seeds = DEFAULT_PARAMS['num_seeds']
                iterations = DEFAULT_PARAMS['iterations']
                n_iter_basic = DEFAULT_PARAMS['n_iter_basic']
                temperature_start = DEFAULT_PARAMS['temperature_start']
                temperature_steps = DEFAULT_PARAMS['temperature_steps']
                lambda_penalty = DEFAULT_PARAMS['lambda_penalty']
                sigma_weight = DEFAULT_PARAMS['sigma_weight']
                sequence_weighting = DEFAULT_PARAMS['sequence_weighting']
                background_model = DEFAULT_PARAMS['background_model']
                n_consensus_iterations = DEFAULT_PARAMS['n_consensus_iterations']
                enable_bootstrap = True
                n_bootstrap_iterations = DEFAULT_PARAMS['n_bootstrap_iterations']
                sample_fraction = DEFAULT_PARAMS['sample_fraction']
                use_trash_cluster = DEFAULT_PARAMS['use_trash_cluster']
                trash_threshold = DEFAULT_PARAMS['trash_threshold']
    
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
                    # filtered, fig_filter = filter_by_CPM(grouped, columns_conditions, cpm_threshold, min_count, plot=plot_filtering)
                    filtered, fig_filter = filter_by_CPM_v2_style(grouped, columns_conditions, cpm_threshold, min_count, plot=plot_filtering)

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
                        # Note: clustering is done on variable peptides, so merge on variable_pep
                        clusters_df = clusters_df.merge(
                            adj_DE.reset_index()[["Clean Peptide", "variable_pep", "log2FoldChange", "padj", "-log10(padj)"]],
                            left_on="Sequence", right_on="variable_pep", how="left"
                        )
                        
                        # Generate improved logos from Clean Peptide (pattern + variable part)
                        try:
                            from Library.Clustering import _logo_from_pwm
                            import logomaker as lm
                            
                            # Find the cluster column
                            cluster_col = None
                            for col in ['Cluster', 'Gn', 'cluster']:
                                if col in clusters_df.columns:
                                    cluster_col = col
                                    break
                            
                            if cluster_col and 'Clean Peptide' in clusters_df.columns:
                                # Create logos from Clean Peptides for each cluster
                                improved_logos = {}
                                for cluster_id in clusters_df[cluster_col].unique():
                                    if pd.notna(cluster_id):
                                        cluster_peptides = clusters_df[clusters_df[cluster_col] == cluster_id]['Clean Peptide'].dropna().tolist()
                                    
                                        if len(cluster_peptides) > 0:
                                            # Convert sequences to PWM for Clean Peptides
                                            max_len = max(len(seq) for seq in cluster_peptides)
                                            sequences_aligned = [seq.ljust(max_len, '-') for seq in cluster_peptides]
                                            
                                            # Create position frequency matrix
                                            aa_list = list('ACDEFGHIKLMNPQRSTVWY')
                                            pfm = np.zeros((max_len, len(aa_list)))
                                            
                                            for i, pos in enumerate(range(max_len)):
                                                pos_counts = {aa: 0 for aa in aa_list}
                                                for seq in sequences_aligned:
                                                    if pos < len(seq) and seq[pos] in aa_list:
                                                        pos_counts[seq[pos]] += 1
                                                
                                                total = sum(pos_counts.values())
                                                if total > 0:
                                                    for j, aa in enumerate(aa_list):
                                                        pfm[i, j] = pos_counts[aa] / total
                                            
                                            # Convert to DataFrame for logomaker
                                            pwm_df = pd.DataFrame(pfm, columns=aa_list)
                                            pwm_df.index = range(len(pwm_df))
                                            
                                            # Generate improved logo
                                            improved_logos[cluster_id] = _logo_from_pwm(pwm_df, f"Cluster {cluster_id}")
                                
                                # Replace original logos with improved ones
                                if improved_logos:
                                    logos = improved_logos
                                
                        except Exception as e:
                            print(f"Warning: Could not generate improved logos: {e}")
                            # Keep original logos
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
            # Ensure proper data types for Arrow serialization
            filtered_display_df = results['filtered'].head(10).copy()
            for col in filtered_display_df.columns:
                if filtered_display_df[col].dtype == 'object':
                    filtered_display_df[col] = filtered_display_df[col].astype(str)
            st.dataframe(filtered_display_df)

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
        # Ensure consistent data types for Arrow serialization
        summary_data['Count'] = summary_data['Count'].astype(str)
        st.dataframe(summary_data, use_container_width=True)
    
    with tab3:
        if not results['clusters_df'].empty:
            # Validation results if available
            if results['validation_results']:
                val_res = results['validation_results']
                st.subheader("🔬 Clustering Quality Measurements")
                
                # Primary metrics
                col_v1, col_v2, col_v3, col_v4 = st.columns(4)
                with col_v1:
                    st.metric("Optimal K", val_res['optimal_k'])
                with col_v2:
                    st.metric("Stability Score", f"{val_res['stability_score']:.3f}")
                with col_v3:
                    st.metric("Silhouette Score", f"{val_res['silhouette_score']:.3f}")
                with col_v4:
                    # Add composite score if available
                    if 'component_scores' in val_res:
                        composite = val_res['component_scores'].get('composite', 0)
                        st.metric("Composite Score", f"{composite:.3f}")
                
                # Additional quality metrics if available
                if 'stability_scores' in val_res or 'silhouette_scores' in val_res:
                    st.subheader("📈 Quality Metrics by K")
                    
                    # Create metrics comparison table
                    metrics_data = []
                    k_range = range(2, val_res['optimal_k'] + 3)  # Show range around optimal
                    
                    for k in k_range:
                        row = {'K': k}
                        if 'stability_scores' in val_res and k in val_res['stability_scores']:
                            row['Stability'] = f"{val_res['stability_scores'][k]:.3f}"
                        else:
                            row['Stability'] = 'N/A'
                            
                        if 'silhouette_scores' in val_res and k in val_res['silhouette_scores']:
                            row['Silhouette'] = f"{val_res['silhouette_scores'][k]:.3f}"  
                        else:
                            row['Silhouette'] = 'N/A'
                            
                        if 'bootstrap_scores' in val_res and k in val_res['bootstrap_scores']:
                            row['Bootstrap Stability'] = f"{val_res['bootstrap_scores'][k]:.3f}"
                        else:
                            row['Bootstrap Stability'] = 'N/A'
                            
                        # Mark optimal K
                        if k == val_res['optimal_k']:
                            row['K'] = f"{k} ⭐"
                            
                        metrics_data.append(row)
                    
                    if metrics_data:
                        metrics_df = pd.DataFrame(metrics_data)
                        st.dataframe(metrics_df, use_container_width=True, hide_index=True)
                
                # Quality interpretation with more detail
                st.subheader("🎯 Quality Assessment")
                
                col_q1, col_q2 = st.columns(2)
                
                with col_q1:
                    silhouette_score = val_res['silhouette_score']
                    st.write("**Cluster Separation Quality:**")
                    if silhouette_score > 0.7:
                        st.success("🎯 Excellent - Very well separated clusters")
                    elif silhouette_score > 0.5:
                        st.info("✅ Good - Well defined clusters")
                    elif silhouette_score > 0.25:
                        st.warning("⚠️ Moderate - Some cluster overlap")
                    else:
                        st.error("❌ Poor - Significant cluster overlap")
                    
                    st.write(f"Silhouette Score: **{silhouette_score:.3f}**")
                
                with col_q2:
                    stability_score = val_res['stability_score'] 
                    st.write("**Clustering Stability:**")
                    if stability_score > 0.8:
                        st.success("🏆 Highly Stable - Robust to data variations")
                    elif stability_score > 0.6:
                        st.info("✅ Stable - Reasonably robust")
                    elif stability_score > 0.4:
                        st.warning("⚠️ Moderately Stable - Some sensitivity")
                    else:
                        st.error("❌ Unstable - High sensitivity to variations")
                    
                    st.write(f"Stability Score: **{stability_score:.3f}**")
                
                # Show validation plot if available
                if 'validation_plot' in results['figures'] and results['figures']['validation_plot']:
                    st.subheader("📊 Validation Metrics Plot")
                    st.pyplot(results['figures']['validation_plot'])
            else:
                # Calculate basic quality metrics for fixed K clustering
                st.subheader("🔬 Clustering Quality Measurements")
                
                try:
                    from sklearn.metrics import silhouette_score
                    from Library.consensus_clustering import ConsensusClusteringValidator
                    
                    # Get cluster assignments and sequences
                    sequences = results['clusters_df']['Sequence'].tolist()
                    # Check for different possible cluster column names
                    cluster_col = None
                    for col in ['Cluster', 'Gn', 'cluster']:
                        if col in results['clusters_df'].columns:
                            cluster_col = col
                            break
                    
                    if cluster_col is None:
                        st.info("No cluster assignments found for quality calculation")
                    else:
                        clusters = results['clusters_df'][cluster_col].tolist()
                        
                        # Calculate silhouette score
                        if len(set(clusters)) > 1:  # Need at least 2 clusters
                            validator = ConsensusClusteringValidator()
                            encoded_seqs = validator.encode_sequences(sequences)
                            sil_score = silhouette_score(encoded_seqs, clusters)
                        
                            # Basic metrics display
                            col_m1, col_m2, col_m3 = st.columns(3)
                            with col_m1:
                                st.metric("Number of Clusters", len(set(clusters)))
                            with col_m2:
                                st.metric("Silhouette Score", f"{sil_score:.3f}")
                            with col_m3:
                                st.metric("Sequences Clustered", len(sequences))
                            
                            # Quality interpretation
                            st.subheader("🎯 Quality Assessment")
                            if sil_score > 0.7:
                                st.success("🎯 Excellent cluster separation")
                            elif sil_score > 0.5:
                                st.info("✅ Good cluster separation") 
                            elif sil_score > 0.25:
                                st.warning("⚠️ Moderate cluster separation")
                            else:
                                st.error("❌ Poor cluster separation")
                            
                            st.write(f"**Silhouette Score**: {sil_score:.3f}")
                        else:
                            st.info("Only one cluster found - quality metrics require multiple clusters")
                        
                except Exception as e:
                    st.warning(f"Could not calculate quality metrics: {e}")
            
            # Sequence logos
            st.subheader("🎨 Sequence Logos")
            if results['logos']:
                logo_cols = st.columns(min(len(results['logos']), 3))
                for i, (k, fig) in enumerate(sorted(results['logos'].items())):
                    with logo_cols[i % 3]:
                        st.write(f"**Cluster {k}**")
                        st.pyplot(fig)
            
            # Box plots for differential expression metrics by cluster
            if 'log2FoldChange' in results['clusters_df'].columns and 'padj' in results['clusters_df'].columns:
                st.subheader("📈 Differential Expression by Cluster")
                
                try:
                    import seaborn as sns
                    import matplotlib.pyplot as plt
                    
                    # Find the cluster column
                    cluster_col = None
                    for col in ['Cluster', 'Gn', 'cluster']:
                        if col in results['clusters_df'].columns:
                            cluster_col = col
                            break
                    
                    if cluster_col is None:
                        st.info("No cluster assignments found for box plots")
                    else:
                        # Prepare data for plotting
                        plot_columns = [cluster_col, 'log2FoldChange', 'padj', '-log10(padj)']
                        available_columns = [col for col in plot_columns if col in results['clusters_df'].columns]
                        
                        plot_df = results['clusters_df'][available_columns].copy()
                        plot_df = plot_df.dropna()
                    
                        if len(plot_df) > 0:
                            fig, axes = plt.subplots(1, 3, figsize=(15, 5))
                            
                            # Log2 Fold Change boxplot
                            if 'log2FoldChange' in plot_df.columns:
                                sns.boxplot(data=plot_df, x=cluster_col, y='log2FoldChange', ax=axes[0])
                                axes[0].set_title('Log2 Fold Change by Cluster')
                                axes[0].set_xlabel('Cluster')
                                axes[0].set_ylabel('Log2 Fold Change')
                            
                            # -log10(padj) boxplot
                            if '-log10(padj)' in plot_df.columns:
                                sns.boxplot(data=plot_df, x=cluster_col, y='-log10(padj)', ax=axes[1])
                                axes[1].set_title('-Log10(Adjusted P-value) by Cluster')
                                axes[1].set_xlabel('Cluster')
                                axes[1].set_ylabel('-Log10(Adjusted P-value)')
                            
                            # padj boxplot (log scale)
                            if 'padj' in plot_df.columns:
                                sns.boxplot(data=plot_df, x=cluster_col, y='padj', ax=axes[2])
                                axes[2].set_title('Adjusted P-value by Cluster')
                                axes[2].set_xlabel('Cluster')
                                axes[2].set_ylabel('Adjusted P-value')
                                axes[2].set_yscale('log')
                            
                            plt.tight_layout()
                            st.pyplot(fig)
                            plt.close()
                        else:
                            st.info("No data available for differential expression box plots")
                        
                except Exception as e:
                    st.warning(f"Could not generate box plots: {e}")
            
            # Cluster data table
            st.subheader("📊 Cluster Assignments")
            # Ensure proper data types for Arrow serialization
            cluster_display_df = results['clusters_df'].head(20).copy()
            for col in cluster_display_df.columns:
                if cluster_display_df[col].dtype == 'object':
                    cluster_display_df[col] = cluster_display_df[col].astype(str)
            st.dataframe(cluster_display_df, use_container_width=True)
        else:
            st.warning("No upregulated peptides found for clustering analysis.")
    
    with tab4:
        st.subheader("📋 Complete Analysis Summary")
        
        # Analysis parameters used
        st.write("**Parameters Used:**")
        params = st.session_state['analysis_params']
        param_df = pd.DataFrame([
            ["Peptide Pattern", str(params['regex'])],
            ["CPM Threshold", str(params['cpm_threshold'])],
            ["Min Samples", str(params['min_count'])],
            ["Clustering Method", "Advanced" if params['use_advanced_clustering'] else "Basic"],
            ["Validation", "Enabled" if params['use_consensus_validation'] else "Disabled"],
            ["Motif Length", str(params['motif_length'])]
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