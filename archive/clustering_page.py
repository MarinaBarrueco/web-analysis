"""
Standalone Clustering Page for Advanced Gibbs Clustering
This provides a dedicated interface for running the advanced Gibbs clustering algorithm
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Library.GibbsClusterAdvanced import gibbs_cluster_advanced, GibbsClusterAdvanced
from io import StringIO
import tempfile
import zipfile
import os
from pathlib import Path

def clustering_page():
    st.set_page_config(page_title="Advanced Gibbs Clustering", layout="wide")
    st.title("🧬 Advanced Gibbs Clustering")
    st.markdown("Based on GibbsCluster-2.0 algorithm for peptide sequence clustering")
    
    # Sidebar for parameters
    st.sidebar.header("Clustering Parameters")
    
    # Input method
    input_method = st.sidebar.selectbox(
        "Input Method",
        ["Upload File", "Paste Sequences", "Use Example Data"]
    )
    
    # Load sequences based on input method
    sequences = None
    
    if input_method == "Upload File":
        uploaded_file = st.file_uploader(
            "Upload peptide sequences", 
            type=["csv", "txt", "fasta"],
            help="CSV with 'Sequence' column, text file with one sequence per line, or FASTA format"
        )
        
        if uploaded_file is not None:
            file_type = uploaded_file.name.split('.')[-1].lower()
            
            if file_type == 'csv':
                df = pd.read_csv(uploaded_file)
                if 'Sequence' in df.columns:
                    sequences = df['Sequence'].tolist()
                elif 'sequence' in df.columns:
                    sequences = df['sequence'].tolist()
                else:
                    sequences = df.iloc[:, 0].tolist()
            
            elif file_type == 'fasta':
                content = uploaded_file.read().decode('utf-8')
                sequences = []
                for line in content.split('\n'):
                    if not line.startswith('>') and line.strip():
                        sequences.append(line.strip())
            
            else:  # txt
                content = uploaded_file.read().decode('utf-8')
                sequences = [line.strip() for line in content.split('\n') if line.strip()]
                
    elif input_method == "Paste Sequences":
        sequence_text = st.text_area(
            "Paste sequences (one per line)",
            height=200,
            placeholder="ACDEFGHIK\nLMNPQRSTU\nVWXYZABCD\n..."
        )
        
        if sequence_text:
            sequences = [line.strip() for line in sequence_text.split('\n') if line.strip()]
    
    else:  # Example data
        sequences = [
            "ACDEFGHIK", "ACDEFGHIL", "ACDEFGHIM", "ACDEFGHIN",
            "LMNPQRSTU", "LMNPQRSTV", "LMNPQRSTW", "LMNPQRSTY",
            "VWXYZABCD", "VWXYZABCE", "VWXYZABCF", "VWXYZABCG",
            "RSTUVWXYZ", "RSTUVWXYA", "RSTUVWXYB", "RSTUVWXYC"
        ]
        st.info("Using example peptide sequences")
    
    if sequences:
        st.success(f"Loaded {len(sequences)} sequences")
        
        # Show sequence preview
        with st.expander("Preview Sequences"):
            preview_df = pd.DataFrame({'Sequence': sequences[:20]})  # Show first 20
            st.dataframe(preview_df)
            if len(sequences) > 20:
                st.write(f"... and {len(sequences) - 20} more sequences")
        
        # Clustering parameters in sidebar
        st.sidebar.subheader("Basic Parameters")
        
        cluster_config = st.sidebar.selectbox(
            "Cluster Configuration", 
            ["Single Number", "Range"],
            help="Single number for fixed clusters, range to test multiple"
        )
        
        if cluster_config == "Single Number":
            num_clusters = st.sidebar.number_input("Number of Clusters", min_value=2, max_value=10, value=4)
        else:
            min_clusters = st.sidebar.number_input("Minimum Clusters", min_value=1, max_value=8, value=2)
            max_clusters = st.sidebar.number_input("Maximum Clusters", min_value=2, max_value=10, value=5)
            num_clusters = (min_clusters, max_clusters)
        
        # Determine motif length from sequences
        seq_lengths = [len(seq) for seq in sequences]
        min_length = min(seq_lengths)
        max_length = max(seq_lengths)
        
        if min_length == max_length:
            motif_length = st.sidebar.number_input(
                "Motif Length", 
                min_value=6, 
                max_value=min_length, 
                value=min(min_length, 9),
                help=f"All sequences have length {min_length}"
            )
        else:
            motif_length = st.sidebar.number_input(
                "Motif Length", 
                min_value=6, 
                max_value=min_length, 
                value=min(min_length, 9),
                help=f"Sequence lengths vary from {min_length} to {max_length}"
            )
        
        num_seeds = st.sidebar.number_input(
            "Number of Seeds", 
            min_value=1, 
            max_value=10, 
            value=3,
            help="Multiple random starts for better convergence"
        )
        
        # Advanced parameters in sidebar expander
        with st.sidebar.expander("Advanced Parameters"):
            iterations = st.number_input("Iterations per Temperature", min_value=5, max_value=100, value=10)
            temperature_start = st.number_input("Starting Temperature", min_value=0.1, max_value=5.0, value=1.5)
            temperature_steps = st.number_input("Temperature Steps", min_value=5, max_value=50, value=20)
            lambda_penalty = st.number_input("Lambda Penalty", min_value=0.0, max_value=2.0, value=0.8)
            sigma_weight = st.number_input("Sigma Weight", min_value=0.1, max_value=20.0, value=5.0)
            
            use_trash_cluster = st.checkbox("Use Trash Cluster", value=False)
            trash_threshold = st.number_input("Trash Threshold", min_value=0.0, max_value=1.0, value=0.0) if use_trash_cluster else 0.0
            
            sequence_weighting = st.selectbox(
                "Sequence Weighting",
                options=[0, 1, 2],
                format_func=lambda x: {0: "1/ns", 1: "clustering", 2: "none"}[x],
                index=0
            )
            
            background_model = st.selectbox(
                "Background Model",
                options=[0, 1, 2],
                format_func=lambda x: {0: "flat", 1: "pre-calculated", 2: "from data"}[x],
                index=1
            )
        
        # Run clustering button
        if st.button("🚀 Run Advanced Gibbs Clustering", type="primary"):
            
            with st.spinner("Running clustering analysis..."):
                try:
                    # Progress tracking
                    progress_bar = st.progress(0)
                    status_text = st.empty()
                    
                    status_text.text("Initializing clustering...")
                    progress_bar.progress(10)
                    
                    # Run clustering
                    clusterer = GibbsClusterAdvanced(
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
                    
                    status_text.text("Running Gibbs sampling...")
                    progress_bar.progress(30)
                    
                    clusterer.fit(sequences)
                    
                    progress_bar.progress(80)
                    status_text.text("Generating results...")
                    
                    # Get results
                    clusters_df = clusterer.get_clusters()
                    pwms = clusterer.get_pwms()
                    logos = clusterer.get_logos()
                    summary = clusterer.get_summary()
                    
                    progress_bar.progress(100)
                    status_text.text("Clustering completed!")
                    
                    # Display results
                    st.success("✅ Clustering Analysis Complete!")
                    
                    # Results overview
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Sequences Processed", len(sequences))
                    with col2:
                        if isinstance(num_clusters, tuple):
                            best_k = max(summary['cluster_results'].keys(), 
                                       key=lambda k: summary['cluster_results'][k]['avg_kld'])
                            st.metric("Best Number of Clusters", best_k)
                        else:
                            st.metric("Number of Clusters", num_clusters)
                    with col3:
                        best_kld = max(result['avg_kld'] for result in summary['cluster_results'].values())
                        st.metric("Best KLD Score", f"{best_kld:.4f}")
                    
                    # Clustering quality analysis for ranges
                    if isinstance(num_clusters, tuple):
                        st.subheader("📊 Clustering Quality Analysis")
                        
                        cluster_nums = sorted(summary['cluster_results'].keys())
                        kld_values = [summary['cluster_results'][k]['avg_kld'] for k in cluster_nums]
                        non_empty = [summary['cluster_results'][k]['non_empty_clusters'] for k in cluster_nums]
                        
                        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
                        
                        # KLD plot
                        ax1.bar(cluster_nums, kld_values, color='steelblue', alpha=0.7)
                        ax1.set_xlabel('Number of Clusters')
                        ax1.set_ylabel('Average KLD')
                        ax1.set_title('Clustering Quality (KLD)')
                        ax1.grid(True, alpha=0.3)
                        
                        # Non-empty clusters plot
                        ax2.bar(cluster_nums, non_empty, color='orange', alpha=0.7)
                        ax2.set_xlabel('Number of Clusters')
                        ax2.set_ylabel('Non-empty Clusters')
                        ax2.set_title('Cluster Utilization')
                        ax2.grid(True, alpha=0.3)
                        
                        plt.tight_layout()
                        st.pyplot(fig)
                        
                        # Results table
                        results_data = []
                        for k in cluster_nums:
                            result = summary['cluster_results'][k]
                            results_data.append({
                                'Clusters': k,
                                'Avg KLD': f"{result['avg_kld']:.4f}",
                                'Non-empty': result['non_empty_clusters'],
                                'Sizes': ', '.join(map(str, result['cluster_sizes']))
                            })
                        
                        st.subheader("📋 Detailed Results")
                        st.dataframe(pd.DataFrame(results_data))
                    
                    # Sequence logos
                    st.subheader("🎨 Sequence Logos")
                    
                    # Display logos in columns
                    num_logos = len(logos)
                    cols_per_row = 2
                    
                    for i in range(0, num_logos, cols_per_row):
                        cols = st.columns(cols_per_row)
                        for j in range(cols_per_row):
                            if i + j < num_logos:
                                cluster_id = list(logos.keys())[i + j]
                                with cols[j]:
                                    st.pyplot(logos[cluster_id])
                    
                    # Cluster assignments
                    st.subheader("📋 Cluster Assignments")
                    
                    # Add cluster statistics
                    cluster_stats = clusters_df['Cluster'].value_counts().sort_index()
                    
                    col1, col2 = st.columns([1, 2])
                    with col1:
                        st.write("**Cluster Sizes:**")
                        for cluster, size in cluster_stats.items():
                            st.write(f"Cluster {cluster}: {size} sequences")
                    
                    with col2:
                        # Show cluster assignments table
                        st.dataframe(clusters_df)
                    
                    # Download options
                    st.subheader("💾 Download Results")
                    
                    col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        # Download cluster assignments
                        csv_data = clusters_df.to_csv(index=False)
                        st.download_button(
                            label="📄 Download Cluster Assignments",
                            data=csv_data,
                            file_name="cluster_assignments.csv",
                            mime="text/csv"
                        )
                    
                    with col2:
                        # Download PWMs
                        pwm_data = ""
                        for i, pwm in enumerate(pwms):
                            pwm_data += f"# Cluster {i+1} PWM\n"
                            pwm_data += pwm.to_csv()
                            pwm_data += "\n"
                        
                        st.download_button(
                            label="📊 Download PWMs",
                            data=pwm_data,
                            file_name="position_weight_matrices.csv",
                            mime="text/csv"
                        )
                    
                    with col3:
                        # Download summary
                        import json
                        summary_json = json.dumps(summary, indent=2)
                        
                        st.download_button(
                            label="📈 Download Summary",
                            data=summary_json,
                            file_name="clustering_summary.json",
                            mime="application/json"
                        )
                    
                    # Create zip with all results
                    with tempfile.TemporaryDirectory() as temp_dir:
                        temp_path = Path(temp_dir)
                        
                        # Save all files
                        clusters_df.to_csv(temp_path / "cluster_assignments.csv", index=False)
                        
                        for i, pwm in enumerate(pwms):
                            pwm.to_csv(temp_path / f"cluster_{i+1}_pwm.csv")
                        
                        for cluster_id, fig in logos.items():
                            fig.savefig(temp_path / f"cluster_{cluster_id}_logo.png", 
                                      dpi=300, bbox_inches='tight')
                        
                        with open(temp_path / "summary.json", 'w') as f:
                            json.dump(summary, f, indent=2)
                        
                        # Create zip
                        zip_path = temp_path / "clustering_results.zip"
                        with zipfile.ZipFile(zip_path, 'w') as zipf:
                            for file_path in temp_path.glob("*"):
                                if file_path != zip_path:
                                    zipf.write(file_path, file_path.name)
                        
                        with open(zip_path, 'rb') as f:
                            zip_data = f.read()
                        
                        st.download_button(
                            label="📦 Download All Results (ZIP)",
                            data=zip_data,
                            file_name="clustering_results.zip",
                            mime="application/zip"
                        )
                
                except Exception as e:
                    st.error(f"❌ Clustering failed: {str(e)}")
                    st.exception(e)
    
    else:
        st.info("👆 Please provide peptide sequences to start clustering analysis")

if __name__ == "__main__":
    clustering_page()