"""
Consensus Clustering and Validation Methods for Peptide Motif Discovery

This module implements advanced clustering validation techniques including:
- Consensus clustering across multiple seeds
- Silhouette analysis for cluster quality
- Stability-based cluster number selection
- Bootstrap validation
- Adjusted Rand Index comparison
- Cluster quality metrics

Based on:
- Monti et al. (2003) "Consensus Clustering: A Resampling-Based Method for Class Discovery"
- Rousseeuw (1987) "Silhouettes: a Graphical Aid to the Interpretation and Validation"
- Strehl & Ghosh (2002) "Cluster Ensembles: A Knowledge Reuse Framework"
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict, Tuple, Optional, Union, Any
from collections import defaultdict
import itertools
from sklearn.metrics import adjusted_rand_score, silhouette_score, calinski_harabasz_score
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.stats import zscore
import warnings
warnings.filterwarnings('ignore')

class ConsensusClusteringValidator:
    """
    Advanced clustering validation using consensus methods for peptide sequences
    """
    
    def __init__(self, random_state: int = 42):
        """
        Initialize consensus clustering validator
        
        Parameters:
        -----------
        random_state : int
            Random seed for reproducibility
        """
        self.random_state = random_state
        np.random.seed(random_state)
        self.consensus_matrices = {}
        self.stability_scores = {}
        self.silhouette_scores = {}
        self.clustering_results = {}
        
    def encode_sequences(self, sequences: List[str]) -> np.ndarray:
        """
        Encode peptide sequences into numerical feature matrix
        
        Parameters:
        -----------
        sequences : List[str]
            List of peptide sequences
            
        Returns:
        --------
        np.ndarray : Encoded feature matrix
        """
        AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
        AA_IDX = {aa: i for i, aa in enumerate(AA_ALPHABET)}
        
        if not sequences:
            raise ValueError("No sequences provided")
        
        # Ensure all sequences have the same length for this analysis
        seq_lengths = [len(seq) for seq in sequences]
        if len(set(seq_lengths)) > 1:
            # Pad shorter sequences or trim longer ones to median length
            median_len = int(np.median(seq_lengths))
            sequences = [seq[:median_len].ljust(median_len, 'A') for seq in sequences]
        
        # One-hot encoding
        seq_len = len(sequences[0])
        n_seqs = len(sequences)
        encoded = np.zeros((n_seqs, seq_len * 20))
        
        for i, seq in enumerate(sequences):
            for j, aa in enumerate(seq):
                if aa in AA_IDX:
                    encoded[i, j * 20 + AA_IDX[aa]] = 1
                    
        return encoded
    
    def calculate_sequence_similarity_matrix(self, sequences: List[str]) -> np.ndarray:
        """
        Calculate pairwise sequence similarity matrix using amino acid properties
        
        Parameters:
        -----------
        sequences : List[str]
            List of peptide sequences
            
        Returns:
        --------
        np.ndarray : Similarity matrix
        """
        # Amino acid property-based scoring
        AA_PROPERTIES = {
            'A': [1.8, 0, 0, 0],   # Hydrophobicity, Charge, Size, Aromatic
            'C': [2.5, 0, 1, 0],
            'D': [-3.5, -1, 1, 0],
            'E': [-3.5, -1, 2, 0],
            'F': [2.8, 0, 2, 1],
            'G': [-0.4, 0, 0, 0],
            'H': [-3.2, 1, 2, 1],
            'I': [4.5, 0, 2, 0],
            'K': [-3.9, 1, 2, 0],
            'L': [3.8, 0, 2, 0],
            'M': [1.9, 0, 2, 0],
            'N': [-3.5, 0, 1, 0],
            'P': [-1.6, 0, 1, 0],
            'Q': [-3.5, 0, 2, 0],
            'R': [-4.5, 1, 3, 0],
            'S': [-0.8, 0, 1, 0],
            'T': [-0.7, 0, 1, 0],
            'V': [4.2, 0, 1, 0],
            'W': [-0.9, 0, 3, 1],
            'Y': [-1.3, 0, 3, 1]
        }
        
        n_seqs = len(sequences)
        similarity_matrix = np.zeros((n_seqs, n_seqs))
        
        for i in range(n_seqs):
            for j in range(i, n_seqs):
                if i == j:
                    similarity_matrix[i, j] = 1.0
                else:
                    seq1, seq2 = sequences[i], sequences[j]
                    min_len = min(len(seq1), len(seq2))
                    
                    # Calculate position-wise similarity
                    similarities = []
                    for k in range(min_len):
                        aa1, aa2 = seq1[k], seq2[k]
                        if aa1 in AA_PROPERTIES and aa2 in AA_PROPERTIES:
                            props1 = np.array(AA_PROPERTIES[aa1])
                            props2 = np.array(AA_PROPERTIES[aa2])
                            # Cosine similarity between property vectors
                            sim = np.dot(props1, props2) / (np.linalg.norm(props1) * np.linalg.norm(props2))
                            similarities.append(sim)
                    
                    avg_similarity = np.mean(similarities) if similarities else 0
                    similarity_matrix[i, j] = similarity_matrix[j, i] = avg_similarity
        
        return similarity_matrix
    
    def run_consensus_clustering(self, 
                                sequences: List[str],
                                k_range: range = range(2, 8),
                                n_iterations: int = 100,
                                sample_fraction: float = 0.8) -> Dict[int, np.ndarray]:
        """
        Run consensus clustering across multiple values of k
        
        Parameters:
        -----------
        sequences : List[str]
            List of peptide sequences
        k_range : range
            Range of cluster numbers to test
        n_iterations : int
            Number of consensus iterations
        sample_fraction : float
            Fraction of samples to use in each iteration
            
        Returns:
        --------
        Dict[int, np.ndarray] : Consensus matrices for each k
        """
        print(f"Running consensus clustering for k in {list(k_range)} with {n_iterations} iterations...")
        
        n_seqs = len(sequences)
        encoded_seqs = self.encode_sequences(sequences)
        
        consensus_matrices = {}
        
        for k in k_range:
            print(f"  Processing k={k}...")
            
            # Initialize consensus matrix
            consensus_matrix = np.zeros((n_seqs, n_seqs))
            indicator_matrix = np.zeros((n_seqs, n_seqs))
            
            for iteration in range(n_iterations):
                # Subsample sequences
                n_sample = int(sample_fraction * n_seqs)
                sample_indices = np.random.choice(n_seqs, n_sample, replace=False)
                
                # Cluster the subsample
                sample_data = encoded_seqs[sample_indices]
                
                # Use KMeans for initial clustering (can be replaced with Gibbs results)
                try:
                    kmeans = KMeans(n_clusters=k, random_state=self.random_state + iteration, n_init=10)
                    cluster_labels = kmeans.fit_predict(sample_data)
                    
                    # Update consensus and indicator matrices
                    for i in range(len(sample_indices)):
                        for j in range(i + 1, len(sample_indices)):
                            idx_i, idx_j = sample_indices[i], sample_indices[j]
                            indicator_matrix[idx_i, idx_j] += 1
                            indicator_matrix[idx_j, idx_i] += 1
                            
                            if cluster_labels[i] == cluster_labels[j]:
                                consensus_matrix[idx_i, idx_j] += 1
                                consensus_matrix[idx_j, idx_i] += 1
                                
                except Exception as e:
                    print(f"    Warning: Iteration {iteration} failed for k={k}: {e}")
                    continue
            
            # Normalize consensus matrix
            with np.errstate(divide='ignore', invalid='ignore'):
                consensus_matrix = np.divide(consensus_matrix, indicator_matrix, 
                                           out=np.zeros_like(consensus_matrix), 
                                           where=indicator_matrix!=0)
            
            # Fill diagonal with 1s
            np.fill_diagonal(consensus_matrix, 1.0)
            
            consensus_matrices[k] = consensus_matrix
            self.consensus_matrices[k] = consensus_matrix
            
        return consensus_matrices
    
    def calculate_stability_scores(self, consensus_matrices: Dict[int, np.ndarray]) -> Dict[int, float]:
        """
        Calculate stability scores for each k based on consensus matrices
        
        Parameters:
        -----------
        consensus_matrices : Dict[int, np.ndarray]
            Consensus matrices for different k values
            
        Returns:
        --------
        Dict[int, float] : Stability scores for each k
        """
        stability_scores = {}
        
        for k, consensus_matrix in consensus_matrices.items():
            # Calculate area under CDF of consensus values
            consensus_values = consensus_matrix[np.triu_indices_from(consensus_matrix, k=1)]
            consensus_values = consensus_values[~np.isnan(consensus_values)]
            
            if len(consensus_values) == 0:
                stability_scores[k] = 0.0
                continue
            
            # Calculate cumulative distribution
            sorted_values = np.sort(consensus_values)
            n_values = len(sorted_values)
            cdf_x = np.linspace(0, 1, 100)
            cdf_y = np.searchsorted(sorted_values, cdf_x * np.max(sorted_values)) / n_values
            
            # Area under CDF (stability score)
            stability_score = np.trapz(cdf_y, cdf_x)
            stability_scores[k] = stability_score
            
        self.stability_scores = stability_scores
        return stability_scores
    
    def calculate_silhouette_scores(self, 
                                   sequences: List[str], 
                                   cluster_assignments: Dict[int, np.ndarray]) -> Dict[int, float]:
        """
        Calculate silhouette scores for different clustering solutions
        
        Parameters:
        -----------
        sequences : List[str]
            List of peptide sequences
        cluster_assignments : Dict[int, np.ndarray]
            Cluster assignments for each k
            
        Returns:
        --------
        Dict[int, float] : Silhouette scores for each k
        """
        encoded_seqs = self.encode_sequences(sequences)
        silhouette_scores = {}
        
        for k, labels in cluster_assignments.items():
            if len(np.unique(labels)) > 1:  # Need at least 2 clusters for silhouette
                try:
                    score = silhouette_score(encoded_seqs, labels)
                    silhouette_scores[k] = score
                except Exception as e:
                    print(f"Warning: Could not calculate silhouette score for k={k}: {e}")
                    silhouette_scores[k] = -1.0
            else:
                silhouette_scores[k] = -1.0
                
        self.silhouette_scores = silhouette_scores
        return silhouette_scores
    
    def gibbs_clustering_to_assignments(self, gibbs_results: Dict[int, Any]) -> Dict[int, np.ndarray]:
        """
        Convert Gibbs clustering results to standard cluster assignment format
        
        Parameters:
        -----------
        gibbs_results : Dict[int, Any]
            Results from Gibbs clustering for different k values
            
        Returns:
        --------
        Dict[int, np.ndarray] : Cluster assignments for each k
        """
        cluster_assignments = {}
        
        for k, result in gibbs_results.items():
            if hasattr(result, 'get_clusters'):
                # Advanced Gibbs clusterer
                clusters_df = result.get_clusters()
                if 'Cluster' in clusters_df.columns:
                    assignments = clusters_df['Cluster'].values - 1  # Convert to 0-based
                elif 'Gn' in clusters_df.columns:
                    assignments = clusters_df['Gn'].values - 1  # Convert to 0-based
                else:
                    assignments = np.zeros(len(clusters_df))
            elif isinstance(result, tuple) and len(result) >= 1:
                # Basic Gibbs clusterer returns (clusters_df, pwms, logos)
                clusters_df = result[0]
                if 'Gn' in clusters_df.columns:
                    assignments = clusters_df['Gn'].values - 1  # Convert to 0-based
                else:
                    assignments = np.zeros(len(clusters_df))
            else:
                # Fallback
                assignments = np.zeros(len(gibbs_results[list(gibbs_results.keys())[0]]))
                
            cluster_assignments[k] = assignments
            
        return cluster_assignments
    
    def calculate_adjusted_rand_index(self, 
                                     assignments1: Dict[int, np.ndarray], 
                                     assignments2: Dict[int, np.ndarray]) -> Dict[int, float]:
        """
        Calculate Adjusted Rand Index between two clustering solutions
        
        Parameters:
        -----------
        assignments1, assignments2 : Dict[int, np.ndarray]
            Cluster assignments for comparison
            
        Returns:
        --------
        Dict[int, float] : ARI scores for each k
        """
        ari_scores = {}
        
        common_ks = set(assignments1.keys()) & set(assignments2.keys())
        
        for k in common_ks:
            try:
                ari = adjusted_rand_score(assignments1[k], assignments2[k])
                ari_scores[k] = ari
            except Exception as e:
                print(f"Warning: Could not calculate ARI for k={k}: {e}")
                ari_scores[k] = 0.0
                
        return ari_scores
    
    def bootstrap_clustering_stability(self, 
                                      sequences: List[str], 
                                      clustering_func: callable,
                                      k_range: range = range(2, 8),
                                      n_bootstrap: int = 50,
                                      sample_fraction: float = 0.8) -> Dict[int, float]:
        """
        Assess clustering stability using bootstrap resampling
        
        Parameters:
        -----------
        sequences : List[str]
            List of peptide sequences
        clustering_func : callable
            Function that performs clustering
        k_range : range
            Range of cluster numbers to test
        n_bootstrap : int
            Number of bootstrap iterations
        sample_fraction : float
            Fraction of samples in each bootstrap
            
        Returns:
        --------
        Dict[int, float] : Bootstrap stability scores for each k
        """
        print(f"Running bootstrap stability analysis with {n_bootstrap} iterations...")
        
        n_seqs = len(sequences)
        bootstrap_stability = {}
        
        for k in k_range:
            print(f"  Processing k={k}...")
            ari_scores = []
            
            # Get reference clustering on full dataset
            try:
                ref_result = clustering_func(sequences, k)
                ref_assignments = self._extract_assignments(ref_result)
            except Exception as e:
                print(f"    Warning: Reference clustering failed for k={k}: {e}")
                bootstrap_stability[k] = 0.0
                continue
            
            for i in range(n_bootstrap):
                try:
                    # Bootstrap sample
                    sample_size = int(sample_fraction * n_seqs)
                    boot_indices = np.random.choice(n_seqs, sample_size, replace=True)
                    boot_sequences = [sequences[idx] for idx in boot_indices]
                    
                    # Cluster bootstrap sample
                    boot_result = clustering_func(boot_sequences, k)
                    boot_assignments = self._extract_assignments(boot_result)
                    
                    # Calculate ARI with reference (only for overlapping sequences)
                    unique_indices, inverse_indices = np.unique(boot_indices, return_inverse=True)
                    if len(unique_indices) > 1:
                        ref_subset = ref_assignments[unique_indices]
                        boot_subset = [boot_assignments[inverse_indices[j]] for j in range(len(unique_indices))]
                        
                        if len(set(ref_subset)) > 1 and len(set(boot_subset)) > 1:
                            ari = adjusted_rand_score(ref_subset, boot_subset)
                            ari_scores.append(ari)
                            
                except Exception as e:
                    print(f"    Warning: Bootstrap iteration {i} failed for k={k}: {e}")
                    continue
            
            # Average ARI as stability measure
            bootstrap_stability[k] = np.mean(ari_scores) if ari_scores else 0.0
            
        return bootstrap_stability
    
    def _extract_assignments(self, clustering_result) -> np.ndarray:
        """Helper function to extract cluster assignments from various result formats"""
        if hasattr(clustering_result, 'get_clusters'):
            clusters_df = clustering_result.get_clusters()
            if 'Cluster' in clusters_df.columns:
                return clusters_df['Cluster'].values - 1
            elif 'Gn' in clusters_df.columns:
                return clusters_df['Gn'].values - 1
        elif isinstance(clustering_result, tuple):
            clusters_df = clustering_result[0]
            if 'Gn' in clusters_df.columns:
                return clusters_df['Gn'].values - 1
        
        # Fallback
        return np.zeros(len(clustering_result) if hasattr(clustering_result, '__len__') else 10)
    
    def select_optimal_k(self, 
                        stability_scores: Dict[int, float],
                        silhouette_scores: Dict[int, float],
                        bootstrap_scores: Optional[Dict[int, float]] = None,
                        weights: Tuple[float, float, float] = (0.4, 0.4, 0.2)) -> Tuple[int, Dict[str, float]]:
        """
        Select optimal number of clusters using multiple criteria
        
        Parameters:
        -----------
        stability_scores : Dict[int, float]
            Consensus stability scores
        silhouette_scores : Dict[int, float]
            Silhouette scores
        bootstrap_scores : Dict[int, float], optional
            Bootstrap stability scores
        weights : Tuple[float, float, float]
            Weights for (stability, silhouette, bootstrap) scores
            
        Returns:
        --------
        Tuple[int, Dict[str, float]] : Optimal k and component scores
        """
        common_ks = set(stability_scores.keys()) & set(silhouette_scores.keys())
        if bootstrap_scores:
            common_ks &= set(bootstrap_scores.keys())
        
        if not common_ks:
            raise ValueError("No common k values across all scoring methods")
        
        # Normalize scores to [0, 1]
        def normalize_scores(scores_dict):
            values = np.array(list(scores_dict.values()))
            if np.max(values) == np.min(values):
                return {k: 0.5 for k in scores_dict.keys()}
            return {k: (v - np.min(values)) / (np.max(values) - np.min(values)) 
                   for k, v in scores_dict.items()}
        
        norm_stability = normalize_scores(stability_scores)
        norm_silhouette = normalize_scores(silhouette_scores)
        
        composite_scores = {}
        
        for k in common_ks:
            score = weights[0] * norm_stability[k] + weights[1] * norm_silhouette[k]
            
            if bootstrap_scores:
                norm_bootstrap = normalize_scores(bootstrap_scores)
                score += weights[2] * norm_bootstrap[k]
            
            composite_scores[k] = score
        
        # Select k with highest composite score
        optimal_k = max(composite_scores.keys(), key=lambda k: composite_scores[k])
        
        component_scores = {
            'stability': stability_scores[optimal_k],
            'silhouette': silhouette_scores[optimal_k],
            'composite': composite_scores[optimal_k]
        }
        
        if bootstrap_scores:
            component_scores['bootstrap'] = bootstrap_scores[optimal_k]
        
        return optimal_k, component_scores
    
    def plot_validation_results(self, 
                               stability_scores: Dict[int, float],
                               silhouette_scores: Dict[int, float],
                               bootstrap_scores: Optional[Dict[int, float]] = None,
                               optimal_k: Optional[int] = None,
                               save_path: Optional[str] = None) -> plt.Figure:
        """
        Create comprehensive validation results plot
        
        Parameters:
        -----------
        stability_scores : Dict[int, float]
            Consensus stability scores
        silhouette_scores : Dict[int, float]
            Silhouette scores  
        bootstrap_scores : Dict[int, float], optional
            Bootstrap stability scores
        optimal_k : int, optional
            Optimal k to highlight
        save_path : str, optional
            Path to save the plot
            
        Returns:
        --------
        plt.Figure : The validation results figure
        """
        n_plots = 3 if bootstrap_scores else 2
        fig, axes = plt.subplots(1, n_plots, figsize=(5 * n_plots, 5))
        if n_plots == 2:
            axes = [axes[0], axes[1]]
        
        # Plot 1: Stability scores
        k_values = sorted(stability_scores.keys())
        stability_values = [stability_scores[k] for k in k_values]
        
        axes[0].plot(k_values, stability_values, 'b-o', linewidth=2, markersize=8)
        axes[0].set_xlabel('Number of Clusters (k)')
        axes[0].set_ylabel('Consensus Stability Score')
        axes[0].set_title('Consensus Clustering Stability')
        axes[0].grid(True, alpha=0.3)
        
        if optimal_k and optimal_k in stability_scores:
            axes[0].axvline(x=optimal_k, color='red', linestyle='--', alpha=0.7, label=f'Optimal k={optimal_k}')
            axes[0].legend()
        
        # Plot 2: Silhouette scores
        silhouette_values = [silhouette_scores[k] for k in k_values if k in silhouette_scores]
        silhouette_k_values = [k for k in k_values if k in silhouette_scores]
        
        axes[1].plot(silhouette_k_values, silhouette_values, 'g-s', linewidth=2, markersize=8)
        axes[1].set_xlabel('Number of Clusters (k)')
        axes[1].set_ylabel('Silhouette Score')
        axes[1].set_title('Silhouette Analysis')
        axes[1].grid(True, alpha=0.3)
        
        if optimal_k and optimal_k in silhouette_scores:
            axes[1].axvline(x=optimal_k, color='red', linestyle='--', alpha=0.7, label=f'Optimal k={optimal_k}')
            axes[1].legend()
        
        # Plot 3: Bootstrap scores (if available)
        if bootstrap_scores:
            bootstrap_values = [bootstrap_scores[k] for k in k_values if k in bootstrap_scores]
            bootstrap_k_values = [k for k in k_values if k in bootstrap_scores]
            
            axes[2].plot(bootstrap_k_values, bootstrap_values, 'r-^', linewidth=2, markersize=8)
            axes[2].set_xlabel('Number of Clusters (k)')
            axes[2].set_ylabel('Bootstrap Stability Score')
            axes[2].set_title('Bootstrap Validation')
            axes[2].grid(True, alpha=0.3)
            
            if optimal_k and optimal_k in bootstrap_scores:
                axes[2].axvline(x=optimal_k, color='red', linestyle='--', alpha=0.7, label=f'Optimal k={optimal_k}')
                axes[2].legend()
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def plot_consensus_matrix(self, 
                             consensus_matrix: np.ndarray, 
                             k: int,
                             sequence_names: Optional[List[str]] = None,
                             save_path: Optional[str] = None) -> plt.Figure:
        """
        Plot consensus matrix heatmap
        
        Parameters:
        -----------
        consensus_matrix : np.ndarray
            Consensus matrix to plot
        k : int
            Number of clusters
        sequence_names : List[str], optional
            Names for sequences
        save_path : str, optional
            Path to save the plot
            
        Returns:
        --------
        plt.Figure : The consensus matrix figure
        """
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Create heatmap
        im = ax.imshow(consensus_matrix, cmap='RdYlBu_r', vmin=0, vmax=1)
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Consensus Value', rotation=270, labelpad=20)
        
        # Set labels if provided
        if sequence_names:
            ax.set_xticks(range(len(sequence_names)))
            ax.set_yticks(range(len(sequence_names)))
            ax.set_xticklabels(sequence_names, rotation=45, ha='right')
            ax.set_yticklabels(sequence_names)
        
        ax.set_title(f'Consensus Matrix (k={k})')
        ax.set_xlabel('Sequences')
        ax.set_ylabel('Sequences')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def generate_validation_report(self, 
                                  sequences: List[str],
                                  results: Dict[str, Any],
                                  optimal_k: int,
                                  save_path: Optional[str] = None) -> str:
        """
        Generate comprehensive validation report
        
        Parameters:
        -----------
        sequences : List[str]
            Input sequences
        results : Dict[str, Any]
            Validation results
        optimal_k : int
            Optimal number of clusters
        save_path : str, optional
            Path to save the report
            
        Returns:
        --------
        str : Validation report
        """
        report = f"""
# Consensus Clustering Validation Report

## Dataset Summary
- **Number of sequences**: {len(sequences)}
- **Sequence length range**: {min(len(s) for s in sequences)} - {max(len(s) for s in sequences)}
- **Optimal number of clusters**: {optimal_k}

## Validation Methods Applied

### 1. Consensus Clustering Analysis
- **Method**: Multiple subsampling with {results.get('n_iterations', 'N/A')} iterations
- **Sample fraction**: {results.get('sample_fraction', 'N/A')}
- **Stability score (k={optimal_k})**: {results.get('stability_scores', {}).get(optimal_k, 'N/A'):.4f}

### 2. Silhouette Analysis
- **Silhouette score (k={optimal_k})**: {results.get('silhouette_scores', {}).get(optimal_k, 'N/A'):.4f}
- **Interpretation**: 
  - > 0.7: Strong clustering structure
  - 0.5-0.7: Reasonable clustering structure  
  - 0.25-0.5: Weak clustering structure
  - < 0.25: No substantial clustering structure

### 3. Bootstrap Validation
- **Bootstrap iterations**: {results.get('n_bootstrap', 'N/A')}
- **Bootstrap stability (k={optimal_k})**: {results.get('bootstrap_scores', {}).get(optimal_k, 'N/A'):.4f}

## Cluster Quality Assessment

### Stability Across k Values
"""
        
        # Add stability scores table
        if 'stability_scores' in results:
            report += "\n| k | Stability Score | Silhouette Score |\n|---|---|---|\n"
            for k in sorted(results['stability_scores'].keys()):
                stability = results['stability_scores'][k]
                silhouette = results.get('silhouette_scores', {}).get(k, 'N/A')
                marker = " ← **Optimal**" if k == optimal_k else ""
                report += f"| {k} | {stability:.4f} | {silhouette:.4f if silhouette != 'N/A' else 'N/A'} |{marker}\n"
        
        report += f"""

## Recommendations

### Selected k = {optimal_k}
"""
        
        # Add specific recommendations based on scores
        if 'stability_scores' in results and optimal_k in results['stability_scores']:
            stability_score = results['stability_scores'][optimal_k]
            if stability_score > 0.8:
                report += "- **High stability**: The clustering solution is very robust\n"
            elif stability_score > 0.6:
                report += "- **Moderate stability**: The clustering solution is reasonably robust\n" 
            else:
                report += "- **Low stability**: Consider alternative clustering approaches or parameter tuning\n"
        
        if 'silhouette_scores' in results and optimal_k in results['silhouette_scores']:
            silhouette_score = results['silhouette_scores'][optimal_k]
            if silhouette_score > 0.7:
                report += "- **Strong cluster separation**: Well-defined clusters with good separation\n"
            elif silhouette_score > 0.5:
                report += "- **Reasonable cluster separation**: Acceptable cluster quality\n"
            else:
                report += "- **Weak cluster separation**: Consider fewer clusters or alternative methods\n"
        
        report += f"""

## Technical Details

### Methodology
- **Consensus clustering**: Based on Monti et al. (2003) resampling approach
- **Silhouette analysis**: Rousseeuw (1987) cluster validation method  
- **Bootstrap validation**: Stability assessment through resampling
- **Feature encoding**: One-hot encoding of amino acid sequences
- **Similarity metric**: Amino acid property-based similarity

### Quality Metrics
- **Stability score**: Area under consensus CDF (higher = more stable)
- **Silhouette score**: Average silhouette width (higher = better separation)
- **Bootstrap score**: Average ARI across bootstrap samples (higher = more stable)

---

*Report generated by ConsensusClusteringValidator*
"""
        
        if save_path:
            with open(save_path, 'w') as f:
                f.write(report)
        
        return report