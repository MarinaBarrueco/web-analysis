#!/usr/bin/env python3
"""
Enhanced Gibbs Clustering with Integrated Validation

Unified interface combining Gibbs clustering with consensus validation methods.
Provides automatic cluster number selection, quality assessment, and robust
clustering with comprehensive validation metrics.

Features:
- Integration of basic and advanced Gibbs clustering algorithms
- Automatic optimal cluster number selection using consensus validation
- Comprehensive clustering quality metrics and assessment
- Bootstrap validation for robustness testing
- Parameter filtering for compatibility between different clustering modes
- Unified result format with quality scores and interpretive guidance

Author: Peptide Analysis Pipeline
Version: 2.0
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import List, Dict, Tuple, Optional, Union, Any
import warnings
warnings.filterwarnings('ignore')

from .Clustering import gibbs_cluster
from .GibbsClusterAdvanced import GibbsClusterAdvanced
from .consensus_clustering import ConsensusClusteringValidator

class EnhancedGibbsClustering:
    """
    Enhanced Gibbs clustering with consensus validation and automatic parameter selection
    """
    
    def __init__(self, 
                 motif_length: int = 8,
                 use_advanced: bool = True,
                 validation_enabled: bool = True,
                 random_seed: int = 42):
        """
        Initialize enhanced Gibbs clustering
        
        Parameters:
        -----------
        motif_length : int
            Length of motifs to discover
        use_advanced : bool
            Whether to use advanced Gibbs clustering
        validation_enabled : bool
            Whether to perform consensus validation
        random_seed : int
            Random seed for reproducibility
        """
        self.motif_length = motif_length
        self.use_advanced = use_advanced
        self.validation_enabled = validation_enabled
        self.random_seed = random_seed
        
        # Initialize validator if validation is enabled
        if self.validation_enabled:
            self.validator = ConsensusClusteringValidator(random_state=random_seed)
        
        # Store results
        self.validation_results = {}
        self.clustering_results = {}
        self.optimal_k = None
        self.final_clusters = None
        
    def run_clustering_with_validation(self,
                                      sequences: List[str],
                                      k_range: Union[range, Tuple[int, int]] = (2, 8),
                                      n_consensus_iterations: int = 50,
                                      n_bootstrap_iterations: int = 30,
                                      **clustering_kwargs) -> Dict[str, Any]:
        """
        Run Gibbs clustering with comprehensive validation
        
        Parameters:
        -----------
        sequences : List[str]
            Input peptide sequences
        k_range : Union[range, Tuple[int, int]]
            Range of cluster numbers to test
        n_consensus_iterations : int
            Number of consensus clustering iterations
        n_bootstrap_iterations : int
            Number of bootstrap validation iterations
        **clustering_kwargs : dict
            Additional parameters for Gibbs clustering
            
        Returns:
        --------
        Dict[str, Any] : Comprehensive clustering and validation results
        """
        print(f"🧬 Enhanced Gibbs Clustering with Consensus Validation")
        print(f"   Sequences: {len(sequences)}")
        print(f"   Motif length: {self.motif_length}")
        print(f"   K range: {k_range}")
        print(f"   Advanced clustering: {self.use_advanced}")
        print(f"   Validation enabled: {self.validation_enabled}")
        
        # Convert k_range to range object if needed
        if isinstance(k_range, tuple):
            k_range = range(k_range[0], k_range[1] + 1)
        
        # Step 1: Run Gibbs clustering for each k
        print(f"\n🔄 Step 1: Running Gibbs clustering for k in {list(k_range)}...")
        gibbs_results = self._run_gibbs_clustering_range(sequences, k_range, **clustering_kwargs)
        
        if not self.validation_enabled:
            # If validation is disabled, return results with middle k
            middle_k = list(k_range)[len(k_range) // 2]
            self.optimal_k = middle_k
            self.final_clusters = gibbs_results[middle_k]
            return {
                'clustering_results': gibbs_results,
                'optimal_k': middle_k,
                'final_clusters': gibbs_results[middle_k],
                'validation_enabled': False
            }
        
        # Step 2: Consensus clustering validation
        print(f"\n📊 Step 2: Consensus clustering validation...")
        consensus_matrices = self.validator.run_consensus_clustering(
            sequences, 
            k_range=k_range, 
            n_iterations=n_consensus_iterations,
            sample_fraction=0.8
        )
        
        # Step 3: Calculate stability scores
        print(f"\n📈 Step 3: Calculating stability scores...")
        stability_scores = self.validator.calculate_stability_scores(consensus_matrices)
        
        # Step 4: Extract cluster assignments and calculate silhouette scores
        print(f"\n🎯 Step 4: Silhouette analysis...")
        cluster_assignments = self.validator.gibbs_clustering_to_assignments(gibbs_results)
        silhouette_scores = self.validator.calculate_silhouette_scores(sequences, cluster_assignments)
        
        # Step 5: Bootstrap validation (optional - can be slow)
        bootstrap_scores = None
        if n_bootstrap_iterations > 0:
            print(f"\n🔄 Step 5: Bootstrap validation ({n_bootstrap_iterations} iterations)...")
            try:
                clustering_func = lambda seqs, k: self._run_single_gibbs_clustering(
                    seqs, k, **clustering_kwargs
                )
                bootstrap_scores = self.validator.bootstrap_clustering_stability(
                    sequences, 
                    clustering_func,
                    k_range=k_range,
                    n_bootstrap=n_bootstrap_iterations
                )
            except Exception as e:
                print(f"   Warning: Bootstrap validation failed: {e}")
                bootstrap_scores = None
        
        # Step 6: Select optimal k
        print(f"\n🎖️ Step 6: Selecting optimal k...")
        optimal_k, component_scores = self.validator.select_optimal_k(
            stability_scores, 
            silhouette_scores, 
            bootstrap_scores
        )
        
        print(f"   Selected k = {optimal_k}")
        print(f"   Stability score: {component_scores['stability']:.4f}")
        print(f"   Silhouette score: {component_scores['silhouette']:.4f}")
        if 'bootstrap' in component_scores:
            print(f"   Bootstrap score: {component_scores['bootstrap']:.4f}")
        
        # Store results
        self.optimal_k = optimal_k
        self.final_clusters = gibbs_results[optimal_k]
        self.validation_results = {
            'stability_scores': stability_scores,
            'silhouette_scores': silhouette_scores,
            'bootstrap_scores': bootstrap_scores,
            'consensus_matrices': consensus_matrices,
            'component_scores': component_scores,
            'n_consensus_iterations': n_consensus_iterations,
            'n_bootstrap_iterations': n_bootstrap_iterations,
            'sample_fraction': 0.8
        }
        self.clustering_results = gibbs_results
        
        # Create comprehensive results
        results = {
            'clustering_results': gibbs_results,
            'validation_results': self.validation_results,
            'optimal_k': optimal_k,
            'final_clusters': gibbs_results[optimal_k],
            'component_scores': component_scores,
            'validation_enabled': True
        }
        
        return results
    
    def _run_gibbs_clustering_range(self, 
                                   sequences: List[str], 
                                   k_range: range, 
                                   **kwargs) -> Dict[int, Any]:
        """Run Gibbs clustering for a range of k values"""
        results = {}
        
        for k in k_range:
            print(f"   Processing k={k}...")
            try:
                result = self._run_single_gibbs_clustering(sequences, k, **kwargs)
                results[k] = result
            except Exception as e:
                print(f"   Warning: Clustering failed for k={k}: {e}")
                # Create dummy result to maintain structure
                results[k] = (pd.DataFrame({'Sequence': sequences, 'Gn': [1] * len(sequences)}), [], {})
        
        return results
    
    def _run_single_gibbs_clustering(self, 
                                    sequences: List[str], 
                                    k: int, 
                                    **kwargs) -> Any:
        """Run single Gibbs clustering"""
        
        # Filter out validation-specific parameters that shouldn't go to clustering classes
        validation_params = {
            'enable_bootstrap', 'n_bootstrap_iterations', 'sample_fraction',
            'n_consensus_iterations', 'validation_enabled'
        }
        clustering_kwargs = {key: value for key, value in kwargs.items() 
                           if key not in validation_params}
        
        if self.use_advanced:
            # Use advanced Gibbs clustering
            clusterer = GibbsClusterAdvanced(
                motif_length=self.motif_length,
                num_clusters=k,
                seed=self.random_seed,
                **clustering_kwargs
            )
            clusterer.fit(sequences)
            return clusterer
        else:
            # Filter out advanced-only parameters for basic clustering
            basic_only_params = {
                'num_seeds', 'iterations', 'temperature_start', 'temperature_steps',
                'lambda_penalty', 'sigma_weight', 'sequence_weighting',
                'background_model', 'use_trash_cluster', 'trash_threshold'
            }
            basic_kwargs = {key: value for key, value in clustering_kwargs.items() 
                          if key not in basic_only_params}
            
            # Use basic Gibbs clustering
            return gibbs_cluster(
                sequences,
                motif_length=self.motif_length,
                num_clusters=k,
                seed=self.random_seed,
                **basic_kwargs
            )
    
    def plot_validation_summary(self, save_path: Optional[str] = None) -> plt.Figure:
        """
        Create validation summary plot
        
        Parameters:
        -----------
        save_path : str, optional
            Path to save the plot
            
        Returns:
        --------
        plt.Figure : Validation summary figure
        """
        if not self.validation_enabled or not self.validation_results:
            raise ValueError("Validation must be enabled and results must be available")
        
        fig = self.validator.plot_validation_results(
            self.validation_results['stability_scores'],
            self.validation_results['silhouette_scores'],
            self.validation_results.get('bootstrap_scores'),
            self.optimal_k,
            save_path
        )
        
        return fig
    
    def plot_consensus_matrix(self, 
                             k: Optional[int] = None, 
                             save_path: Optional[str] = None) -> plt.Figure:
        """
        Plot consensus matrix for specified k
        
        Parameters:
        -----------
        k : int, optional
            Number of clusters (defaults to optimal_k)
        save_path : str, optional
            Path to save the plot
            
        Returns:
        --------
        plt.Figure : Consensus matrix figure
        """
        if not self.validation_enabled or not self.validation_results:
            raise ValueError("Validation must be enabled and results must be available")
        
        if k is None:
            k = self.optimal_k
        
        if k not in self.validation_results['consensus_matrices']:
            raise ValueError(f"Consensus matrix not available for k={k}")
        
        consensus_matrix = self.validation_results['consensus_matrices'][k]
        
        fig = self.validator.plot_consensus_matrix(
            consensus_matrix, 
            k, 
            save_path=save_path
        )
        
        return fig
    
    def generate_validation_report(self, 
                                  sequences: List[str],
                                  save_path: Optional[str] = None) -> str:
        """
        Generate comprehensive validation report
        
        Parameters:
        -----------
        sequences : List[str]
            Input sequences
        save_path : str, optional
            Path to save the report
            
        Returns:
        --------
        str : Validation report
        """
        if not self.validation_enabled or not self.validation_results:
            raise ValueError("Validation must be enabled and results must be available")
        
        report = self.validator.generate_validation_report(
            sequences,
            self.validation_results,
            self.optimal_k,
            save_path
        )
        
        return report
    
    def get_final_clustering_results(self) -> Tuple[pd.DataFrame, List, Dict]:
        """
        Get the final clustering results for the optimal k
        
        Returns:
        --------
        Tuple[pd.DataFrame, List, Dict] : (clusters_df, pwms, logos)
        """
        if self.final_clusters is None:
            raise ValueError("No clustering results available. Run clustering first.")
        
        if self.use_advanced:
            # Advanced clusterer
            if hasattr(self.final_clusters, 'get_clusters'):
                clusters_df = self.final_clusters.get_clusters()
                pwms = self.final_clusters.get_pwms() if hasattr(self.final_clusters, 'get_pwms') else []
                logos = self.final_clusters.get_logos() if hasattr(self.final_clusters, 'get_logos') else {}
                return clusters_df, pwms, logos
            else:
                # Fallback
                return pd.DataFrame(), [], {}
        else:
            # Basic clusterer returns tuple
            if isinstance(self.final_clusters, tuple) and len(self.final_clusters) >= 3:
                return self.final_clusters[:3]
            else:
                return pd.DataFrame(), [], {}
    
    def get_validation_summary(self) -> Dict[str, float]:
        """
        Get summary of validation metrics
        
        Returns:
        --------
        Dict[str, float] : Validation summary metrics
        """
        if not self.validation_enabled or not self.validation_results:
            return {'validation_enabled': False}
        
        summary = {
            'validation_enabled': True,
            'optimal_k': self.optimal_k,
            'stability_score': self.validation_results['component_scores']['stability'],
            'silhouette_score': self.validation_results['component_scores']['silhouette'],
            'composite_score': self.validation_results['component_scores']['composite']
        }
        
        if 'bootstrap' in self.validation_results['component_scores']:
            summary['bootstrap_score'] = self.validation_results['component_scores']['bootstrap']
        
        return summary

def enhanced_gibbs_cluster_with_validation(sequences: List[str],
                                          motif_length: int = 8,
                                          k_range: Union[range, Tuple[int, int]] = (2, 6),
                                          use_advanced: bool = True,
                                          validation_enabled: bool = True,
                                          n_consensus_iterations: int = 50,
                                          n_bootstrap_iterations: int = 20,
                                          random_seed: int = 42,
                                          **clustering_kwargs) -> Tuple[pd.DataFrame, List, Dict, Dict]:
    """
    Convenience function for enhanced Gibbs clustering with validation
    
    Parameters:
    -----------
    sequences : List[str]
        Input peptide sequences
    motif_length : int
        Length of motifs to discover
    k_range : Union[range, Tuple[int, int]]
        Range of cluster numbers to test
    use_advanced : bool
        Whether to use advanced Gibbs clustering
    validation_enabled : bool
        Whether to perform consensus validation
    n_consensus_iterations : int
        Number of consensus clustering iterations
    n_bootstrap_iterations : int
        Number of bootstrap validation iterations  
    random_seed : int
        Random seed for reproducibility
    **clustering_kwargs : dict
        Additional parameters for Gibbs clustering
        
    Returns:
    --------
    Tuple[pd.DataFrame, List, Dict, Dict] : (clusters_df, pwms, logos, validation_results)
    """
    
    # Initialize enhanced clustering
    enhanced_clustering = EnhancedGibbsClustering(
        motif_length=motif_length,
        use_advanced=use_advanced,
        validation_enabled=validation_enabled,
        random_seed=random_seed
    )
    
    # Run clustering with validation
    results = enhanced_clustering.run_clustering_with_validation(
        sequences,
        k_range=k_range,
        n_consensus_iterations=n_consensus_iterations,
        n_bootstrap_iterations=n_bootstrap_iterations,
        **clustering_kwargs
    )
    
    # Get final clustering results
    clusters_df, pwms, logos = enhanced_clustering.get_final_clustering_results()
    
    # Get validation summary
    validation_summary = enhanced_clustering.get_validation_summary()
    
    return clusters_df, pwms, logos, validation_summary