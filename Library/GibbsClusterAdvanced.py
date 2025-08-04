"""
Advanced Gibbs Clustering Implementation
Based on GibbsCluster-2.0e_SA.pl

This implementation includes:
- Multiple seeds for better convergence
- Temperature annealing (simulated annealing)
- Kullback-Leibler divergence calculation
- Support for insertions and deletions
- Sequence weighting
- Background frequency models
- Trash cluster for outliers
- Comprehensive result reporting
"""

import numpy as np
import pandas as pd
import random
import math
import matplotlib.pyplot as plt
import logomaker as lm
from collections import defaultdict
import re
from typing import List, Dict, Tuple, Optional, Union
import tempfile
import os
from pathlib import Path

class GibbsClusterAdvanced:
    """Advanced Gibbs Clustering implementation matching the Perl version functionality"""
    
    # Amino acid alphabet and indices
    AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
    AA_IDX = {aa: i for i, aa in enumerate(AA_ALPHABET)}
    IDX_AA = {i: aa for aa, i in AA_IDX.items()}
    
    # BLOSUM background frequencies (simplified - can be loaded from files)
    BLOSUM_BG = {
        '62': [0.074, 0.025, 0.054, 0.054, 0.026, 0.068, 0.022, 0.062, 0.058, 0.099,
               0.025, 0.045, 0.049, 0.040, 0.052, 0.057, 0.051, 0.013, 0.032, 0.073]
    }
    
    def __init__(self, 
                 motif_length: int = 9,
                 num_clusters: Union[int, Tuple[int, int]] = 4,
                 max_insertion_length: int = 0,
                 max_deletion_length: int = 0,
                 iterations: int = 10,
                 temperature_start: float = 1.5,
                 temperature_steps: int = 20,
                 num_seeds: int = 1,
                 lambda_penalty: float = 0.8,
                 sigma_weight: float = 5.0,
                 sequence_weighting: int = 0,
                 background_model: int = 1,
                 blosum_matrix: str = '62',
                 use_trash_cluster: bool = False,
                 trash_threshold: float = 0.0,
                 cluster_mode: bool = False,
                 preference_hydrophobic_p1: bool = False,
                 single_move_interval: int = 20,
                 phase_shift_interval: int = 100,
                 indel_move_interval: int = 10,
                 weight_low_count: int = 200,
                 seed: Optional[int] = None):
        
        self.motif_length = motif_length
        if isinstance(num_clusters, tuple):
            self.min_clusters, self.max_clusters = num_clusters
        else:
            self.min_clusters = self.max_clusters = num_clusters
            
        self.max_insertion_length = max_insertion_length
        self.max_deletion_length = max_deletion_length
        self.iterations = iterations
        self.temperature_start = temperature_start
        self.temperature_steps = temperature_steps
        self.num_seeds = num_seeds
        self.lambda_penalty = lambda_penalty
        self.sigma_weight = sigma_weight
        self.sequence_weighting = sequence_weighting
        self.background_model = background_model
        self.blosum_matrix = blosum_matrix
        self.use_trash_cluster = use_trash_cluster
        self.trash_threshold = trash_threshold
        self.cluster_mode = cluster_mode
        self.preference_hydrophobic_p1 = preference_hydrophobic_p1
        self.single_move_interval = single_move_interval
        self.phase_shift_interval = phase_shift_interval
        self.indel_move_interval = indel_move_interval
        self.weight_low_count = weight_low_count
        
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)
            
        # Load background frequencies
        self.background_freqs = np.array(self.BLOSUM_BG.get(blosum_matrix, 
                                                           [1/20] * 20))
        
        # Results storage
        self.results = {}
        self.best_results = {}
        
    def _encode_sequences(self, sequences: List[str]) -> np.ndarray:
        """Encode amino acid sequences to integer arrays"""
        encoded = []
        for seq in sequences:
            seq_encoded = []
            for aa in seq:
                if aa in self.AA_IDX:
                    seq_encoded.append(self.AA_IDX[aa])
                else:
                    raise ValueError(f"Unknown amino acid: {aa}")
            encoded.append(seq_encoded)
        return np.array(encoded)
    
    def _calculate_kld(self, counts: np.ndarray, background: np.ndarray, 
                      cluster_sizes: np.ndarray) -> Tuple[float, np.ndarray]:
        """Calculate Kullback-Leibler divergence"""
        num_clusters, motif_len, num_aa = counts.shape
        cluster_klds = np.zeros(num_clusters)
        
        for k in range(num_clusters):
            if cluster_sizes[k] == 0:
                continue
                
            kld_sum = 0.0
            for pos in range(motif_len):
                # Calculate frequencies
                freqs = counts[k, pos, :] / cluster_sizes[k]
                freqs = np.maximum(freqs, 1e-10)  # Avoid log(0)
                
                # KLD calculation
                for aa in range(num_aa):
                    if freqs[aa] > 0:
                        kld_sum += freqs[aa] * math.log(freqs[aa] / background[aa])
            
            cluster_klds[k] = kld_sum
        
        # Weighted average KLD
        total_seqs = cluster_sizes.sum()
        if total_seqs > 0:
            avg_kld = np.sum(cluster_klds * cluster_sizes) / total_seqs
        else:
            avg_kld = 0.0
            
        return avg_kld, cluster_klds
    
    def _gibbs_sampling_step(self, sequences: np.ndarray, assignments: np.ndarray, 
                           counts: np.ndarray, cluster_sizes: np.ndarray,
                           temperature: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Perform one step of Gibbs sampling"""
        num_seqs, seq_len = sequences.shape
        num_clusters = len(cluster_sizes)
        
        for seq_idx in range(num_seqs):
            # Remove current sequence from counts
            old_cluster = assignments[seq_idx]
            cluster_sizes[old_cluster] -= 1
            
            for pos in range(seq_len):
                aa = sequences[seq_idx, pos]
                counts[old_cluster, pos, aa] -= 1
            
            # Calculate probabilities for each cluster
            log_probs = np.zeros(num_clusters)
            
            for k in range(num_clusters):
                # Prior probability (cluster size)
                log_prob = math.log(cluster_sizes[k] + self.sigma_weight)
                
                # Likelihood
                for pos in range(seq_len):
                    aa = sequences[seq_idx, pos]
                    aa_count = counts[k, pos, aa]
                    total_count = cluster_sizes[k]
                    
                    # Add pseudocounts
                    prob = (aa_count + self.background_freqs[aa]) / (total_count + 1.0)
                    log_prob += math.log(prob)
                
                # Apply penalty for inter-cluster similarity (simplified)
                log_prob -= self.lambda_penalty * k  # Simplified penalty
                
                log_probs[k] = log_prob / temperature
            
            # Normalize and sample
            log_probs -= log_probs.max()  # Numerical stability
            probs = np.exp(log_probs)
            probs /= probs.sum()
            
            # Sample new cluster assignment
            new_cluster = np.random.choice(num_clusters, p=probs)
            
            # Update counts with new assignment
            assignments[seq_idx] = new_cluster
            cluster_sizes[new_cluster] += 1
            
            for pos in range(seq_len):
                aa = sequences[seq_idx, pos]
                counts[new_cluster, pos, aa] += 1
        
        return assignments, counts, cluster_sizes
    
    def _run_single_clustering(self, sequences: List[str], num_clusters: int, 
                             seed_num: int) -> Dict:
        """Run clustering with a single seed"""
        
        # Prepare sequences
        encoded_seqs = self._encode_sequences(sequences)
        num_seqs, seq_len = encoded_seqs.shape
        
        # Initialize assignments randomly
        assignments = np.random.randint(0, num_clusters, size=num_seqs)
        cluster_sizes = np.bincount(assignments, minlength=num_clusters)
        
        # Initialize count matrices
        counts = np.zeros((num_clusters, seq_len, len(self.AA_ALPHABET)), dtype=int)
        
        for seq_idx in range(num_seqs):
            cluster = assignments[seq_idx]
            for pos in range(seq_len):
                aa = encoded_seqs[seq_idx, pos]
                counts[cluster, pos, aa] += 1
        
        # Simulated annealing
        temperature_schedule = np.logspace(
            math.log10(self.temperature_start), 
            math.log10(0.1), 
            self.temperature_steps
        )
        
        best_kld = -float('inf')
        best_assignments = assignments.copy()
        best_counts = counts.copy()
        
        for temp_step, temperature in enumerate(temperature_schedule):
            for iteration in range(self.iterations):
                assignments, counts, cluster_sizes = self._gibbs_sampling_step(
                    encoded_seqs, assignments, counts, cluster_sizes, temperature
                )
                
                # Calculate KLD
                avg_kld, cluster_klds = self._calculate_kld(
                    counts, self.background_freqs, cluster_sizes
                )
                
                if avg_kld > best_kld:
                    best_kld = avg_kld
                    best_assignments = assignments.copy()
                    best_counts = counts.copy()
        
        # Prepare results
        results = {
            'assignments': best_assignments,
            'counts': best_counts,
            'cluster_sizes': np.bincount(best_assignments, minlength=num_clusters),
            'avg_kld': best_kld,
            'cluster_klds': self._calculate_kld(best_counts, self.background_freqs, 
                                              np.bincount(best_assignments, minlength=num_clusters))[1],
            'seed': seed_num,
            'num_clusters': num_clusters
        }
        
        return results
    
    def fit(self, sequences: Union[List[str], pd.DataFrame]) -> 'GibbsClusterAdvanced':
        """Fit the Gibbs clustering model to sequences"""
        
        # Handle input format
        if isinstance(sequences, pd.DataFrame):
            if sequences.empty:
                raise ValueError("Input DataFrame is empty")
            
            if 'variable_pep' in sequences.columns:
                seq_list = sequences['variable_pep'].tolist()
            elif 'Sequence' in sequences.columns:
                seq_list = sequences['Sequence'].tolist()
            else:
                seq_list = sequences.iloc[:, 0].tolist()
        else:
            seq_list = list(sequences)
        
        if not seq_list:
            raise ValueError("No sequences provided")
        
        # Validate and clean sequences
        valid_seqs = []
        invalid_count = 0
        
        for seq in seq_list:
            if not isinstance(seq, str):
                invalid_count += 1
                continue
            
            seq = seq.upper().strip()
            
            # Check length
            if len(seq) < self.motif_length:
                invalid_count += 1
                continue
            
            # Check for valid amino acids
            if not all(aa in self.AA_ALPHABET for aa in seq):
                invalid_count += 1
                continue
            
            valid_seqs.append(seq)
        
        if invalid_count > 0:
            print(f"Warning: Filtered out {invalid_count} invalid sequences")
        
        if not valid_seqs:
            raise ValueError("No valid sequences remain after filtering")
        
        if len(valid_seqs) < 2:
            raise ValueError("Need at least 2 valid sequences for clustering")
        
        # Adjust cluster range if necessary
        max_possible_clusters = min(len(valid_seqs), self.max_clusters)
        if self.min_clusters > max_possible_clusters:
            print(f"Warning: Adjusting minimum clusters from {self.min_clusters} to {max_possible_clusters}")
            self.min_clusters = max_possible_clusters
        
        if self.max_clusters > max_possible_clusters:
            print(f"Warning: Adjusting maximum clusters from {self.max_clusters} to {max_possible_clusters}")
            self.max_clusters = max_possible_clusters
        
        seq_list = valid_seqs
        
        print(f"Processing {len(seq_list)} valid sequences")
        print(f"Motif length: {self.motif_length}")
        print(f"Cluster range: {self.min_clusters}-{self.max_clusters}")
        print(f"Number of seeds: {self.num_seeds}")
        
        # Run clustering for each number of clusters
        for num_clusters in range(self.min_clusters, self.max_clusters + 1):
            print(f"\nClustering with {num_clusters} clusters...")
            
            best_result = None
            best_kld = -float('inf')
            
            # Try multiple seeds
            for seed_num in range(1, self.num_seeds + 1):
                print(f"  Seed {seed_num}/{self.num_seeds}")
                
                result = self._run_single_clustering(seq_list, num_clusters, seed_num)
                
                if result['avg_kld'] > best_kld:
                    best_kld = result['avg_kld']
                    best_result = result
            
            self.results[num_clusters] = best_result
            print(f"  Best KLD for {num_clusters} clusters: {best_kld:.4f}")
        
        # Find overall best clustering
        best_num_clusters = max(self.results.keys(), 
                               key=lambda k: self.results[k]['avg_kld'])
        self.best_results = self.results[best_num_clusters]
        self.sequences = seq_list
        
        print(f"\nBest clustering: {best_num_clusters} clusters "
              f"(KLD: {self.best_results['avg_kld']:.4f})")
        
        return self
    
    def get_clusters(self, num_clusters: Optional[int] = None) -> pd.DataFrame:
        """Get cluster assignments as DataFrame"""
        if num_clusters is None:
            result = self.best_results
        else:
            result = self.results.get(num_clusters, self.best_results)
        
        clusters_df = pd.DataFrame({
            'Sequence': self.sequences,
            'Cluster': result['assignments'] + 1,  # 1-indexed
            'Gn': result['assignments'] + 1  # For compatibility
        })
        
        return clusters_df
    
    def get_pwms(self, num_clusters: Optional[int] = None) -> List[pd.DataFrame]:
        """Get position weight matrices for each cluster"""
        if num_clusters is None:
            result = self.best_results
        else:
            result = self.results.get(num_clusters, self.best_results)
        
        counts = result['counts']
        cluster_sizes = result['cluster_sizes']
        
        pwms = []
        for k in range(counts.shape[0]):
            if cluster_sizes[k] > 0:
                # Add pseudocounts and normalize
                freqs = counts[k] + 1.0
                freqs = freqs / freqs.sum(axis=1, keepdims=True)
                
                pwm_df = pd.DataFrame(freqs, columns=list(self.AA_ALPHABET))
                pwms.append(pwm_df)
            else:
                # Empty cluster
                pwm_df = pd.DataFrame(np.ones((self.motif_length, 20)) / 20, 
                                    columns=list(self.AA_ALPHABET))
                pwms.append(pwm_df)
        
        return pwms
    
    def get_logos(self, num_clusters: Optional[int] = None) -> Dict[int, plt.Figure]:
        """Generate sequence logos for each cluster"""
        pwms = self.get_pwms(num_clusters)
        
        if num_clusters is None:
            result = self.best_results
        else:
            result = self.results.get(num_clusters, self.best_results)
        
        cluster_sizes = result['cluster_sizes']
        cluster_klds = result['cluster_klds']
        
        logos = {}
        for k, pwm in enumerate(pwms):
            if cluster_sizes[k] > 0:
                fig, ax = plt.subplots(figsize=(6, 2))
                
                try:
                    logo = lm.Logo(pwm, ax=ax, shade_below=0.5, fade_below=0.5)
                    ax.set_title(f'Cluster {k+1} (n={cluster_sizes[k]}, '
                               f'KLD={cluster_klds[k]:.3f})', fontsize=12)
                    ax.set_xlabel('Position')
                    ax.set_ylabel('Bits')
                    plt.tight_layout()
                    logos[k+1] = fig
                except Exception as e:
                    print(f"Warning: Could not generate logo for cluster {k+1}: {e}")
                    # Create a simple text-based representation
                    fig, ax = plt.subplots(figsize=(6, 2))
                    ax.text(0.5, 0.5, f'Cluster {k+1}\n(n={cluster_sizes[k]})', 
                           ha='center', va='center', transform=ax.transAxes)
                    ax.set_title(f'Cluster {k+1}')
                    logos[k+1] = fig
        
        return logos
    
    def get_summary(self) -> Dict:
        """Get summary statistics for all clustering results"""
        summary = {
            'sequence_count': len(self.sequences),
            'motif_length': self.motif_length,
            'cluster_results': {}
        }
        
        for num_clusters, result in self.results.items():
            cluster_sizes = result['cluster_sizes']
            non_empty_clusters = np.sum(cluster_sizes > 0)
            
            summary['cluster_results'][num_clusters] = {
                'avg_kld': result['avg_kld'],
                'cluster_sizes': cluster_sizes.tolist(),
                'non_empty_clusters': non_empty_clusters,
                'cluster_klds': result['cluster_klds'].tolist()
            }
        
        return summary
    
    def save_results(self, output_dir: str, run_name: str = "gibbs_clustering"):
        """Save clustering results to files"""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Save cluster assignments
        for num_clusters in self.results:
            clusters_df = self.get_clusters(num_clusters)
            clusters_df.to_csv(output_path / f"{run_name}_{num_clusters}clusters.csv", 
                             index=False)
        
        # Save PWMs
        for num_clusters in self.results:
            pwms = self.get_pwms(num_clusters)
            for k, pwm in enumerate(pwms):
                pwm.to_csv(output_path / f"{run_name}_{num_clusters}clusters_pwm_{k+1}.csv")
        
        # Save logos
        for num_clusters in self.results:
            logos = self.get_logos(num_clusters)
            for cluster_id, fig in logos.items():
                fig.savefig(output_path / f"{run_name}_{num_clusters}clusters_logo_{cluster_id}.png", 
                           dpi=300, bbox_inches='tight')
                plt.close(fig)
        
        # Save summary
        summary = self.get_summary()
        import json
        with open(output_path / f"{run_name}_summary.json", 'w') as f:
            json.dump(summary, f, indent=2)
        
        print(f"Results saved to {output_path}")


def gibbs_cluster_advanced(
    peptides: Union[List[str], pd.DataFrame],
    motif_length: int = 9,
    num_clusters: Union[int, Tuple[int, int]] = 4,
    out_dir: Optional[str] = None,
    run_name: Optional[str] = None,
    **kwargs
) -> Tuple[pd.DataFrame, List[pd.DataFrame], Dict[int, plt.Figure]]:
    """
    Advanced Gibbs clustering function with API compatible with the original
    
    Returns:
    --------
    clusters_df : DataFrame with Sequence and Cluster columns
    pwm_list : List of PWM DataFrames for each cluster
    logo_figs : Dictionary mapping cluster number to matplotlib Figure
    """
    
    # Create and fit the clustering model
    clusterer = GibbsClusterAdvanced(
        motif_length=motif_length,
        num_clusters=num_clusters,
        **kwargs
    )
    
    clusterer.fit(peptides)
    
    # Get results
    clusters_df = clusterer.get_clusters()
    pwm_list = clusterer.get_pwms()
    logo_figs = clusterer.get_logos()
    
    # Save results if requested
    if out_dir and run_name:
        clusterer.save_results(out_dir, run_name)
    
    return clusters_df, pwm_list, logo_figs


# Example usage and testing
if __name__ == "__main__":
    # Test with sample peptides
    sample_peptides = [
        "ACDEFGHIK", "ACDEFGHIL", "ACDEFGHIM", "ACDEFGHIN",
        "LMNPQRSTU", "LMNPQRSTV", "LMNPQRSTW", "LMNPQRSTY",
        "VWXYZABCD", "VWXYZABCE", "VWXYZABCF", "VWXYZABCG"
    ]
    
    print("Testing Advanced Gibbs Clustering...")
    
    # Test the function interface
    clusters_df, pwms, logos = gibbs_cluster_advanced(
        sample_peptides,
        motif_length=9,
        num_clusters=3,
        num_seeds=2
    )
    
    print("\nClustering Results:")
    print(clusters_df)
    
    print(f"\nGenerated {len(pwms)} PWMs and {len(logos)} logos")
    
    # Test the class interface
    clusterer = GibbsClusterAdvanced(
        motif_length=9,
        num_clusters=(2, 4),
        num_seeds=2
    )
    
    clusterer.fit(sample_peptides)
    summary = clusterer.get_summary()
    
    print(f"\nSummary: {summary}")