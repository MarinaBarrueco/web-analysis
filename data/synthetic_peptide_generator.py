#!/usr/bin/env python3
"""
Synthetic Peptide Data Generator

Creates realistic peptide library datasets with tunable parameters for testing
the analysis pipeline. Generates peptides following specific patterns with
controllable differential expression patterns and noise characteristics.

Usage:
    python synthetic_peptide_generator.py --help
    python synthetic_peptide_generator.py --output synthetic_data.csv --n_peptides 10000
"""

import numpy as np
import pandas as pd
import argparse
import random
from typing import List, Tuple, Dict
import re

class SyntheticPeptideGenerator:
    """Generate synthetic peptide data with realistic characteristics"""
    
    def __init__(self, seed: int = 42):
        """Initialize generator with random seed"""
        np.random.seed(seed)
        random.seed(seed)
        
        # Amino acid frequencies (approximate natural frequencies)
        self.aa_frequencies = {
            'A': 0.082, 'R': 0.057, 'N': 0.043, 'D': 0.054, 'C': 0.013,
            'Q': 0.039, 'E': 0.063, 'G': 0.071, 'H': 0.022, 'I': 0.059,
            'L': 0.096, 'K': 0.058, 'M': 0.024, 'F': 0.039, 'P': 0.047,
            'S': 0.068, 'T': 0.053, 'W': 0.013, 'Y': 0.032, 'V': 0.068
        }
        
        self.amino_acids = list(self.aa_frequencies.keys())
        self.aa_weights = list(self.aa_frequencies.values())
    
    def generate_random_peptide(self, length: int, pattern: str = None) -> str:
        """Generate a random peptide sequence of given length"""
        if pattern:
            return self._generate_pattern_peptide(pattern)
        else:
            return ''.join(np.random.choice(self.amino_acids, size=length, p=self.aa_weights))
    
    def _generate_pattern_peptide(self, pattern: str) -> str:
        """Generate peptide following a specific pattern (regex-like)"""
        # Convert pattern to actual sequence
        result = ""
        i = 0
        while i < len(pattern):
            char = pattern[i]
            if char == '.':
                # Any amino acid
                result += np.random.choice(self.amino_acids, p=self.aa_weights)
            elif char == '{' and i+2 < len(pattern) and pattern[i+2] == '}':
                # Repeat pattern like {7}
                n_repeats = int(pattern[i+1])
                for _ in range(n_repeats):
                    result += np.random.choice(self.amino_acids, p=self.aa_weights)
                i += 2  # Skip the {n}
            else:
                # Fixed amino acid
                result += char
            i += 1
        return result
    
    def generate_motif_variants(self, core_motif: str, n_variants: int, 
                               mutation_rate: float = 0.1) -> List[str]:
        """Generate variants of a core motif with controlled mutations"""
        variants = [core_motif]  # Include original
        
        for _ in range(n_variants - 1):
            variant = list(core_motif)
            # Apply mutations
            for j in range(len(variant)):
                if np.random.random() < mutation_rate:
                    variant[j] = np.random.choice(self.amino_acids, p=self.aa_weights)
            variants.append(''.join(variant))
        
        return variants
    
    def generate_count_data(self, peptides: List[str], 
                           control_samples: int = 3,
                           treatment_samples: int = 3,
                           base_expression: float = 1000,
                           differential_peptides: Dict[str, float] = None,
                           noise_level: float = 0.2,
                           dropout_rate: float = 0.1) -> pd.DataFrame:
        """
        Generate realistic count data with differential expression patterns
        
        Parameters:
        -----------
        peptides : List[str]
            List of peptide sequences
        control_samples : int
            Number of control samples
        treatment_samples : int  
            Number of treatment samples
        base_expression : float
            Base expression level (mean counts)
        differential_peptides : Dict[str, float]
            Peptides with fold changes {peptide: fold_change}
        noise_level : float
            Biological noise (coefficient of variation)
        dropout_rate : float
            Probability of zero counts (dropout events)
        """
        
        # Initialize data structure
        sample_names = ([f"Control_{i+1}" for i in range(control_samples)] + 
                       [f"Treatment_{i+1}" for i in range(treatment_samples)])
        
        data = {'peptide': peptides}
        
        if differential_peptides is None:
            differential_peptides = {}
        
        # Generate counts for each sample
        for i, sample_name in enumerate(sample_names):
            is_treatment = i >= control_samples
            counts = []
            
            for peptide in peptides:
                # Base expression with log-normal distribution
                base_count = np.random.lognormal(
                    mean=np.log(base_expression), 
                    sigma=noise_level
                )
                
                # Apply differential expression
                if peptide in differential_peptides and is_treatment:
                    fold_change = differential_peptides[peptide]
                    base_count *= fold_change
                
                # Add technical noise
                count = np.random.poisson(base_count)
                
                # Apply dropout
                if np.random.random() < dropout_rate:
                    count = 0
                
                counts.append(int(count))
            
            data[sample_name] = counts
        
        return pd.DataFrame(data)
    
    def generate_clustered_dataset(self, 
                                  pattern: str = "A.C.{7}C",
                                  n_total_peptides: int = 10000,
                                  n_differential: int = 500,
                                  n_clusters: int = 4,
                                  control_samples: int = 3,
                                  treatment_samples: int = 3,
                                  fold_change_range: Tuple[float, float] = (2.0, 10.0),
                                  base_expression: float = 100,
                                  noise_level: float = 0.3) -> pd.DataFrame:
        """
        Generate a complete synthetic dataset with clustered motifs
        
        Parameters:
        -----------
        pattern : str
            Peptide library pattern (regex-like)
        n_total_peptides : int
            Total number of peptides to generate
        n_differential : int
            Number of differentially expressed peptides
        n_clusters : int
            Number of distinct motif clusters in differential peptides
        fold_change_range : Tuple[float, float]
            Range of fold changes for differential peptides
        """
        
        print(f"🧬 Generating synthetic dataset with {n_total_peptides} peptides...")
        
        # Generate background peptides (non-differential)
        background_peptides = []
        for _ in range(n_total_peptides - n_differential):
            peptide = self._generate_pattern_peptide(pattern)
            background_peptides.append(peptide)
        
        # Generate clustered differential peptides
        differential_peptides = {}
        cluster_size = n_differential // n_clusters
        
        for cluster_id in range(n_clusters):
            print(f"  📊 Generating cluster {cluster_id + 1}/{n_clusters}...")
            
            # Create a core motif for this cluster
            core_motif = self._generate_pattern_peptide(pattern)
            
            # Generate variants of this motif
            cluster_peptides = self.generate_motif_variants(
                core_motif, 
                cluster_size,
                mutation_rate=0.15  # 15% mutation rate for realistic clusters
            )
            
            # Assign fold changes to cluster peptides
            for peptide in cluster_peptides:
                fold_change = np.random.uniform(*fold_change_range)
                differential_peptides[peptide] = fold_change
        
        # Combine all peptides
        all_peptides = background_peptides + list(differential_peptides.keys())
        
        # Remove duplicates while preserving differential info
        unique_peptides = []
        final_differential = {}
        seen = set()
        
        for peptide in all_peptides:
            if peptide not in seen:
                unique_peptides.append(peptide)
                seen.add(peptide)
                if peptide in differential_peptides:
                    final_differential[peptide] = differential_peptides[peptide]
        
        print(f"  ✅ Generated {len(unique_peptides)} unique peptides")
        print(f"  📈 {len(final_differential)} differential peptides across {n_clusters} clusters")
        
        # Generate count data
        print("  🔢 Generating count data...")
        dataset = self.generate_count_data(
            unique_peptides,
            control_samples=control_samples,
            treatment_samples=treatment_samples,
            base_expression=base_expression,
            differential_peptides=final_differential,
            noise_level=noise_level
        )
        
        print(f"  ✅ Dataset generated: {dataset.shape[0]} peptides × {dataset.shape[1]-1} samples")
        
        return dataset

def main():
    parser = argparse.ArgumentParser(description='Generate synthetic peptide data for testing')
    
    parser.add_argument('--output', '-o', 
                       default='synthetic_peptide_data.csv',
                       help='Output CSV file name')
    
    parser.add_argument('--n_peptides', '-n', 
                       type=int, default=10000,
                       help='Total number of peptides to generate')
    
    parser.add_argument('--n_differential', '-d',
                       type=int, default=500,
                       help='Number of differentially expressed peptides')
    
    parser.add_argument('--n_clusters', '-c',
                       type=int, default=4,
                       help='Number of motif clusters')
    
    parser.add_argument('--pattern', '-p',
                       default='A.C.{7}C',
                       help='Peptide library pattern (regex-like)')
    
    parser.add_argument('--control_samples',
                       type=int, default=3,
                       help='Number of control samples')
    
    parser.add_argument('--treatment_samples',
                       type=int, default=3,
                       help='Number of treatment samples')
    
    parser.add_argument('--fold_change_min',
                       type=float, default=2.0,
                       help='Minimum fold change for differential peptides')
    
    parser.add_argument('--fold_change_max',
                       type=float, default=10.0,
                       help='Maximum fold change for differential peptides')
    
    parser.add_argument('--base_expression',
                       type=float, default=100,
                       help='Base expression level (mean counts)')
    
    parser.add_argument('--noise_level',
                       type=float, default=0.3,
                       help='Biological noise level (0.1-0.5)')
    
    parser.add_argument('--seed',
                       type=int, default=42,
                       help='Random seed for reproducibility')
    
    args = parser.parse_args()
    
    # Create generator
    generator = SyntheticPeptideGenerator(seed=args.seed)
    
    # Generate dataset
    dataset = generator.generate_clustered_dataset(
        pattern=args.pattern,
        n_total_peptides=args.n_peptides,
        n_differential=args.n_differential,
        n_clusters=args.n_clusters,
        control_samples=args.control_samples,
        treatment_samples=args.treatment_samples,
        fold_change_range=(args.fold_change_min, args.fold_change_max),
        base_expression=args.base_expression,
        noise_level=args.noise_level
    )
    
    # Save dataset
    dataset.to_csv(args.output, index=False)
    print(f"💾 Dataset saved to: {args.output}")
    
    # Print summary statistics
    print("\n📊 Dataset Summary:")
    print(f"  • Total peptides: {len(dataset)}")
    print(f"  • Samples: {len([col for col in dataset.columns if col != 'peptide'])}")
    print(f"  • Control samples: {args.control_samples}")
    print(f"  • Treatment samples: {args.treatment_samples}")
    print(f"  • Expected differential peptides: {args.n_differential}")
    print(f"  • Expected clusters: {args.n_clusters}")
    print(f"  • Pattern: {args.pattern}")
    
    # Show data preview
    print(f"\n🔍 Data Preview:")
    print(dataset.head())

if __name__ == "__main__":
    main()