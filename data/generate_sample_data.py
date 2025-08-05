#!/usr/bin/env python3
"""
Generate sample synthetic datasets for testing the pipeline
"""

from synthetic_peptide_generator import SyntheticPeptideGenerator
import os

def generate_test_datasets():
    """Generate several test datasets with different characteristics"""
    
    generator = SyntheticPeptideGenerator(seed=42)
    
    # Dataset 1: Small test dataset
    print("🧪 Generating small test dataset...")
    small_dataset = generator.generate_clustered_dataset(
        pattern="A.C.{7}C",
        n_total_peptides=1000,
        n_differential=100,
        n_clusters=3,
        control_samples=3,
        treatment_samples=3,
        fold_change_range=(2.0, 8.0),
        base_expression=200,
        noise_level=0.2
    )
    small_dataset.to_csv("test_small.csv", index=False)
    print("💾 Saved: test_small.csv\n")
    
    # Dataset 2: Medium dataset with high differential signal
    print("🧬 Generating medium dataset with strong signal...")
    medium_dataset = generator.generate_clustered_dataset(
        pattern="A.C.{7}C",
        n_total_peptides=5000,
        n_differential=300,
        n_clusters=4,
        control_samples=3,
        treatment_samples=3,
        fold_change_range=(3.0, 15.0),
        base_expression=150,
        noise_level=0.25
    )
    medium_dataset.to_csv("test_medium.csv", index=False)
    print("💾 Saved: test_medium.csv\n")
    
    # Dataset 3: Large dataset with subtle effects
    print("📊 Generating large dataset with subtle effects...")
    large_dataset = generator.generate_clustered_dataset(
        pattern="A.C.{7}C",
        n_total_peptides=10000,
        n_differential=500,
        n_clusters=5,
        control_samples=4,
        treatment_samples=4,
        fold_change_range=(1.5, 5.0),  # More subtle changes
        base_expression=100,
        noise_level=0.3
    )
    large_dataset.to_csv("test_large.csv", index=False)
    print("💾 Saved: test_large.csv\n")
    
    # Dataset 4: Different peptide pattern
    print("🔬 Generating dataset with different pattern...")
    pattern_dataset = generator.generate_clustered_dataset(
        pattern=".{2}C.{5}C.{2}",  # Different pattern: XX-C-XXXXX-C-XX
        n_total_peptides=3000,
        n_differential=200,
        n_clusters=3,
        control_samples=3,
        treatment_samples=3,
        fold_change_range=(2.5, 10.0),
        base_expression=180,
        noise_level=0.25
    )
    pattern_dataset.to_csv("test_pattern.csv", index=False)
    print("💾 Saved: test_pattern.csv\n")
    
    print("✅ All test datasets generated!")
    print("\n📋 Summary:")
    print("  • test_small.csv: 1K peptides, 100 differential, 3 clusters")
    print("  • test_medium.csv: 5K peptides, 300 differential, 4 clusters")  
    print("  • test_large.csv: 10K peptides, 500 differential, 5 clusters")
    print("  • test_pattern.csv: 3K peptides, different pattern")
    print("\n🚀 Try these in the dashboard by uploading the CSV files!")

if __name__ == "__main__":
    # Change to data directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    generate_test_datasets()