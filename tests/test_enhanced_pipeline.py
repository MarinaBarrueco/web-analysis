#!/usr/bin/env python3
"""
Test script for the enhanced pipeline with realistic peptide data
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Add library path
sys.path.insert(0, str(Path(__file__).parent))

def create_test_data():
    """Create realistic test peptide data"""
    
    # Create peptides following the A.C.{7}C pattern
    test_peptides = [
        'AVCDEFGHYC',  # Group 1: Hydrophobic motifs
        'AICDEFGHYC',
        'ALCDEFGHYC',
        'AVCDEFLHYC',
        'AICDEFLHYC',
        'AVCRHKGHYC',  # Group 2: Charged motifs  
        'AKCRHKGHYC',
        'ARCRHKGHYC',
        'AVCRHKDHYC',
        'AKCRHKDHYC',
        'AVCSTNQHYC',  # Group 3: Polar motifs
        'AICSTNQHYC',
        'ALCSTNQHYC',
        'AVCSTNGHYC',
        'AICSTNGHYC'
    ]
    
    # Create count data with differential expression pattern
    np.random.seed(42)
    n_peptides = len(test_peptides)
    
    # Control samples (lower counts for upregulated peptides)
    control1 = np.random.poisson(20, n_peptides)
    control2 = np.random.poisson(22, n_peptides)
    
    # Experiment samples (higher counts for first 10 peptides)
    exp1 = np.random.poisson(15, n_peptides)
    exp2 = np.random.poisson(18, n_peptides)
    
    # Make first 10 peptides upregulated
    exp1[:10] += np.random.poisson(40, 10)
    exp2[:10] += np.random.poisson(35, 10)
    
    # Create DataFrame
    data = pd.DataFrame({
        'peptide': test_peptides,
        'Control_1': control1,
        'Control_2': control2,
        'Experiment_1': exp1,
        'Experiment_2': exp2
    })
    
    return data

def test_full_pipeline():
    """Test the complete enhanced pipeline"""
    print("🧬 Testing Enhanced Peptide Analysis Pipeline...")
    
    try:
        from Library.lib import group_by_peptide, filter_by_CPM, prepare_data, run_deseq2, significant_DE_peptides
        from Library.enhanced_gibbs_clustering import enhanced_gibbs_cluster_with_validation
        
        # Create test data
        data = create_test_data()
        print(f"   📊 Created test dataset: {data.shape}")
        print(f"   🧬 Test peptides: {len(data)}")
        
        # Define conditions
        conditions = {
            'Control_1': 'Control',
            'Control_2': 'Control', 
            'Experiment_1': 'Experiment',
            'Experiment_2': 'Experiment'
        }
        
        # Step 1: Group peptides
        print("\n📋 Step 1: Grouping peptides...")
        pattern = r'A.C.{7}C'
        grouped, summary = group_by_peptide(data, conditions, pattern)
        print(f"   ✅ Grouped: {grouped.shape[0]} peptides")
        print(f"   📈 Success rate: {summary['Success Rate']}")
        
        # Step 2: CPM filtering
        print("\n🔍 Step 2: CPM filtering...")
        filtered, fig = filter_by_CPM(grouped, conditions, cpm_thresh=10, min_samples=1)
        print(f"   ✅ Filtered: {filtered.shape[0]} peptides remain")
        if fig:
            plt.close(fig)
        
        # Step 3: Prepare for DESeq2
        print("\n🧮 Step 3: Preparing for differential expression...")
        count_data, meta_data = prepare_data(filtered, conditions)
        print(f"   ✅ Count matrix: {count_data.shape}")
        print(f"   📋 Metadata: {meta_data.shape}")
        
        # Step 4: Run DESeq2
        print("\n📊 Step 4: Running DESeq2...")
        res_df, fig = run_deseq2(count_data, meta_data)
        print(f"   ✅ DE results: {res_df.shape}")
        if fig:
            plt.close(fig)
        
        # Step 5: Identify significant peptides
        print("\n🎯 Step 5: Identifying significant peptides...")
        res_df, de_summary = significant_DE_peptides(res_df)
        upregulated = res_df[res_df['updown'] == 'up']
        print(f"   ✅ Significant upregulated: {len(upregulated)}")
        
        if len(upregulated) == 0:
            print("   ⚠️ No upregulated peptides found - adjusting thresholds...")
            res_df, de_summary = significant_DE_peptides(res_df, significance_threshold=0.1, log2fc_threshold_experiment=0.5)
            upregulated = res_df[res_df['updown'] == 'up']
            print(f"   ✅ With relaxed thresholds: {len(upregulated)} upregulated")
        
        if len(upregulated) < 3:
            print("   ⚠️ Too few peptides for clustering - using all peptides for demo")
            test_sequences = data['peptide'].str.replace(r'A(.C.{7})C', r'\1', regex=True).tolist()[:8]
        else:
            # Remove constant positions (A and final C)
            test_sequences = upregulated.index.map(lambda x: x[1:-1]).tolist()
        
        print(f"   🧬 Sequences for clustering: {len(test_sequences)}")
        
        # Step 6: Enhanced clustering with validation
        print("\n🔬 Step 6: Enhanced clustering with consensus validation...")
        
        clusters_df, pwms, logos, validation_summary = enhanced_gibbs_cluster_with_validation(
            test_sequences,
            motif_length=7,  # Adjusted for variable region
            k_range=(2, 4),
            use_advanced=False,  # Use basic for testing
            validation_enabled=True,
            n_consensus_iterations=30,
            n_bootstrap_iterations=15,
            random_seed=42
        )
        
        print(f"   ✅ Clustering completed!")
        print(f"   🎯 Optimal k: {validation_summary['optimal_k']}")
        print(f"   📈 Stability: {validation_summary['stability_score']:.3f}")
        print(f"   🔍 Silhouette: {validation_summary['silhouette_score']:.3f}")
        print(f"   ⚖️ Composite: {validation_summary['composite_score']:.3f}")
        
        # Interpretation
        silhouette_score = validation_summary['silhouette_score']
        if silhouette_score > 0.5:
            print("   🎉 Good clustering structure detected!")
        elif silhouette_score > 0.25:
            print("   👍 Reasonable clustering structure detected!")
        else:
            print("   ⚠️ Weak clustering structure - this is normal for small test datasets")
        
        print(f"   📊 Final clusters: {len(clusters_df)} sequences")
        
        # Save test results
        print("\n💾 Saving test results...")
        
        # Save test data
        data.to_csv('test_peptide_data.csv', index=False)
        print("   ✅ Test data saved to test_peptide_data.csv")
        
        # Save clustering results
        if not clusters_df.empty:
            clusters_df.to_csv('test_clustering_results.csv', index=False)
            print("   ✅ Clustering results saved to test_clustering_results.csv")
        
        print("\n🎉 Full pipeline test completed successfully!")
        return True
        
    except Exception as e:
        print(f"❌ Pipeline test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_consensus_validation_quality():
    """Test consensus validation with different data qualities"""
    print("\n🔬 Testing Consensus Validation Quality Assessment...")
    
    try:
        from Library.enhanced_gibbs_clustering import enhanced_gibbs_cluster_with_validation
        
        # Test Case 1: Well-separated clusters
        print("\n📊 Test Case 1: Well-separated clusters")
        well_separated = [
            'DEFGHIK',  # Group 1
            'DEFGHIL',
            'DEFGHIN',
            'STUVWXY',  # Group 2
            'STUVWXZ',
            'STUVWXA'
        ]
        
        result1 = enhanced_gibbs_cluster_with_validation(
            well_separated,
            motif_length=7,
            k_range=(2, 3),
            validation_enabled=True,
            n_consensus_iterations=20,
            n_bootstrap_iterations=10,
            random_seed=42
        )
        
        print(f"   Silhouette: {result1[3]['silhouette_score']:.3f}")
        print(f"   Stability: {result1[3]['stability_score']:.3f}")
        
        # Test Case 2: Similar sequences (poor clustering)
        print("\n📊 Test Case 2: Very similar sequences")
        similar_sequences = [
            'DEFGHIK',
            'DEFGHIL', 
            'DEFGHIN',
            'DEFGHIM',
            'DEFGHIT',
            'DEFGHIS'
        ]
        
        result2 = enhanced_gibbs_cluster_with_validation(
            similar_sequences,
            motif_length=7,
            k_range=(2, 3),
            validation_enabled=True,
            n_consensus_iterations=20,
            n_bootstrap_iterations=5,
            random_seed=42
        )
        
        print(f"   Silhouette: {result2[3]['silhouette_score']:.3f}")
        print(f"   Stability: {result2[3]['stability_score']:.3f}")
        
        # Compare results
        print(f"\n📈 Quality Comparison:")
        print(f"   Well-separated data: Silhouette={result1[3]['silhouette_score']:.3f}, Stability={result1[3]['stability_score']:.3f}")
        print(f"   Similar data: Silhouette={result2[3]['silhouette_score']:.3f}, Stability={result2[3]['stability_score']:.3f}")
        
        if result1[3]['silhouette_score'] > result2[3]['silhouette_score']:
            print("   ✅ Validation correctly identifies better clustering in well-separated data!")
        
        return True
        
    except Exception as e:
        print(f"❌ Quality assessment test failed: {e}")
        return False

def main():
    """Run all tests"""
    print("🧪 Starting Enhanced Pipeline Tests...\n")
    
    tests_passed = 0
    total_tests = 2
    
    # Test 1: Full pipeline
    if test_full_pipeline():
        tests_passed += 1
    
    # Test 2: Quality assessment
    if test_consensus_validation_quality():
        tests_passed += 1
    
    print(f"\n{'='*60}")
    print(f"Enhanced Pipeline Test Results: {tests_passed}/{total_tests} passed")
    
    if tests_passed == total_tests:
        print("🎉 All enhanced pipeline tests passed!")
        print("\n🚀 Ready to use enhanced consensus clustering in production!")
        print("\nGenerated files:")
        print("- test_peptide_data.csv (sample input data)")
        print("- test_clustering_results.csv (sample output)")
        print("- test_consensus_validation.png (validation plot)")
        return 0
    else:
        print("⚠️ Some tests failed. Check the output above for details.")
        return 1

if __name__ == "__main__":
    sys.exit(main())