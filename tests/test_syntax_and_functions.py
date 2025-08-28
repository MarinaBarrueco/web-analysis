#!/usr/bin/env python3
"""
Comprehensive test script to verify syntax and basic functionality
of the peptide analysis pipeline
"""

import sys
import traceback
import pandas as pd
import numpy as np
from pathlib import Path

def test_imports():
    """Test all critical imports"""
    try:
        # Core library imports
        from Library.lib import (
            group_by_peptide, filter_by_CPM, prepare_data, run_deseq2,
            significant_DE_peptides, run_vst, plot_heatmap, volcano_plot,
            constant_positions, remove_constant, box_plot, merge_data
        )
        
        # Clustering imports
        from Library.Clustering import gibbs_cluster
        from Library.GibbsClusterAdvanced import gibbs_cluster_advanced, GibbsClusterAdvanced
        
        # Validation and monitoring imports
        from Library.validation import (
            validate_dataframe, validate_peptide_sequences, 
            validate_experimental_design, validate_clustering_parameters,
            create_validation_summary
        )
        from Library.monitoring import AnalysisMonitor, safe_analysis_execution
        
        print("✅ All imports successful")
        return True
    except Exception as e:
        print(f"❌ Import failed: {e}")
        traceback.print_exc()
        return False

def test_data_validation():
    """Test data validation functions"""
    try:
        from Library.validation import validate_dataframe, validate_peptide_sequences
        
        # Test DataFrame validation
        df = pd.DataFrame({
            'peptide': ['ACDEFGHIK', 'LMNPQRSTU', 'VWXYZZABC'],
            'sample1': [10, 20, 30],
            'sample2': [15, 25, 35]
        })
        
        is_valid, errors = validate_dataframe(df, required_columns=['peptide'])
        assert is_valid, f"DataFrame validation failed: {errors}"
        
        # Test sequence validation
        sequences = ['ACDEFGHIK', 'ACDEFGHIL', 'INVALID123']
        valid_seqs, errors = validate_peptide_sequences(sequences, min_length=5)
        assert len(valid_seqs) == 2, f"Expected 2 valid sequences, got {len(valid_seqs)}"
        assert len(errors) == 1, f"Expected 1 error, got {len(errors)}"
        
        print("✅ Data validation tests passed")
        return True
    except Exception as e:
        print(f"❌ Data validation test failed: {e}")
        traceback.print_exc()
        return False

def test_clustering_functions():
    """Test clustering functions with minimal data"""
    try:
        from Library.Clustering import gibbs_cluster
        from Library.GibbsClusterAdvanced import GibbsClusterAdvanced
        
        # Test data - all sequences same length
        sequences = ['ACDEFGHI', 'ACDEFGHL', 'LMNPQRST', 'LMNPQRSV']
        
        # Test basic clustering
        clusters_df, pwms, logos = gibbs_cluster(
            sequences, 
            motif_length=8, 
            num_clusters=2,
            n_iter=10,  # Minimal iterations for testing
            return_logos=False  # Skip logo generation for speed
        )
        
        assert len(clusters_df) == 4, f"Expected 4 sequences, got {len(clusters_df)}"
        assert 'Sequence' in clusters_df.columns, "Missing Sequence column"
        assert 'Gn' in clusters_df.columns, "Missing Gn column"
        
        # Test advanced clustering (minimal config)
        advanced_clusterer = GibbsClusterAdvanced(
            motif_length=8,
            num_clusters=2,
            num_seeds=1,
            iterations=5,
            temperature_steps=5
        )
        
        advanced_clusterer.fit(sequences)
        clusters_df2 = advanced_clusterer.get_clusters()
        assert len(clusters_df2) == 4, f"Expected 4 sequences in advanced clustering, got {len(clusters_df2)}"
        
        print("✅ Clustering function tests passed")
        return True
    except Exception as e:
        print(f"❌ Clustering test failed: {e}")
        traceback.print_exc()
        return False

def test_library_functions():
    """Test core library functions"""
    try:
        from Library.lib import group_by_peptide, filter_by_CPM, constant_positions
        
        # Test data
        df = pd.DataFrame({
            'peptide': ['ACDEFGHIK', 'ACDEFGHIL', 'LMNPQRSTW', 'LMNPQRSTV'],
            'control1': [10, 15, 5, 8],
            'control2': [12, 18, 7, 10],
            'exp1': [20, 25, 15, 18],
            'exp2': [22, 28, 17, 20]
        })
        
        conditions = {
            'control1': 'Control',
            'control2': 'Control', 
            'exp1': 'Experiment',
            'exp2': 'Experiment'
        }
        
        # Test grouping
        grouped, summary = group_by_peptide(df, conditions, r'[A-Z]{9}')
        assert not grouped.empty, "Grouping resulted in empty DataFrame"
        assert 'Clean Peptide' in grouped.columns, "Missing Clean Peptide column"
        
        # Test CPM filtering
        filtered, fig = filter_by_CPM(grouped, conditions, cpm_thresh=1.0, min_samples=1)
        assert not filtered.empty, "CPM filtering resulted in empty DataFrame"
        
        # Test constant positions
        const_pos = constant_positions('A.C.{7}')
        assert isinstance(const_pos, list), "constant_positions should return a list"
        
        print("✅ Library function tests passed")
        return True
    except Exception as e:
        print(f"❌ Library function test failed: {e}")
        traceback.print_exc()
        return False

def test_monitoring():
    """Test monitoring system"""
    try:
        from Library.monitoring import AnalysisMonitor
        
        monitor = AnalysisMonitor()
        
        # Test step monitoring
        with monitor.monitor_step("test_step"):
            # Simulate some work
            data = pd.DataFrame({'a': range(100), 'b': range(100, 200)})
            quality = monitor.check_data_quality(data, "test_step")
            assert 'step' in quality, "Quality metrics missing step info"
        
        # Test report generation
        report = monitor.generate_safety_report()
        assert 'analysis_summary' in report, "Report missing analysis summary"
        assert 'performance_metrics' in report, "Report missing performance metrics"
        
        print("✅ Monitoring system tests passed")
        return True
    except Exception as e:
        print(f"❌ Monitoring test failed: {e}")
        traceback.print_exc()
        return False

def main():
    """Run all tests"""
    print("🧪 Running comprehensive syntax and functionality tests...\n")
    
    tests = [
        ("Import Tests", test_imports),
        ("Data Validation Tests", test_data_validation),
        ("Clustering Function Tests", test_clustering_functions),
        ("Library Function Tests", test_library_functions),
        ("Monitoring System Tests", test_monitoring)
    ]
    
    passed = 0
    failed = 0
    
    for test_name, test_func in tests:
        print(f"Running {test_name}...")
        try:
            if test_func():
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"❌ {test_name} crashed: {e}")
            failed += 1
        print()
    
    print("=" * 50)
    print(f"Test Results: {passed} passed, {failed} failed")
    
    if failed == 0:
        print("🎉 All tests passed! The pipeline is ready for use.")
        return 0
    else:
        print("⚠️ Some tests failed. Please review the errors above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())