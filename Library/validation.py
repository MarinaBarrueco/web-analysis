"""
Data Validation and Safety Utilities for Peptide Analysis Pipeline
"""

import pandas as pd
import numpy as np
import re
from typing import List, Dict, Union, Tuple, Optional

# Valid amino acid alphabet
VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

def validate_dataframe(df: pd.DataFrame, required_columns: List[str] = None, 
                      min_rows: int = 1) -> Tuple[bool, List[str]]:
    """
    Validate a DataFrame for basic requirements.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame to validate
    required_columns : list
        List of required column names
    min_rows : int
        Minimum number of rows required
        
    Returns:
    --------
    tuple: (is_valid, error_messages)
    """
    errors = []
    
    if df is None:
        errors.append("DataFrame is None")
        return False, errors
    
    if df.empty:
        errors.append("DataFrame is empty")
        return False, errors
    
    if len(df) < min_rows:
        errors.append(f"DataFrame has {len(df)} rows, minimum required: {min_rows}")
    
    if required_columns:
        missing_cols = [col for col in required_columns if col not in df.columns]
        if missing_cols:
            errors.append(f"Missing required columns: {missing_cols}")
    
    return len(errors) == 0, errors

def validate_peptide_sequences(sequences: Union[List[str], pd.Series], 
                             min_length: int = 1, 
                             max_length: int = 100) -> Tuple[List[str], List[str]]:
    """
    Validate and clean peptide sequences.
    
    Parameters:
    -----------
    sequences : list or pd.Series
        Input sequences to validate
    min_length : int
        Minimum sequence length
    max_length : int
        Maximum sequence length
        
    Returns:
    --------
    tuple: (valid_sequences, error_messages)
    """
    if isinstance(sequences, pd.Series):
        sequences = sequences.tolist()
    
    valid_seqs = []
    errors = []
    
    for i, seq in enumerate(sequences):
        if not isinstance(seq, str):
            errors.append(f"Sequence {i} is not a string: {type(seq)}")
            continue
        
        # Clean sequence
        seq = seq.upper().strip()
        
        # Check length
        if len(seq) < min_length:
            errors.append(f"Sequence {i} too short: {len(seq)} < {min_length}")
            continue
        
        if len(seq) > max_length:
            errors.append(f"Sequence {i} too long: {len(seq)} > {max_length}")
            continue
        
        # Check for valid amino acids
        invalid_chars = set(seq) - VALID_AMINO_ACIDS
        if invalid_chars:
            errors.append(f"Sequence {i} contains invalid characters: {invalid_chars}")
            continue
        
        valid_seqs.append(seq)
    
    return valid_seqs, errors

def validate_experimental_design(df: pd.DataFrame, 
                                conditions: Dict[str, str]) -> Tuple[bool, List[str]]:
    """
    Validate experimental design setup.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input data
    conditions : dict
        Mapping of column names to conditions
        
    Returns:
    --------
    tuple: (is_valid, error_messages)
    """
    errors = []
    
    # Check if condition columns exist in dataframe
    missing_cols = [col for col in conditions.keys() if col not in df.columns]
    if missing_cols:
        errors.append(f"Condition columns not found in data: {missing_cols}")
        return False, errors
    
    # Check for valid condition labels
    valid_conditions = {"Control", "Experiment"}
    invalid_conditions = set(conditions.values()) - valid_conditions
    if invalid_conditions:
        errors.append(f"Invalid condition labels: {invalid_conditions}. Must be 'Control' or 'Experiment'")
    
    # Check for balanced design
    condition_counts = pd.Series(list(conditions.values())).value_counts()
    
    if "Control" not in condition_counts:
        errors.append("No Control samples specified")
    elif condition_counts["Control"] < 2:
        errors.append(f"Only {condition_counts['Control']} Control sample(s). Recommend at least 2 for statistics")
    
    if "Experiment" not in condition_counts:
        errors.append("No Experiment samples specified")
    elif condition_counts["Experiment"] < 2:
        errors.append(f"Only {condition_counts['Experiment']} Experiment sample(s). Recommend at least 2 for statistics")
    
    # Check for non-negative values in condition columns
    for col in conditions.keys():
        if (df[col] < 0).any():
            errors.append(f"Column {col} contains negative values")
    
    return len(errors) == 0, errors

def validate_clustering_parameters(num_clusters: Union[int, Tuple[int, int]], 
                                 motif_length: int,
                                 num_sequences: int) -> Tuple[bool, List[str]]:
    """
    Validate clustering parameters.
    
    Parameters:
    -----------
    num_clusters : int or tuple
        Number of clusters or range
    motif_length : int  
        Motif length
    num_sequences : int
        Number of input sequences
        
    Returns:
    --------
    tuple: (is_valid, error_messages)
    """
    errors = []
    
    # Validate motif length
    if motif_length <= 0:
        errors.append("Motif length must be positive")
    elif motif_length > 50:
        errors.append("Motif length unusually large (>50). This may cause memory issues")
    
    # Validate number of clusters
    if isinstance(num_clusters, tuple):
        min_clusters, max_clusters = num_clusters
        if min_clusters <= 0:
            errors.append("Minimum clusters must be positive")
        if max_clusters <= 0:
            errors.append("Maximum clusters must be positive")
        if min_clusters > max_clusters:
            errors.append("Minimum clusters cannot exceed maximum clusters")
        if max_clusters > num_sequences:
            errors.append(f"Maximum clusters ({max_clusters}) exceeds number of sequences ({num_sequences})")
    else:
        if num_clusters <= 0:
            errors.append("Number of clusters must be positive")
        if num_clusters > num_sequences:
            errors.append(f"Number of clusters ({num_clusters}) exceeds number of sequences ({num_sequences})")
    
    # Check for reasonable clustering
    if num_sequences < 5:
        errors.append("Very few sequences for clustering (<5). Results may not be meaningful")
    
    return len(errors) == 0, errors

def validate_regex_pattern(pattern: str) -> Tuple[bool, Optional[str]]:
    """
    Validate a regular expression pattern.
    
    Parameters:
    -----------
    pattern : str
        Regular expression pattern to validate
        
    Returns:
    --------
    tuple: (is_valid, error_message)
    """
    try:
        re.compile(pattern)
        return True, None
    except re.error as e:
        return False, f"Invalid regex pattern: {e}"

def check_data_quality(df: pd.DataFrame, columns: List[str]) -> Dict[str, any]:
    """
    Perform basic data quality checks.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input data
    columns : list
        Columns to check
        
    Returns:
    --------
    dict: Quality metrics and warnings
    """
    quality_report = {
        "total_rows": len(df),
        "warnings": [],
        "statistics": {}
    }
    
    for col in columns:
        if col not in df.columns:
            quality_report["warnings"].append(f"Column {col} not found")
            continue
        
        col_data = df[col]
        stats = {
            "missing_values": col_data.isna().sum(),
            "zero_values": (col_data == 0).sum(),
            "negative_values": (col_data < 0).sum() if col_data.dtype in ['int64', 'float64'] else 0,
            "unique_values": col_data.nunique(),
            "data_type": str(col_data.dtype)
        }
        
        if col_data.dtype in ['int64', 'float64']:
            stats.update({
                "mean": col_data.mean(),
                "median": col_data.median(),
                "std": col_data.std(),
                "min": col_data.min(),
                "max": col_data.max()
            })
        
        quality_report["statistics"][col] = stats
        
        # Generate warnings
        if stats["missing_values"] > 0:
            quality_report["warnings"].append(f"Column {col} has {stats['missing_values']} missing values")
        
        if stats["negative_values"] > 0:
            quality_report["warnings"].append(f"Column {col} has {stats['negative_values']} negative values")
        
        if col_data.dtype in ['int64', 'float64'] and stats["std"] == 0:
            quality_report["warnings"].append(f"Column {col} has zero variance")
    
    return quality_report

def safe_numeric_operation(operation_func, *args, default_value=0, **kwargs):
    """
    Safely perform numeric operations with error handling.
    
    Parameters:
    -----------
    operation_func : callable
        Function to execute
    *args : arguments for the function
    default_value : any
        Value to return if operation fails
    **kwargs : keyword arguments for the function
        
    Returns:
    --------
    Result of operation or default_value if failed
    """
    try:
        result = operation_func(*args, **kwargs)
        
        # Check for invalid results
        if isinstance(result, (int, float)):
            if np.isnan(result) or np.isinf(result):
                return default_value
        elif isinstance(result, (pd.Series, pd.DataFrame, np.ndarray)):
            if result.empty if hasattr(result, 'empty') else len(result) == 0:
                return default_value
        
        return result
    except Exception as e:
        print(f"Warning: Numeric operation failed: {e}. Using default value: {default_value}")
        return default_value

def create_validation_summary(df: pd.DataFrame, conditions: Dict[str, str], 
                            pattern: str) -> Dict[str, any]:
    """
    Create a comprehensive validation summary for the analysis pipeline.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input data
    conditions : dict
        Experimental conditions
    pattern : str
        Regex pattern for peptide matching
        
    Returns:
    --------  
    dict: Comprehensive validation summary
    """
    summary = {
        "timestamp": pd.Timestamp.now(),
        "data_validation": {},
        "sequence_validation": {},
        "experimental_design": {},
        "pattern_validation": {},
        "overall_status": "UNKNOWN",
        "critical_errors": [],
        "warnings": []
    }
    
    # Basic data validation
    is_valid, errors = validate_dataframe(df, required_columns=["peptide"], min_rows=1)
    summary["data_validation"] = {"valid": is_valid, "errors": errors}
    if not is_valid:
        summary["critical_errors"].extend(errors)
    
    # Experimental design validation
    if is_valid:
        is_valid, errors = validate_experimental_design(df, conditions)
        summary["experimental_design"] = {"valid": is_valid, "errors": errors}
        if not is_valid:
            summary["critical_errors"].extend(errors)
    
    # Pattern validation
    is_valid, error = validate_regex_pattern(pattern)
    summary["pattern_validation"] = {"valid": is_valid, "error": error}
    if not is_valid:
        summary["critical_errors"].append(error)
    
    # Sequence validation (if applicable)
    if "peptide" in df.columns:
        valid_seqs, seq_errors = validate_peptide_sequences(df["peptide"], min_length=5)
        summary["sequence_validation"] = {
            "total_sequences": len(df),
            "valid_sequences": len(valid_seqs),
            "invalid_sequences": len(seq_errors),
            "errors": seq_errors[:10]  # Limit errors shown
        }
        
        if len(valid_seqs) == 0:
            summary["critical_errors"].append("No valid peptide sequences found")
        elif len(seq_errors) > len(df) * 0.5:
            summary["warnings"].append(f"High proportion of invalid sequences: {len(seq_errors)}/{len(df)}")
    
    # Data quality check
    if conditions:
        quality_report = check_data_quality(df, list(conditions.keys()) + ["peptide"])
        summary["data_quality"] = quality_report
        summary["warnings"].extend(quality_report["warnings"])
    
    # Overall status
    if summary["critical_errors"]:
        summary["overall_status"] = "FAILED"
    elif summary["warnings"]:
        summary["overall_status"] = "WARNING"
    else:
        summary["overall_status"] = "PASSED"
    
    return summary