# Peptide Analysis Pipeline

A comprehensive pipeline for differential expression analysis and clustering of peptide sequences from library screening experiments.

## 🚀 Quick Start

### Requirements
Install dependencies:
```bash
pip install -r requirements.txt
```

### Run Analysis
- **Dashboard Interface**: `streamlit run main_dashboard.py`
- **Basic Interface**: `streamlit run main.py`

## 📁 Project Structure

```
web-analysis/
├── README.md                 # This file
├── requirements.txt          # Python dependencies
├── main.py                   # Basic Streamlit interface
├── main_dashboard.py         # Enhanced dashboard interface
├── Library/                  # Core analysis modules
│   ├── lib.py               # Main analysis functions
│   ├── Clustering.py        # Basic Gibbs clustering
│   ├── GibbsClusterAdvanced.py  # Advanced clustering
│   ├── consensus_clustering.py  # Validation methods
│   ├── enhanced_gibbs_clustering.py  # Integrated clustering
│   ├── validation.py        # Input validation
│   └── monitoring.py        # Analysis monitoring
├── docs/                     # Documentation
│   ├── THEORETICAL_EXPLANATION.md     # Mathematical foundations
│   ├── ANALYSIS_OPTIONS_GUIDE.md      # Parameter documentation
│   ├── BIOCHEMICAL_ANALYSIS_REVIEW.md # Scientific review
│   └── README_GibbsClustering.md      # Gibbs clustering details
├── data/                     # Sample data and generators
├── outputs/                  # Analysis results
├── tests/                    # Test scripts and validation
├── archive/                  # Legacy files
└── gibbscluster-2.0/        # External Gibbs clustering tool
```

## 🧬 Analysis Pipeline

### 1. Data Processing
- **Pattern Matching**: Extract peptides matching library design (e.g., `A.C.{7}C`)
- **CPM Filtering**: Remove low-abundance sequences using counts per million
- **Quality Control**: Validate data integrity and experimental design

### 2. Differential Expression
- **DESeq2 Analysis**: Statistical testing for condition differences
- **Multiple Testing**: Benjamini-Hochberg FDR correction
- **Effect Size**: Log2 fold change thresholds for biological significance

### 3. Clustering Analysis
- **Gibbs Sampling**: Discover conserved sequence motifs
- **Consensus Validation**: Stability assessment through subsampling
- **Quality Metrics**: Silhouette analysis and bootstrap validation
- **Visualization**: Generate sequence logos and cluster plots

## 🔧 Key Features

- **Flexible Input**: CSV files with peptide sequences and count data
- **Robust Statistics**: DESeq2-based differential expression analysis
- **Advanced Clustering**: Consensus methods for reliable motif discovery
- **Quality Control**: Comprehensive validation and error checking
- **Interactive Dashboard**: User-friendly Streamlit interface
- **Comprehensive Documentation**: Theoretical foundations and practical guides

## 📊 Usage Examples

### Basic Analysis
```python
from Library.lib import *

# Load and process data
data = pd.read_csv("peptides.csv")
grouped, summary = group_by_peptide(data, conditions, regex="A.C.{7}C")
filtered, fig = filter_by_CPM(grouped, conditions, cpm_threshold=5, min_count=2)

# Differential expression
count_data, meta_data = prepare_data(filtered, conditions)
results_df, ma_plot = run_deseq2(count_data, meta_data)
```

### Advanced Clustering
```python
from Library.enhanced_gibbs_clustering import enhanced_gibbs_cluster_with_validation

# Run clustering with validation
clusters_df, pwms, logos, validation = enhanced_gibbs_cluster_with_validation(
    sequences, 
    motif_length=8,
    k_range=(2, 6),
    validation_enabled=True
)
```

## 📚 Documentation

- **[Theoretical Explanation](docs/THEORETICAL_EXPLANATION.md)**: Mathematical foundations and statistical methods
- **[Analysis Options Guide](docs/ANALYSIS_OPTIONS_GUIDE.md)**: Detailed parameter documentation
- **[Biochemical Review](docs/BIOCHEMICAL_ANALYSIS_REVIEW.md)**: Scientific assessment and best practices

## 🧪 Testing

Run validation tests:
```bash
python tests/test_enhanced_pipeline.py
python tests/test_consensus_clustering.py
```

## 📈 Results

Analysis generates:
- **Differential expression results**: Fold changes, p-values, significance classifications
- **Cluster assignments**: Sequence groupings with motif information
- **Quality metrics**: Validation scores and confidence assessments
- **Visualizations**: MA plots, volcano plots, heatmaps, sequence logos
- **Reports**: Comprehensive analysis summaries

## ⚙️ Configuration

Key parameters can be adjusted:
- **CPM threshold**: Minimum abundance for sequence retention (default: 5)
- **Fold change threshold**: Minimum effect size for significance (default: 1.5)
- **Cluster number**: Fixed or range-based selection (default: auto-selection)
- **Validation iterations**: Consensus and bootstrap repetitions (default: 50/30)

## 📄 License

This project is for research use. Please cite appropriately if used in publications.

## 🤝 Contributing

For questions or contributions, please refer to the documentation in the `docs/` folder.