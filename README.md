# 🧬 Peptide Analysis Dashboard

A comprehensive interactive pipeline for differential expression analysis and clustering of peptide sequences from library screening experiments. Features an intuitive Streamlit dashboard with automated quality control, consensus clustering validation, and publication-ready visualizations.

## 🚀 Quick Start

### Installation
1. Clone or download this repository
2. Install dependencies:
```bash
pip install -r requirements.txt
```

### Launch Dashboard
```bash
streamlit run main_dashboard.py
```

The dashboard will open in your web browser with an interactive interface featuring:
- **Step-by-step workflow guidance**
- **Built-in help and documentation**
- **Real-time analysis validation**  
- **Publication-ready visualizations**

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

## 🧬 Dashboard Workflow

The interactive dashboard guides users through a complete analysis pipeline:

### 📁 **Step 1: Data Upload**
- Upload CSV files with peptide sequences and count data
- Automatic data validation and format checking
- Preview data structure and quality metrics

### ⚙️ **Step 2: Configure Parameters**
- **Simplified Clustering Modes**:
  - **Simple**: Basic clustering with fixed K (recommended for beginners)
  - **Advanced**: Enhanced features with fixed K
  - **Auto-Selection**: Automatic optimal K discovery
- **Pattern Matching**: Flexible regex for library designs
- **Quality Filters**: CPM thresholds and sample requirements

### 🔬 **Step 3: Experimental Design**
- Interactive sample assignment (Control/Experiment/Exclude)
- Real-time validation with warnings and suggestions
- Automatic design balance checking

### 📊 **Step 4: Analysis Execution**
1. **Pattern Matching & Filtering**: Quality-controlled sequence processing
2. **Differential Expression**: DESeq2-based statistical analysis with multiple testing correction
3. **Significance Assessment**: Automated classification of up/down/non-significant peptides  
4. **Clustering Analysis**: Advanced Gibbs sampling with consensus validation
5. **Quality Metrics**: Silhouette scoring, stability assessment, and bootstrap validation

### 📈 **Step 5: Results Exploration**
- **Interactive Tabs**: Filtering, Differential Expression, Clustering, Summary
- **Publication-Ready Plots**: MA plots, volcano plots, box plots by cluster
- **Enhanced Visualizations**: Clean sequence logos showing full peptide motifs
- **Quality Dashboards**: Comprehensive clustering metrics and interpretations

## 🔧 Key Features

### 🎯 **User Experience**
- **Intuitive Dashboard**: Step-by-step guided workflow with progress tracking
- **Built-in Help**: Comprehensive instructions and tooltips throughout
- **Real-time Validation**: Immediate feedback on data quality and parameter choices
- **Error Handling**: Graceful failure modes with helpful error messages

### 📊 **Analysis Capabilities**
- **Flexible Input**: Support for various CSV formats and experimental designs
- **Robust Statistics**: DESeq2-based differential expression with proper normalization
- **Advanced Clustering**: Three modes (Simple/Advanced/Auto-Selection) for different user needs
- **Quality Metrics**: Comprehensive clustering validation with interpretive guidance

### 🎨 **Visualizations**
- **Clean Sequence Logos**: Noise-filtered motifs showing full peptide patterns
- **Interactive Plots**: MA plots, volcano plots, and differential expression box plots
- **Quality Dashboards**: Clustering metrics with color-coded assessments
- **Publication Ready**: High-resolution figures with proper styling

### 🔬 **Scientific Rigor**
- **Consensus Validation**: Subsampling-based stability assessment
- **Bootstrap Analysis**: Additional robustness testing for clustering
- **Multiple Testing Correction**: Proper FDR control for differential expression
- **Comprehensive Documentation**: Theoretical foundations and best practices

## 💾 **Input Data Format**

Your CSV file should have this structure:

```csv
peptide,Sample1_Control,Sample2_Control,Sample1_Treatment,Sample2_Treatment
AYCPFRSWPGCGG,1250,890,2340,1980
AFCSLWRSGDCGG,890,1100,0,45
AACALWRSTPCGG,2100,1850,3200,2890
```

**Requirements:**
- First column: `peptide` containing sequence strings
- Remaining columns: Numeric count data for each sample
- Column names can be anything descriptive for your experiment
- Missing values should be 0 or empty (will be treated as 0)

## 🎛️ **Dashboard Interface**

### Getting Started
1. Launch the dashboard: `streamlit run main_dashboard.py`
2. Click **"📖 How to Use This Dashboard"** for detailed instructions
3. Upload your CSV file and follow the guided workflow
4. Download results as CSV files or publication-ready plots

### Key Interface Elements
- **Progress Tracker**: Shows your current step in the analysis pipeline
- **Parameter Validation**: Live feedback on parameter compatibility
- **Interactive Results**: Tabs for filtering, differential expression, clustering, and summary
- **Quality Assessments**: Color-coded metrics with interpretive guidance
- **Download Options**: Export results and visualizations

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