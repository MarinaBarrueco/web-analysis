# Advanced Gibbs Clustering Implementation

This implementation provides a comprehensive Python version of the GibbsCluster-2.0 algorithm for peptide sequence clustering, based on the original Perl implementation (`GibbsCluster-2.0e_SA.pl`).

## Features

### Core Algorithm Features
- **Gibbs Sampling**: Probabilistic clustering using Gibbs sampling with simulated annealing
- **Multiple Seeds**: Run multiple random initializations for better convergence
- **Temperature Annealing**: Gradual cooling schedule for optimal clustering
- **Kullback-Leibler Divergence**: Quality metric for cluster assessment
- **Sequence Weighting**: Different weighting schemes for sequences
- **Background Models**: Multiple background frequency models (flat, pre-calculated, data-derived)

### Advanced Features
- **Cluster Range Testing**: Automatically test different numbers of clusters
- **Trash Cluster**: Option to remove outlier sequences
- **Insertion/Deletion Support**: Handle variable-length sequences (framework implemented)
- **Interactive Web Interface**: Streamlit-based GUI for easy use
- **Comprehensive Reporting**: Detailed results with visualizations

### Output Features
- **Sequence Logos**: Automatic generation using logomaker
- **Position Weight Matrices**: Export PWMs for each cluster
- **Quality Metrics**: KLD scores and cluster statistics
- **Multiple Export Formats**: CSV, JSON, PNG logos, ZIP archives

## File Structure

```
/Users/jq23948/Desktop/web-analysis/
├── Library/
│   ├── GibbsClusterAdvanced.py    # Main advanced implementation
│   ├── Clustering.py              # Original simple implementation
│   └── lib.py                     # Supporting functions
├── main.py                        # Integrated analysis pipeline
├── clustering_page.py             # Standalone clustering interface
├── requirements.txt               # Python dependencies
└── README_GibbsClustering.md     # This documentation
```

## Usage

### 1. Integrated Pipeline (main.py)

The advanced clustering is integrated into the main differential expression analysis pipeline:

```bash
streamlit run main.py
```

Features in integrated mode:
- Upload peptide expression data
- Perform differential expression analysis
- Run advanced clustering on upregulated peptides
- Compare clustering quality across different cluster numbers

### 2. Standalone Clustering (clustering_page.py)

For dedicated clustering analysis:

```bash
streamlit run clustering_page.py
```

Features in standalone mode:
- Upload sequences via CSV, text, or FASTA files
- Paste sequences directly
- Use example data for testing
- Full parameter control
- Comprehensive result download options

### 3. Programmatic Usage

```python
from Library.GibbsClusterAdvanced import gibbs_cluster_advanced, GibbsClusterAdvanced

# Simple function interface
sequences = ["ACDEFGHIK", "LMNPQRSTU", "VWXYZABCD"]
clusters_df, pwms, logos = gibbs_cluster_advanced(
    sequences,
    motif_length=9,
    num_clusters=3,
    num_seeds=5
)

# Advanced class interface
clusterer = GibbsClusterAdvanced(
    motif_length=9,
    num_clusters=(2, 5),  # Test 2-5 clusters
    num_seeds=3,
    iterations=15,
    temperature_start=2.0,
    lambda_penalty=0.8
)

clusterer.fit(sequences)
summary = clusterer.get_summary()
```

## Parameters

### Basic Parameters
- **motif_length**: Length of the conserved motif to identify
- **num_clusters**: Number of clusters (int) or range (tuple)
- **num_seeds**: Number of random initializations
- **iterations**: Gibbs sampling iterations per temperature step
- **temperature_start**: Initial temperature for simulated annealing
- **temperature_steps**: Number of temperature reduction steps

### Advanced Parameters
- **lambda_penalty**: Penalty for inter-cluster similarity (0.0-2.0)
- **sigma_weight**: Weight factor for small clusters (0.1-20.0)
- **sequence_weighting**: Weighting scheme (0=1/ns, 1=clustering, 2=none)
- **background_model**: Background frequencies (0=flat, 1=pre-calculated, 2=from data)
- **use_trash_cluster**: Enable outlier removal
- **trash_threshold**: Threshold for outlier classification

## Algorithm Details

### Gibbs Sampling Process
1. **Initialize**: Random cluster assignments
2. **Sample**: For each sequence, calculate cluster probabilities
3. **Update**: Reassign sequence to new cluster based on probabilities
4. **Anneal**: Gradually reduce temperature for fine-tuning
5. **Converge**: Repeat until stable clustering

### Probability Calculation
For each sequence and cluster combination:
```
P(cluster|sequence) ∝ P(cluster) × P(sequence|cluster)
```

Where:
- P(cluster) includes size bias and penalty terms
- P(sequence|cluster) uses position-specific amino acid frequencies
- Temperature scaling allows exploration vs exploitation trade-off

### Quality Assessment
- **KLD Score**: Measures information content of clusters
- **Cluster Utilization**: Percentage of non-empty clusters
- **Size Distribution**: Balance of cluster sizes

## Comparison with Original

### Implemented Features from Perl Version
✅ **Core algorithm**: Gibbs sampling with simulated annealing  
✅ **Multiple seeds**: Better convergence through random restarts  
✅ **Temperature schedule**: Logarithmic cooling  
✅ **KLD calculation**: Cluster quality assessment  
✅ **Sequence weighting**: Multiple weighting schemes  
✅ **Background models**: Pre-calculated and data-derived frequencies  
✅ **Result export**: Comprehensive output formats  

### Enhanced Features
✅ **Web interface**: User-friendly Streamlit GUI  
✅ **Integrated pipeline**: Combined with differential expression  
✅ **Modern visualizations**: Matplotlib/logomaker graphics  
✅ **Batch processing**: Handle multiple cluster numbers  
✅ **Progress tracking**: Real-time status updates  

### Framework for Future Features
🔄 **Insertion/deletion handling**: Structure implemented, algorithm pending  
🔄 **Trash cluster**: Basic implementation, needs refinement  
🔄 **Phase shift moves**: Framework for variable positioning  
🔄 **Parallel processing**: Multi-threading support structure  

## Performance Considerations

### Memory Usage
- Scales with: sequences × motif_length × num_clusters × amino_acids (20)
- Typical usage: 1000 sequences × 9 positions × 4 clusters × 20 AA = ~720KB

### Runtime
- Depends on: num_seeds × temperature_steps × iterations × sequences
- Typical runtime: 2-30 seconds for 100-1000 sequences

### Recommendations
- Start with 3-5 seeds for exploration
- Use 10-20 iterations for most datasets
- Increase temperature_steps for difficult clustering problems
- Test cluster ranges (2-5) to find optimal number

## Installation

1. **Clone/download** the web-analysis directory
2. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```
3. **Run application**:
   ```bash
   streamlit run main.py          # Integrated pipeline
   streamlit run clustering_page.py  # Standalone clustering
   ```

## Troubleshooting

### Common Issues
1. **Import errors**: Ensure all dependencies are installed
2. **Memory errors**: Reduce num_seeds or cluster range for large datasets
3. **Convergence issues**: Increase iterations or temperature_steps
4. **Logo generation fails**: Check logomaker installation

### Performance Tips
1. Start with smaller datasets to test parameters
2. Use single cluster number before testing ranges
3. Monitor memory usage with large sequence sets
4. Save intermediate results for long-running analyses

## Citation

If you use this implementation, please cite the original GibbsCluster paper:

> Andreatta, M., Alvarez, B., & Nielsen, M. (2017). GibbsCluster: unsupervised clustering and alignment of peptide sequences. Nucleic acids research, 45(W1), W458-W463.

## Support

For issues with this implementation:
1. Check parameter settings and input format
2. Review error messages and logs
3. Test with example data first
4. Report bugs with minimal reproducible examples