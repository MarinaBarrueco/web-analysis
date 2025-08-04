# Analysis Options Guide - Detailed Parameter Documentation

## 🎯 Overview

This guide provides comprehensive documentation of all analysis parameters and options in the peptide analysis pipeline. Each parameter is explained with its biological significance, recommended values, and practical considerations.

---

## 📋 Section 1: Data Input and Validation

### 1.1 CSV File Requirements

**Parameter**: Input CSV file
**Description**: Primary data file containing peptide sequences and count data
**Format Requirements**:
- **First column**: Must be named 'peptide' and contain peptide sequences
- **Additional columns**: Count data for different samples/conditions
- **Encoding**: UTF-8 recommended
- **Headers**: First row must contain column names

**Example Structure**:
```csv
peptide,Control_1,Control_2,Experiment_1,Experiment_2
AVCDEFGHYC,25,30,45,50
AICDEFGHYC,20,22,38,42
...
```

**Common Issues**:
- Missing or misnamed 'peptide' column
- Non-numeric count data
- Empty cells (handled automatically as 0)
- Invalid amino acid characters (filtered out with warnings)

---

## 🔍 Section 2: Pattern Recognition Parameters

### 2.1 Peptide Library Pattern (Regular Expression)

**Parameter**: `regex` (default: `A.C.{7}C`)
**Description**: Regular expression defining the expected peptide library structure
**Purpose**: Filters and extracts peptides matching the designed library format

**Syntax Guide**:
- **Literal characters**: `A`, `C`, `D` = specific amino acids required
- **Wildcard**: `.` = any natural amino acid (A-Z except B, J, O, U, X, Z)
- **Quantifiers**: `{n}` = exactly n repetitions, `{n,m}` = n to m repetitions
- **Character classes**: `[AILMFWV]` = any hydrophobic amino acid

**Common Patterns**:
- `A.C.{7}C` = A, any AA, C, 7 random AAs, C (11 total length)
- `C.{8}C` = C, 8 random AAs, C (10 total length)
- `[FWY].{6}[FWY]` = Aromatic, 6 random, aromatic (8 total length)
- `K.{4}K` = Lysine, 4 random, lysine (6 total length)

**Validation**:
- **Success Rate**: Percentage of reads matching pattern (displayed in results)
- **Warning Signs**: Success rate <50% may indicate pattern mismatch
- **Troubleshooting**: Check actual sequences vs. expected pattern

### 2.2 Amino Acid Substitutions

**Parameter**: Automatic `*` → `Q` substitution
**Description**: Replaces stop codons (*) with glutamine (Q)
**Rationale**: Stop codons from sequencing errors converted to chemically similar amino acid
**Alternative**: Can be modified in code for other substitution strategies

---

## 📊 Section 3: Filtering and Normalization Parameters

### 3.1 CPM (Counts Per Million) Filtering

**Parameter**: `cpm_threshold` (default: 5)
**Description**: Minimum CPM value required for peptide retention
**Range**: 0-1000 (typical: 1-20)
**Purpose**: Removes low-abundance peptides likely representing noise

**Biological Interpretation**:
- **CPM = 5**: Peptide represents ≥0.0005% of library
- **Higher thresholds**: More stringent, fewer false positives
- **Lower thresholds**: More permissive, may include noise

**Recommendations by Application**:
- **High-stringency screens**: CPM ≥ 10
- **Discovery screens**: CPM ≥ 1-5  
- **Rare peptide analysis**: CPM ≥ 0.5-1
- **Quality control**: Monitor distribution before/after filtering

### 3.2 Minimum Sample Requirement

**Parameter**: `min_samples` (default: 1)
**Description**: Minimum number of samples that must exceed CPM threshold
**Range**: 1 to total number of samples
**Purpose**: Ensures consistent detection across replicates

**Selection Guidelines**:
- **Single replicate**: min_samples = 1
- **Biological replicates**: min_samples = 2-3 (depending on total replicates)
- **Technical + biological**: Consider total replicate structure
- **Conservative approach**: Require detection in majority of samples

**Statistical Impact**:
- **Higher requirements**: Reduced false positives, increased confidence
- **Lower requirements**: Higher sensitivity, potential noise inclusion
- **Balance**: Consider experimental goals and statistical power

### 3.3 Plot Filtering Visualization

**Parameter**: `plot_filtering` (default: True)
**Description**: Generate before/after filtering comparison plots
**Visualizations**:
- **Count distributions**: Histograms of raw counts
- **CPM distributions**: Normalized count distributions
- **Correlation plots**: Sample-to-sample relationships

**Interpretation**:
- **Bimodal distributions**: Clear separation of signal vs. noise
- **Correlations**: High correlation indicates good replicate consistency
- **Outliers**: May indicate technical issues or biological variation

---

## 🧮 Section 4: Experimental Design Parameters

### 4.1 Condition Assignment

**Parameter**: Column condition mapping
**Options**: "Control", "Experiment", "Exclude"
**Description**: Defines experimental groups for statistical comparison

**Design Considerations**:
- **Minimum replicates**: ≥2 per condition (≥3 recommended)
- **Balanced design**: Equal sample sizes preferred but not required
- **Controls**: Must include appropriate baseline conditions
- **Exclusions**: Remove problematic samples without losing data

**Statistical Implications**:
- **Underpowered designs**: May miss true differences
- **Unbalanced designs**: Handled by DESeq2 but may reduce power
- **Batch effects**: Consider including batch as covariate if relevant

---

## 🔬 Section 5: Clustering Parameters

### 5.1 Basic vs. Advanced Clustering

**Parameter**: `use_advanced_clustering` (default: True)
**Description**: Choose between basic and advanced Gibbs clustering algorithms

**Basic Clustering Features**:
- **Simple implementation**: Fast, minimal parameters
- **Fixed parameters**: Default Dirichlet priors, standard sampling
- **Output**: Cluster assignments, basic PWMs, sequence logos

**Advanced Clustering Features**:
- **Multiple seeds**: Improved convergence assessment
- **Temperature annealing**: Better global optimization
- **Parameter tuning**: Extensive customization options
- **Quality metrics**: KLD scores, convergence monitoring

**Selection Criteria**:
- **Exploratory analysis**: Basic clustering sufficient
- **Publication quality**: Advanced clustering recommended
- **Large datasets**: Advanced clustering more robust
- **Parameter optimization**: Advanced clustering necessary

### 5.2 Cluster Number Selection

#### 5.2.1 Manual Specification

**Parameter**: `num_clusters` (default: 4)
**Range**: 2 to N/2 (where N = number of sequences)
**Description**: Fixed number of clusters for analysis

**Selection Guidelines**:
- **Prior knowledge**: Use known number of binding sites/motifs
- **Exploratory**: Start with 3-5 clusters
- **Large datasets**: May support more clusters (up to 10-15)
- **Small datasets**: Limit to 2-4 clusters

#### 5.2.2 Range-Based Search

**Parameters**: `min_clusters`, `max_clusters`
**Description**: Test multiple k values and select optimal
**Advantage**: Data-driven selection, comprehensive evaluation
**Disadvantage**: Increased computational time

**Recommended Ranges**:
- **Small datasets** (<50 sequences): 2-5 clusters
- **Medium datasets** (50-200 sequences): 2-8 clusters
- **Large datasets** (>200 sequences): 2-12 clusters

### 5.3 Advanced Clustering Parameters

#### 5.3.1 Sampling Parameters

**Parameter**: `num_seeds` (default: 3)
**Description**: Number of independent sampling runs
**Purpose**: Assess convergence and find global optimum
**Range**: 1-10 (typical: 3-5)

**Parameter**: `iterations` (default: 10)
**Description**: Gibbs sampling iterations per temperature step
**Purpose**: Sufficient sampling at each temperature
**Range**: 5-50 (typical: 10-20)

**Parameter**: `motif_length` (default: 8)
**Description**: Length of conserved motif to discover
**Purpose**: Defines resolution of pattern discovery
**Considerations**: Should match biological expectations

#### 5.3.2 Temperature Annealing

**Parameter**: `temperature_start` (default: 1.5)
**Description**: Initial sampling temperature
**Purpose**: High temperature allows exploration of solution space
**Range**: 0.5-3.0 (typical: 1.0-2.0)

**Parameter**: `temperature_steps` (default: 20)
**Description**: Number of temperature reduction steps
**Purpose**: Gradual cooling for better convergence
**Range**: 10-50 (typical: 15-25)

**Annealing Schedule**: Temperature decreases exponentially
```
T(step) = T_start × (T_final/T_start)^(step/total_steps)
```

#### 5.3.3 Regularization Parameters

**Parameter**: `lambda_penalty` (default: 0.8)
**Description**: Penalty for similarity between clusters
**Purpose**: Encourages discovery of distinct motifs
**Range**: 0.0-2.0 (typical: 0.5-1.0)

**Parameter**: `sigma_weight` (default: 5.0)
**Description**: Weight penalty for small clusters
**Purpose**: Prevents creation of tiny, uninformative clusters
**Range**: 1.0-10.0 (typical: 3.0-7.0)

#### 5.3.4 Sequence Weighting

**Parameter**: `sequence_weighting` (default: 0)
**Options**: 
- **0**: "1/ns" - Weight by inverse number of similar sequences
- **1**: "clustering" - Weight based on cluster assignment confidence
- **2**: "none" - Equal weight for all sequences

**Purpose**: Handle redundancy and improve motif quality
**Recommendation**: Use option 0 for most applications

#### 5.3.5 Background Models

**Parameter**: `background_model` (default: 1)
**Options**:
- **0**: "flat" - Uniform amino acid frequencies (5% each)
- **1**: "pre-calculated" - BLOSUM-based evolutionary frequencies
- **2**: "from data" - Calculate from input sequences

**Biological Significance**:
- **Flat**: Assumes no amino acid bias
- **BLOSUM**: Reflects natural protein composition
- **Data-derived**: Accounts for library-specific biases

#### 5.3.6 Outlier Handling

**Parameter**: `use_trash_cluster` (default: False)
**Description**: Create additional cluster for outlier sequences
**Purpose**: Prevent outliers from distorting main clusters

**Parameter**: `trash_threshold` (default: 0.0)
**Description**: Probability threshold for trash cluster assignment
**Range**: 0.0-0.5 (typical: 0.1-0.3 when enabled)

---

## 📊 Section 6: Consensus Validation Parameters

### 6.1 Validation Control

**Parameter**: `use_consensus_validation` (default: True)
**Description**: Enable comprehensive cluster validation
**Benefits**: Robust cluster number selection, quality assessment
**Cost**: Increased computational time

### 6.2 Consensus Clustering

**Parameter**: `n_consensus_iterations` (default: 50)
**Description**: Number of subsampling iterations for stability assessment
**Purpose**: Build consensus matrix through repeated subsampling
**Range**: 20-200 (typical: 30-100)

**Trade-offs**:
- **More iterations**: Better stability estimates, longer runtime
- **Fewer iterations**: Faster analysis, less precise estimates

**Parameter**: Sample fraction (fixed: 0.8)
**Description**: Fraction of sequences used in each iteration
**Purpose**: Create diversity in subsamples while maintaining structure

### 6.3 Bootstrap Validation

**Parameter**: `enable_bootstrap` (default: True)
**Description**: Additional stability assessment through bootstrap resampling
**Purpose**: Complement consensus clustering with different validation approach

**Parameter**: `n_bootstrap_iterations` (default: 30)
**Description**: Number of bootstrap samples for stability testing
**Range**: 10-100 (typical: 20-50)

### 6.4 Automatic K Selection

**Parameter**: `auto_k_selection` (default: True)
**Description**: Automatically select optimal cluster number
**Method**: Composite scoring of stability, silhouette, and bootstrap metrics
**Weights**: 40% stability, 40% silhouette, 20% bootstrap (default)

**Quality Thresholds**:
- **Excellent**: Stability >0.8, Silhouette >0.7
- **Good**: Stability >0.6, Silhouette >0.5  
- **Moderate**: Stability >0.4, Silhouette >0.25
- **Poor**: Below moderate thresholds

---

## 📈 Section 7: Statistical Parameters

### 7.1 Differential Expression Thresholds

**Parameter**: `significance_threshold` (default: 0.05)
**Description**: FDR (adjusted p-value) cutoff for statistical significance
**Range**: 0.001-0.2 (typical: 0.01-0.1)
**Interpretation**: Expected false discovery rate among significant results

**Parameter**: `log2fc_threshold` (default: 1.5)
**Description**: Minimum absolute log2 fold change for biological significance
**Range**: 0.5-3.0 (typical: 1.0-2.0)
**Conversion**: log2FC=1.5 equals 2.83-fold change

**Combined Criteria**: Both statistical AND biological significance required
- **Stringent**: padj <0.01, |log2FC| >2.0
- **Standard**: padj <0.05, |log2FC| >1.5
- **Permissive**: padj <0.1, |log2FC| >1.0

### 7.2 Multiple Testing Considerations

**Method**: Benjamini-Hochberg FDR correction (automatic)
**Purpose**: Control expected proportion of false discoveries
**Implementation**: Built into DESeq2 workflow

**Interpretation**:
- **FDR = 0.05**: Expect 5% of significant results to be false positives
- **Conservative approach**: Use lower FDR (0.01-0.02)
- **Discovery approach**: Use higher FDR (0.1-0.2)

---

## 🎨 Section 8: Visualization Parameters

### 8.1 Plot Generation

**Parameter**: `plot_filtering` (default: True)
**Description**: Generate CPM filtering comparison plots
**Outputs**: Before/after histograms, correlation plots
**Purpose**: Quality control and method validation

### 8.2 Heatmap Options

**Scope**: Automatically generated for significant upregulated peptides
**Method**: VST (Variance Stabilizing Transformation) normalized values
**Color scheme**: Viridis (perceptually uniform, colorblind-friendly)
**Clustering**: Hierarchical clustering of both peptides and samples

### 8.3 Volcano Plot Customization

**Automatic generation**: Based on differential expression results
**Color coding**:
- **Red**: Significantly upregulated
- **Blue**: Significantly downregulated  
- **Grey**: Not significant
**Annotation**: Most significant peptides labeled automatically

### 8.4 Sequence Logos

**Generation**: Automatic for each discovered cluster
**Information content**: Height represents conservation strength
**Color scheme**: Standard biochemical amino acid colors
**Format**: Publication-ready vector graphics

---

## 💾 Section 9: Output Options

### 9.1 Result Files

**Clustering results**: CSV with cluster assignments and statistics
**Differential expression**: CSV with fold changes and p-values
**Validation metrics**: Comprehensive quality assessment report
**Visualizations**: High-resolution plots in PNG/PDF format

### 9.2 Quality Reports

**Validation summary**: Stability scores, silhouette analysis
**Method documentation**: Parameters used, interpretation guide
**Troubleshooting**: Common issues and recommended solutions

---

## 🔧 Section 10: Troubleshooting Guide

### 10.1 Common Parameter Issues

**Low success rate** (<30%):
- Check regex pattern matches actual sequences
- Verify peptide column format and content
- Consider sequence quality and length distribution

**No significant peptides**:
- Lower statistical thresholds (higher FDR, lower fold change)
- Check experimental design and sample labeling
- Verify sufficient library diversity and coverage

**Poor clustering quality**:
- Increase consensus validation iterations
- Try different cluster number ranges
- Check sequence length uniformity
- Consider amino acid property-based filtering

### 10.2 Performance Optimization

**Large datasets** (>1000 peptides):
- Reduce consensus iterations for initial exploration
- Use basic clustering for rapid analysis
- Consider subsetting data for parameter optimization

**Small datasets** (<50 peptides):
- Lower cluster numbers (2-4)
- Increase regularization (higher sigma_weight)
- Use shorter motif lengths if appropriate

### 10.3 Result Interpretation

**Validation scores**:
- **High stability, low silhouette**: May indicate overlapping clusters
- **Low stability, high silhouette**: May indicate insufficient sampling
- **Both low**: Consider different cluster numbers or parameters

**Biological validation**:
- Compare motifs to known binding sites
- Check amino acid property enrichment
- Validate key findings experimentally

---

## 📚 Best Practices Summary

### Parameter Selection Strategy:
1. **Start with defaults** for initial exploration
2. **Examine validation metrics** to assess quality
3. **Adjust based on results** and biological knowledge
4. **Document all parameters** for reproducibility

### Quality Assurance:
1. **Monitor success rates** and data quality metrics
2. **Use consensus validation** for robust results
3. **Cross-validate** findings with biological knowledge
4. **Report uncertainty** and confidence intervals

### Computational Efficiency:
1. **Use appropriate precision** for the research question
2. **Balance quality vs. speed** based on dataset size
3. **Parallelize when possible** for large analyses
4. **Save intermediate results** for parameter exploration

This comprehensive guide ensures informed parameter selection and optimal results from the peptide analysis pipeline.