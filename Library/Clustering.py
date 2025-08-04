# ----------------------------------------------------------------------
# Pure-Python Gibbs clustering + logo generation
# ----------------------------------------------------------------------
from collections import defaultdict
import numpy as np, pandas as pd, random, math, matplotlib.pyplot as plt
import logomaker as lm                       # NEW

AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX      = {a:i for i,a in enumerate(AA_ALPHABET)}
IDX_AA      = {i:a for a,i in AA_IDX.items()}

def _encode(peps, L):
    enc = np.empty((len(peps), L), dtype=np.int8)
    for i,p in enumerate(peps):
        if len(p) != L:
            raise ValueError(f"Peptide {p!r} not length {L}")
        enc[i] = [AA_IDX[a] for a in p]
    return enc

def _counts_to_pwm(counts, alpha=1.0):
    """counts[k, pos, aa] → list[DataFrame] of probability PWMs"""
    # Handle both 2D and 3D input arrays
    if len(counts.shape) == 2:
        # 2D case: single cluster (pos, aa)
        mat = counts + alpha
        mat = mat / mat.sum(axis=1, keepdims=True)
        df = pd.DataFrame(mat, columns=list(AA_ALPHABET))
        return [df]
    elif len(counts.shape) == 3:
        # 3D case: multiple clusters (k, pos, aa)
        g, L, _ = counts.shape
        pwm_list = []
        for k in range(g):
            mat = counts[k] + alpha        # smooth
            mat = mat / mat.sum(axis=1, keepdims=True)  # normalize across amino acids (axis=1)
            df = pd.DataFrame(mat, columns=list(AA_ALPHABET))
            pwm_list.append(df)
        return pwm_list
    else:
        raise ValueError(f"Counts array must be 2D or 3D, got shape {counts.shape}")

def _logo_from_pwm(pwm_df, title="logo"):
    fig, ax = plt.subplots(figsize=(4,1.5))
    lm.Logo(pwm_df, ax=ax, shade_below=.5, fade_below=.5, stack_order='small_on_top')
    ax.set_xticks(range(pwm_df.shape[0]))
    ax.set_xlabel("Position")
    ax.set_ylabel("")
    ax.set_title(title, fontsize=10)
    plt.tight_layout()
    return fig

def gibbs_cluster(
    peptides,
    motif_length,
    num_clusters,
    out_dir=None, run_name=None,          # kept for API compat; unused
    *,
    n_iter=1000, burn_in=200,
    alpha=1.0, beta=1.0, seed=42,
    return_logos=True                     # NEW
):
    """
    Basic Gibbs clustering for peptide sequences
    
    Parameters:
    -----------
    peptides : list or pd.DataFrame
        Input peptide sequences
    motif_length : int
        Length of motifs to discover
    num_clusters : int
        Number of clusters to create
    n_iter : int
        Number of Gibbs sampling iterations
    alpha, beta : float
        Dirichlet hyperparameters
    seed : int
        Random seed for reproducibility
    return_logos : bool
        Whether to generate sequence logos
        
    Returns:
    --------
    tuple: (clusters_df, pwm_list, logo_figs)
        clusters_df: DataFrame with sequence assignments
        pwm_list: List of position weight matrices
        logo_figs: Dictionary of matplotlib figures
    """
    
    # Input validation
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    # Handle different input formats
    if isinstance(peptides, pd.DataFrame):
        if peptides.empty:
            raise ValueError("Input DataFrame is empty")
        # Try different column name options
        if 'variable_pep' in peptides.columns:
            seqs = list(peptides['variable_pep'])
        elif 'Sequence' in peptides.columns:
            seqs = list(peptides['Sequence'])
        else:
            seqs = list(peptides.iloc[:, 0])
    else:
        seqs = list(peptides)
    
    # Validate sequences
    if not seqs:
        raise ValueError("No sequences provided")
    
    # Remove empty or invalid sequences
    valid_seqs = []
    for seq in seqs:
        if isinstance(seq, str) and len(seq) >= motif_length:
            # Check for valid amino acids
            if all(aa in AA_ALPHABET for aa in seq.upper()):
                valid_seqs.append(seq.upper())
            else:
                print(f"Warning: Sequence contains invalid amino acids: {seq}")
        else:
            print(f"Warning: Sequence too short or invalid: {seq}")
    
    if not valid_seqs:
        raise ValueError("No valid sequences remain after filtering")
    
    seqs = valid_seqs
    N = len(seqs)
    L = motif_length
    
    # Validate parameters
    if num_clusters <= 0:
        raise ValueError("Number of clusters must be positive")
    if num_clusters > N:
        print(f"Warning: More clusters ({num_clusters}) than sequences ({N}). Reducing to {N}")
        num_clusters = N
    if motif_length <= 0:
        raise ValueError("Motif length must be positive")
    if alpha <= 0 or beta <= 0:
        raise ValueError("Alpha and beta must be positive")
    
    print(f"Clustering {N} sequences into {num_clusters} clusters (motif length: {L})")
    
    try:
        X = _encode(seqs, L)
    except Exception as e:
        raise ValueError(f"Error encoding sequences: {e}")

    # Initialize clustering
    z = np.random.randint(0, num_clusters, size=N)
    Nk = np.bincount(z, minlength=num_clusters)
    C = np.zeros((num_clusters, L, 20), dtype=int)
    
    for n in range(N):
        k = z[n]
        for j in range(L):
            C[k, j, X[n, j]] += 1

    # Gibbs sampling iterations
    print("Running Gibbs sampling...")
    for it in range(n_iter):
        if it % 100 == 0:
            print(f"  Iteration {it}/{n_iter}")
            
        for n in range(N):
            k_old = z[n]
            Nk[k_old] -= 1
            for j in range(L):
                C[k_old, j, X[n, j]] -= 1

            logp = np.empty(num_clusters)
            for k in range(num_clusters):
                # Avoid log(0) by ensuring beta > 0
                lp = math.log(max(Nk[k] + beta, 1e-10))
                denom = max(Nk[k] + 20*alpha, 1e-10)
                
                for j in range(L):
                    a = X[n, j]
                    numerator = max(C[k, j, a] + alpha, 1e-10)
                    lp += math.log(numerator / denom)
                logp[k] = lp
            
            # Numerical stability
            logp -= logp.max()
            p = np.exp(logp)
            p_sum = p.sum()
            if p_sum > 0:
                p /= p_sum
            else:
                p = np.ones(num_clusters) / num_clusters
            
            k_new = np.random.choice(num_clusters, p=p)

            z[n] = k_new
            Nk[k_new] += 1
            for j in range(L):
                C[k_new, j, X[n, j]] += 1

    # Prepare results
    clusters_df = pd.DataFrame({"Sequence": seqs, "Gn": z + 1})
    
    try:
        pwm_list = _counts_to_pwm(C, alpha)
    except Exception as e:
        print(f"Warning: Error creating PWMs: {e}")
        pwm_list = []
    
    logo_figs = {}
    if return_logos and pwm_list:
        try:
            for k, pwm in enumerate(pwm_list, start=1):
                logo_figs[k] = _logo_from_pwm(pwm, title=f"Cluster {k}")
        except Exception as e:
            print(f"Warning: Error creating logos: {e}")

    # Print clustering summary
    cluster_sizes = np.bincount(z, minlength=num_clusters)
    print("Clustering complete!")
    print("Cluster sizes:", cluster_sizes.tolist())
    
    return clusters_df, pwm_list, logo_figs
