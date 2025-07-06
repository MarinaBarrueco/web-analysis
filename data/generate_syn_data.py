import pandas as pd
import numpy as np
import string

# Set random seed for reproducibility
np.random.seed(42)

# Generate synthetic dataset
num_rows = 1000
amino_acids = list('ACDEFGHIKLMNPQRSTVWY')  # Standard amino acids

data = []
for i in range(num_rows):
    # Create a peptide matching the pattern A.C.{7}C
    seq = f"A{np.random.choice(amino_acids)}C{''.join(np.random.choice(amino_acids, 7))}C"
    # Simulate counts for 3 experimental and 3 control samples
    exp_counts = np.random.poisson(lam=100, size=3)   # higher mean for experiment
    ctrl_counts = np.random.poisson(lam=20, size=3)   # lower mean for control
    data.append([seq, *exp_counts, *ctrl_counts])

# Create DataFrame
columns = [
    "peptide",
    "Experiment 1", "Experiment 2", "Experiment 3",
    "Control 1", "Control 2", "Control 3"
]
df_large = pd.DataFrame(data, columns=columns)

# Save to CSV
df_large.to_csv("data/synthetic_peptide_data.csv", index=False)
print("Synthetic dataset generated and saved to 'data/synthetic_peptide_data.csv'.")

