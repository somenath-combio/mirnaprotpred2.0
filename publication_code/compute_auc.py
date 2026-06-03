import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, precision_score, recall_score, f1_score

# Load ground truth
df = pd.read_csv("virmirna_external/virmirna_85_with_targets.csv")
labels = {f"{row['miRNA']}|{row['label']}|{idx}": row['label']
          for idx, row in df.iterrows()}

# Parse miRanda hits
miranda_scores = {}
with open("benchmark/miranda_output_85.txt") as f:
    for line in f:
        if line.startswith(">>"):
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                mirna_header = parts[0].replace(">>","").strip()
                target_header = parts[1].strip()
                
                mirna_parts = mirna_header.split("|")
                target_parts = target_header.split("|")
                
                if len(mirna_parts) == 3 and len(target_parts) == 2:
                    mirna_idx = int(mirna_parts[2])
                    target_idx = int(target_parts[0].replace("target_", ""))
                    
                    if mirna_idx == target_idx:
                        score = float(parts[2])
                        max_score = float(parts[4])
                        energy = float(parts[3])
                        
                        miranda_scores[mirna_header] = (score, max_score, energy)

# Build prediction vectors
y_true = []
y_score_align = []
y_score_max = []
y_score_energy = []

for key, label in labels.items():
    y_true.append(label)
    if key in miranda_scores:
        y_score_align.append(miranda_scores[key][0])
        y_score_max.append(miranda_scores[key][1])
        y_score_energy.append(-miranda_scores[key][2])
    else:
        y_score_align.append(0.0)
        y_score_max.append(0.0)
        y_score_energy.append(0.0)

y_true = np.array(y_true)
y_score_align = np.array(y_score_align)
y_score_max = np.array(y_score_max)
y_score_energy = np.array(y_score_energy)

# Metrics for Alignment Score
y_pred_align = (y_score_align > 0).astype(int)
auc_align = roc_auc_score(y_true, y_score_align)
prec_align = precision_score(y_true, y_pred_align, zero_division=0)
rec_align = recall_score(y_true, y_pred_align, zero_division=0)
f1_align = f1_score(y_true, y_pred_align, zero_division=0)

# Metrics for Max Score
y_pred_max = (y_score_max > 0).astype(int)
auc_max = roc_auc_score(y_true, y_score_max)
prec_max = precision_score(y_true, y_pred_max, zero_division=0)
rec_max = recall_score(y_true, y_pred_max, zero_division=0)
f1_max = f1_score(y_true, y_pred_max, zero_division=0)

# Metrics for Energy
auc_energy = roc_auc_score(y_true, y_score_energy)

print("=== miRanda on 85-pair test set ===")
print("Using Alignment Score (Column 2):")
print(f"  AUC       : {auc_align:.4f}")
print(f"  Precision : {prec_align:.4f}")
print(f"  Recall    : {rec_align:.4f}")
print(f"  F1        : {f1_align:.4f}")
print(f"  Hits found: {(y_score_align>0).sum()}/{len(y_true)}")
print()
print("Using Max Score (Column 4, parsed in original script):")
print(f"  AUC       : {auc_max:.4f}")
print(f"  Precision : {prec_max:.4f}")
print(f"  Recall    : {rec_max:.4f}")
print(f"  F1        : {f1_max:.4f}")
print(f"  Hits found: {(y_score_max>0).sum()}/{len(y_true)}")
print()
print("Using Thermodynamic Energy (-MFE, Column 3):")
print(f"  AUC       : {auc_energy:.4f}")
