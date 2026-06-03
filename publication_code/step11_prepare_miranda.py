import pandas as pd
import numpy as np
import random
import os

random.seed(42)
np.random.seed(42)

# ALWAYS read from normalized file to avoid duplicate appending bug
csv_path = "virmirna_external/virmirna_85_normalized.csv"

df_pos = pd.read_csv(csv_path)
df_pos['label'] = 1

# Generate hard negatives by shuffling target sequence
neg_rows = []
for idx, row in df_pos.iterrows():
    t_seq = row.get('target_sequence', row.get('CTS_sequence', ''))
    cts_seq = str(t_seq).upper().replace('T', 'U')
    
    seq_list = list(cts_seq)
    random.shuffle(seq_list)
    neg_seq = "".join(seq_list).replace('U', 'T')
    
    neg_row = row.copy()
    if 'target_sequence' in neg_row:
        neg_row['target_sequence'] = neg_seq
    if 'CTS_sequence' in neg_row:
        neg_row['CTS_sequence'] = neg_seq
    neg_row['label'] = 0
    neg_rows.append(neg_row)

df_neg = pd.DataFrame(neg_rows)
df_all = pd.concat([df_pos, df_neg], ignore_index=True)

# Save the unified dataset with both positive and negative target sequences
df_all.to_csv("virmirna_external/virmirna_85_with_targets.csv", index=False)

# Write miRNA FASTA with unique indexes
with open("benchmark/mirnas_85.fa", "w") as f:
    for idx, row in df_all.iterrows():
        seq = str(row['miRNA_sequence']).upper().replace('U','T')
        f.write(f">{row['miRNA']}|{row['label']}|{idx}\n{seq}\n")

# Write Target FASTA with unique indexes
with open("benchmark/targets_85.fa", "w") as f:
    for idx, row in df_all.iterrows():
        t_seq = row.get('target_sequence', row.get('CTS_sequence', ''))
        seq = str(t_seq).upper().replace('U','T')
        f.write(f">target_{idx}|{row['label']}\n{seq}\n")

print("Created 170 pairs (85 positive, 85 negative) in virmirna_85_with_targets.csv with unique indexes")
print("FASTA counts:")
print("miRNAs:", len(df_all))
print("targets:", len(df_all))
