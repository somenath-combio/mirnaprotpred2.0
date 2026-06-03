import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, precision_score, recall_score, f1_score

df = pd.read_csv("virmirna_external/virmirna_85_with_targets.csv")

def targetscan_score_user(mirna, cts):
    comp = {'A':'T','T':'A','G':'C','C':'G','U':'A'}
    mirna = mirna.upper().replace('U','T')
    cts   = cts.upper().replace('U','T')
    if len(mirna) < 8: return 0

    seed2_8 = mirna[1:8]
    seed2_7 = mirna[1:7]
    rc = lambda s: ''.join(comp.get(b,'N') for b in reversed(s))

    rc_2_8 = rc(seed2_8)
    rc_2_7 = rc(seed2_7)

    # User's logic: checking cts[idx-1] == 'A'
    if rc_2_8 in cts:
        idx = cts.index(rc_2_8)
        if idx > 0 and cts[idx-1] == 'A':
            return 4
        return 3

    if rc_2_7 in cts:
        idx = cts.index(rc_2_7)
        if idx > 0 and cts[idx-1] == 'A':
            return 2
        return 1

    return 0

def targetscan_score_biological(mirna, cts):
    comp = {'A':'T','T':'A','G':'C','C':'G','U':'A'}
    mirna = mirna.upper().replace('U','T')
    cts   = cts.upper().replace('U','T')
    if len(mirna) < 8: return 0

    seed2_8 = mirna[1:8]
    seed2_7 = mirna[1:7]
    rc = lambda s: ''.join(comp.get(b,'N') for b in reversed(s))

    rc_2_8 = rc(seed2_8)
    rc_2_7 = rc(seed2_7)

    # Biological TargetScan logic: opposite miRNA position 1 is downstream in 5'-to-3' target sequence
    if rc_2_8 in cts:
        idx = cts.index(rc_2_8)
        # Nucleotide downstream of the 7nt match
        next_idx = idx + 7
        if next_idx < len(cts) and cts[next_idx] == 'A':
            return 4
        return 3

    if rc_2_7 in cts:
        idx = cts.index(rc_2_7)
        # Nucleotide downstream of the 6nt match
        next_idx = idx + 6
        if next_idx < len(cts) and cts[next_idx] == 'A':
            return 2
        return 1

    return 0

# Run evaluation for both versions
for name, scoring_func in [("User's Indexing (cts[idx-1])", targetscan_score_user), 
                           ("Biological Indexing (cts[idx+len])", targetscan_score_biological)]:
    y_true, y_score = [], []
    for _, row in df.iterrows():
        t_seq = row.get('target_sequence', row.get('CTS_sequence', ''))
        score = scoring_func(str(row['miRNA_sequence']), str(t_seq))
        y_true.append(row['label'])
        y_score.append(score)

    y_true  = np.array(y_true)
    y_score = np.array(y_score)
    y_pred  = (y_score > 0).astype(int)

    auc  = roc_auc_score(y_true, y_score)
    prec = precision_score(y_true, y_pred, zero_division=0)
    rec  = recall_score(y_true, y_pred, zero_division=0)
    f1   = f1_score(y_true, y_pred, zero_division=0)

    print(f"=== TargetScan seed logic: {name} ===")
    print(f"AUC       : {auc:.4f}")
    print(f"Precision : {prec:.4f}")
    print(f"Recall    : {rec:.4f}")
    print(f"F1        : {f1:.4f}")
    print(f"Seed hits : {(y_score>0).sum()}/{len(y_true)}")
    print(f"  8mer    : {(y_score==4).sum()}")
    print(f"  7mer-m8 : {(y_score==3).sum()}")
    print(f"  7mer-A1 : {(y_score==2).sum()}")
    print(f"  6mer    : {(y_score==1).sum()}")
    print(f"  None    : {(y_score==0).sum()}")
    print("-" * 50)
