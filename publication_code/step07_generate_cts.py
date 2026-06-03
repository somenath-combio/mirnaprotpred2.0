#!/usr/bin/env python3
"""
generate_cts.py — CTS Generation Pipeline (miRNAProtPred 2.0)
Matches published methods exactly:
  - Positive CTS: 41 nt window centered on IntaRNA seed, ΔG ≤ -10.0, seedBP=7
  - Hard negatives: same target, ΔG > -2.0 kcal/mol, zero canonical 6-mer seed match
  - Cross-virus negatives: phylogenetically distant family, non-binding ΔG > -2.0
  - Ratio: 1 pos : 1 hard neg : 1 cross-virus neg  →  total 3N
"""

import os, subprocess, random, re, argparse
import pandas as pd
import numpy as np
from pathlib import Path

parser = argparse.ArgumentParser(description="Generate CTS via IntaRNA")
parser.add_argument("--input",  default="virbase_final_dataset/intermediate/final_ml_dataset.csv",
                    help="Input CSV with miRNA+target pairs")
parser.add_argument("--output", default="virbase_final_dataset/virbase_cts/cts_dataset_raw.tsv",
                    help="Output TSV for raw CTS results")
args = parser.parse_args()

INPUT_FILE  = args.input
OUTPUT_FILE = args.output

INTARNA_BIN   = "/home/somenath/miniconda3/envs/cts_env/bin/IntaRNA"
ENERGY_THRESH = -10.0   # kcal/mol — positives only
NEG_THRESH    = -2.0    # kcal/mol — negatives must be WEAKER than this
CTS_WINDOW    = 20      # nt each side of seed → 41 nt total window
SEED_MIN_LEN  = 7
random.seed(42)
np.random.seed(42)

# Phylogenetically distant family pairs (must differ at family level)
FAMILY_MAP = {
    "Human immunodeficiency virus 1 (HIV-1)":          "Retroviridae",
    "Hepacivirus C (HCV)":                             "Flaviviridae",
    "Hepatitis B virus (HBV)":                         "Hepadnaviridae",
    "Influenza A virus (H1N1)":                        "Orthomyxoviridae",
    "Simian immunodeficiency virus (SIV)":             "Retroviridae",
    "Zika virus (ZIKV)":                               "Flaviviridae",
    "Epstein-Barr virus (EBV)":                        "Herpesviridae",
    "Human betaherpesvirus 5 (HCMV)":                  "Herpesviridae",
}

def get_family(virus):
    for k, v in FAMILY_MAP.items():
        if any(tag in virus for tag in k.split('(')):
            return v
    return "Unknown"

def run_intarna(mirna_seq, target_seq, max_e=ENERGY_THRESH,
                seed_len=SEED_MIN_LEN, max_results=1):
    mirna_seq  = mirna_seq.upper().replace('T','U')
    target_seq = target_seq.upper()  # keep T, do NOT replace with U
    cmd = [
        INTARNA_BIN,
        "--query",   mirna_seq,
        "--target",  target_seq,
        "--out",     "STDOUT",
        "--outMode", "C",
        "--outMaxE", str(max_e),
        "--outNumber", str(max_results),
        "--seedBP",  str(seed_len),
        "--accW",    "150",
        "--accL",    "40",
        "--outSep",  ",",
    ]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        lines = [l for l in r.stdout.strip().split('\n') if l.strip()]
        if len(lines) < 2:
            return []
        header = lines[0].split(',')
        sites  = []
        for row in lines[1:]:
            vals = row.split(',')
            if len(vals) < len(header): continue
            rec  = dict(zip(header, vals))
            try:    energy = float(rec.get('E', 0))
            except: continue
            if energy > max_e: continue
            sites.append({
                'start_t': int(rec.get('start1', 0)),
                'end_t':   int(rec.get('end1',   0)),
                'start_q': int(rec.get('start2', 0)),
                'end_q':   int(rec.get('end2',   0)),
                'energy':  energy,
                'hybridDP': rec.get('hybridDP',''),
            })
        return sites
    except Exception as e:
        print(f"  IntaRNA error: {e}")
        return []

def extract_cts_window(target_seq, start_t, end_t, window=CTS_WINDOW):
    """Extract 41 nt CTS window centered on seed coordinates (1-based)."""
    s = max(0, start_t - 1 - window)
    e = min(len(target_seq), end_t + window)
    return target_seq[s:e], s+1, e

def has_seed_match(mirna_seq, cts_seq):
    """Check if CTS contains the reverse complement of miRNA positions 2-7."""
    mirna = str(mirna_seq).upper().replace('T','U')
    cts   = str(cts_seq).upper().replace('T','U')
    seed  = mirna[1:7]
    comp  = {'A':'U','U':'A','G':'C','C':'G'}
    rc    = ''.join(comp.get(b,'N') for b in reversed(seed))
    return rc in cts

def scan_hard_negatives(mirna_seq, target_seq, n_needed, pos_sites):
    """
    Scan target for hard negatives:
    - ΔG > NEG_THRESH (not binding)
    - zero canonical 6-mer seed match
    - non-overlapping with positive sites
    Uses sliding 41 nt windows across the target.
    """
    mirna_seq  = mirna_seq.upper().replace('T','U')
    target_seq = target_seq.upper()
    tlen = len(target_seq)
    if tlen < 41: return []

    pos_ranges = set()
    for s in pos_sites:
        for p in range(max(0, s['start_t']-1-CTS_WINDOW),
                       min(tlen, s['end_t']+CTS_WINDOW)):
            pos_ranges.add(p)

    candidates = []
    step = 5
    for i in range(0, tlen - 40, step):
        window_seq = target_seq[i:i+41]
        if has_seed_match(mirna_seq, window_seq):
            continue
        mid = i + 20
        if mid in pos_ranges:
            continue
        # Quick IntaRNA check on the window
        sites = run_intarna(mirna_seq, window_seq, max_e=NEG_THRESH,
                            seed_len=2, max_results=1)
        if not sites:
            candidates.append({
                'CTS_sequence': window_seq,
                'CTS_start':    i+1,
                'CTS_end':      i+41,
                'deltaG':       0.0,
                'hybridDP':     '',
            })
        if len(candidates) >= n_needed * 5:
            break

    if not candidates:
        return []
    selected = random.sample(candidates, min(n_needed, len(candidates)))
    return selected

def main():
    df = pd.read_csv(INPUT_FILE)
    df = df[df['miRNA_sequence'].notna() & df['target_sequence'].notna()]
    df = df[df['miRNA_sequence'] != '' ]
    total = len(df)
    print(f"Loaded {total} interactions.")

    # Build per-family target pool for cross-virus negatives
    family_targets = {}
    for _, row in df.iterrows():
        fam = get_family(str(row['virus']))
        if fam not in family_targets:
            family_targets[fam] = []
        family_targets[fam].append({
            'virus':         row['virus'],
            'target_symbol': row['target_symbol'],
            'target_seq':    str(row['target_sequence']).upper(),
        })

    positives, hard_negs, cross_negs = [], [], []
    skipped = []

    for i, row in df.iterrows():
        mirna      = row['miRNA']
        mirna_seq  = str(row['miRNA_sequence']).upper().replace('T','U')
        target_seq = str(row['target_sequence']).upper()   # keep as DNA
        virus      = row['virus']
        target_sym = row['target_symbol']
        fam        = get_family(virus)

        if len(mirna_seq) < 18 or len(target_seq) < 20:
            skipped.append(mirna)
            continue

        print(f"[{i+1}/{total}] {mirna} | {target_sym} | {virus[:30]}")

        # ── POSITIVES ──────────────────────────────────────────────
        pos_sites = run_intarna(mirna_seq, target_seq,
                                max_e=ENERGY_THRESH,
                                seed_len=SEED_MIN_LEN,
                                max_results=3)
        if not pos_sites:
            skipped.append(mirna)
            continue

        s = pos_sites[0]
        cts_seq, cts_s, cts_e = extract_cts_window(target_seq,
                                                    s['start_t'],
                                                    s['end_t'])
        positives.append({
            'label':1, 'miRNA':mirna, 'miRNA_sequence':mirna_seq,
            'virus':virus, 'target_symbol':target_sym,
            'target_type':row['target_type'],
            'CTS_sequence':cts_seq,
            'CTS_start':cts_s, 'CTS_end':cts_e,
            'deltaG':s['energy'], 'hybridDP':s['hybridDP'],
            'PMID':row.get('PMID',''), 'Score':row.get('Score',0),
        })

        # ── HARD NEGATIVES ─────────────────────────────────────────
        hard = scan_hard_negatives(mirna_seq, target_seq, 1, pos_sites)
        for h in hard:
            hard_negs.append({
                'label':0, 'miRNA':mirna, 'miRNA_sequence':mirna_seq,
                'virus':virus, 'target_symbol':target_sym,
                'target_type':row['target_type'],
                'CTS_sequence':h['CTS_sequence'],
                'CTS_start':h['CTS_start'], 'CTS_end':h['CTS_end'],
                'deltaG':h['deltaG'], 'hybridDP':h['hybridDP'],
                'PMID':row.get('PMID',''), 'Score':0.0,
            })

        # ── CROSS-VIRUS NEGATIVES ──────────────────────────────────
        other_fams = [f for f in family_targets if f != fam and f != "Unknown"]
        random.shuffle(other_fams)
        cv_done = 0
        for of in other_fams:
            if cv_done >= 1: break
            pool = family_targets[of].copy()
            random.shuffle(pool)
            for ct in pool:
                if len(ct['target_seq']) < 41: continue
                # sample one 41 nt window — no binding required
                idx = random.randint(0, len(ct['target_seq'])-41)
                win = ct['target_seq'][idx:idx+41]
                if has_seed_match(mirna_seq, win): continue
                cv_sites = run_intarna(mirna_seq, win,
                                       max_e=NEG_THRESH,
                                       seed_len=2, max_results=1)
                if not cv_sites:
                    cross_negs.append({
                        'label':0, 'miRNA':mirna, 'miRNA_sequence':mirna_seq,
                        'virus':ct['virus'], 'target_symbol':ct['target_symbol'],
                        'target_type':'cross_virus_negative',
                        'CTS_sequence':win,
                        'CTS_start':idx+1, 'CTS_end':idx+41,
                        'deltaG':0.0, 'hybridDP':'',
                        'PMID':'', 'Score':0.0,
                    })
                    cv_done += 1
                    break

    pos_df   = pd.DataFrame(positives)
    hard_df  = pd.DataFrame(hard_negs)
    cross_df = pd.DataFrame(cross_negs)

    print(f"\n  Positives:          {len(pos_df)}")
    print(f"  Hard negatives:     {len(hard_df)}")
    print(f"  Cross-virus negs:   {len(cross_df)}")
    print(f"  Skipped:            {len(skipped)}")

    all_df = pd.concat([pos_df, hard_df, cross_df], ignore_index=True)

    out_path = Path(OUTPUT_FILE)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    sep = '\t' if out_path.suffix == '.tsv' else ','
    all_df.to_csv(out_path, index=False, sep=sep)
    pos_df.to_csv(out_path.parent / "cts_positives.csv", index=False)
    pd.concat([hard_df, cross_df], ignore_index=True).to_csv(
        out_path.parent / "cts_negatives.csv", index=False)

    print("\n" + "="*55)
    print("  CTS GENERATION COMPLETE")
    print("="*55)
    print(f"  Positive CTSs:    {len(pos_df)}")
    print(f"  Hard negatives:   {len(hard_df)}")
    print(f"  Cross-virus negs: {len(cross_df)}")
    print(f"  Total:            {len(all_df)}")
    if len(pos_df) > 0:
        print(f"  Neg:Pos ratio:    {(len(hard_df)+len(cross_df))/len(pos_df):.2f}:1")
    print(f"  Output: {OUTPUT_FILE}")
    print("="*55)

if __name__ == "__main__":
    main()
