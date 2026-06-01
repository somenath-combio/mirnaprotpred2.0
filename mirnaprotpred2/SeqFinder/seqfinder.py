#!/usr/bin/env python3
"""
seqfinder.py — Stage 1 of miRNAProtPred 2.0 pipeline
───────────────────────────────────────────────────────
Thermodynamic scan of a viral genome for all physically possible
miRNA binding sites. Uses seed matching + IntaRNA delta_G.

Usage:
  python3 seqfinder.py \
    --mirna   UAGCUUAUCAGACUGAUGUUGA \
    --genome  hcv_genome.fasta \
    --top     50 \
    --out     candidates.csv

Output CSV columns:
  start, end, window_seq, delta_G, hybridDP,
  seed_match, supp_match, cts_gc, n_base_pairs
"""

import argparse
import subprocess
import csv
import sys
from pathlib import Path


def assign_confidence(prob, dg):
    """
    3-tier confidence for novel virus generalization.
    No retraining needed for unseen virus families.

    High   : prob >= 0.50              (matches VIRBase experimental profile)
    Medium : prob >= 0.25 AND dG<=-9.0 (thermodynamically supported)
    Low    : prob <  0.25 OR  dG > -9.0 (computational prediction only)

    Trained on: HCV, HBV, HIV, DENV, Influenza, EBV (VIRBase).
    Zero-shot validated: SARS-CoV-2 3-UTR, hsa-miR-3941 (PMID:34198800).
    """
    if prob >= 0.50:
        return "High"
    if prob >= 0.25 and dg <= -9.0:
        return "Medium"
    return "Low"

INTARNA_BIN = "/home/somenath/miniconda3/envs/cts_env/bin/IntaRNA"

# ── helpers ────────────────────────────────────────────────────────

def rc(seq):
    comp = {'A':'T','T':'A','G':'C','C':'G','U':'A','N':'N'}
    return ''.join(comp.get(b,'N') for b in reversed(seq.upper()))

def seed_match(mirna, cts):
    """7-mer seed complement check (positions 2-8 of miRNA)."""
    if len(mirna) < 8: return 0
    seed = mirna[1:8].upper().replace('U','T')
    seed_rc = rc(seed)
    return 1 if seed_rc in cts.upper().replace('U','T') else 0

def supp_match(mirna, cts):
    """Supplementary pairing (positions 12-16)."""
    if len(mirna) < 16: return 0
    supp = mirna[11:16].upper().replace('U','T')
    supp_rc = rc(supp)
    return 1 if supp_rc in cts.upper().replace('U','T') else 0

def gc_content(seq):
    seq = seq.upper().replace('U','C')
    return sum(1 for b in seq if b in 'GC') / max(len(seq), 1)

def run_intarna(mirna, target, timeout=20):
    cmd = [
        INTARNA_BIN,
        "-q", mirna.upper().replace('T','U'),
        "-t", target.upper().replace('T','U'),
        "--outMode", "C",
        "--outCsvCols", "E,hybridDP,start1,start2,end1,end2",
        "--threads", "1",
        "--seedBP", "7",
    ]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        lines = [l for l in r.stdout.strip().split('\n') if l]
        if len(lines) < 2: return None
        p = lines[1].strip().split(';')
        if len(p) < 2: return None
        dg_str = p[0].strip()
        if not dg_str.lstrip('-').replace('.','').isdigit(): return None
        hyb = p[1].strip() if len(p) > 1 else ""
        nbp = hyb.split('&')[0].count('(') if '&' in hyb else 0
        return float(dg_str), hyb, nbp
    except Exception:
        return None

def load_fasta(path):
    """Returns concatenated sequence from a FASTA file."""
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'): continue
            seq.append(line.upper().replace('U','T'))
    return ''.join(seq)

# ── main scan ──────────────────────────────────────────────────────

def scan_genome(mirna_seq, genome_seq, win_len=50, step=10,
                dg_threshold=-5.0, top_n=50, seed_required=False):
    """
    Slide a window across the genome and collect all candidate
    sites where IntaRNA finds binding with delta_G <= dg_threshold.

    Parameters
    ----------
    mirna_seq      : str  — mature miRNA sequence (T or U)
    genome_seq     : str  — full viral mRNA/genome (T or U)
    win_len        : int  — window length to scan (default 50 nt)
    step           : int  — step size for sliding window (default 10 nt)
    dg_threshold   : float — maximum delta_G to retain (default -5.0)
    top_n          : int  — return only top N sites by delta_G
    seed_required  : bool — if True, only return sites with seed_match=1

    Returns
    -------
    List of dicts, sorted by delta_G ascending (most negative first)
    """
    candidates = []
    genome_len = len(genome_seq)
    total_windows = (genome_len - win_len) // step

    print(f"  Genome length  : {genome_len:,} nt")
    print(f"  Window length  : {win_len} nt")
    print(f"  Step size      : {step} nt")
    print(f"  Total windows  : {total_windows:,}")
    print(f"  delta_G cutoff : <= {dg_threshold}")
    print(f"  Seed required  : {seed_required}")
    print(f"  Scanning", end="", flush=True)

    scanned = 0
    for start in range(0, genome_len - win_len + 1, step):
        end = start + win_len
        window = genome_seq[start:end]
        scanned += 1
        if scanned % 500 == 0:
            print(".", end="", flush=True)

        # Fast pre-filter: seed match check (avoids IntaRNA call)
        sm = seed_match(mirna_seq, window)
        if seed_required and sm == 0:
            continue

        result = run_intarna(mirna_seq, window)
        if result is None: continue
        dg, hyb, nbp = result
        if dg > dg_threshold: continue

        candidates.append({
            'start'      : start,
            'end'        : end,
            'window_seq' : window,
            'delta_G'    : round(dg, 3),
            'hybridDP'   : hyb,
            'n_base_pairs': nbp,
            'seed_match' : sm,
            'supp_match' : supp_match(mirna_seq, window),
            'cts_gc'     : round(gc_content(window), 4),
        })

    print(f"\n  Scanned: {scanned:,} windows | Candidates: {len(candidates)}")

    # Sort by delta_G ascending (most stable binding first)
    candidates.sort(key=lambda x: x['delta_G'])
    return candidates[:top_n]


def stage2_score(candidates, model_pkl, train_csv):
    """
    Run Stage 2: score candidates through miRNAProtPred 2.0.
    Requires: the saved RandomForest model and training CSV for feature column names.
    """
    import pickle, pandas as pd, numpy as np

    with open(model_pkl, 'rb') as f:
        model = pickle.load(f)

    ref = pd.read_csv(train_csv, nrows=0)
    meta = ['label','miRNA','virus','target_symbol']
    DROP = ['delta_G','delta_G_norm','n_base_pairs','struct_matches',
            'struct_wobbles','struct_mismatches','struct_bulges',
            'struct_loops','site_position_norm']
    feat_cols = [c for c in ref.columns if c not in meta and c not in DROP]

    def dinuc(seq):
        bases = 'ACGT'; seq = seq.upper().replace('U','T')
        d = {a+b: 0 for a in bases for b in bases}
        total = max(len(seq)-1, 1)
        for i in range(len(seq)-1):
            k = seq[i:i+2]
            if k in d: d[k] += 1
        return {f'dn_{k}': v/total for k,v in d.items()}

    def trinuc(seq):
        seq = seq.upper().replace('U','T')
        cats = {'AAA':0,'AAG':0,'GAA':0,'GAG':0,
                'CCC':0,'CCT':0,'TCC':0,'TCT':0}
        total = max(len(seq)-2, 1)
        for i in range(len(seq)-2):
            t = seq[i:i+3]
            if t in cats: cats[t] += 1
        return {f'tn_{k}': v/total for k,v in cats.items()}

    rows = []
    for c in candidates:
        row = {}
        row['cts_gc'] = c['cts_gc']
        row['cts_au'] = sum(1 for b in c['window_seq'].upper() if b in 'ATU') / max(len(c['window_seq']),1)
        row['mirna_gc'] = 0.0   # will be filled below
        row['seed_match'] = c['seed_match']
        row['supp_match'] = c['supp_match']
        row['motif_identity'] = c['seed_match'] * c['supp_match']
        row.update(dinuc(c['window_seq']))
        row.update(trinuc(c['window_seq']))
        rows.append(row)

    df = pd.DataFrame(rows)
    for col in feat_cols:
        if col not in df.columns: df[col] = 0
    df = df[feat_cols].fillna(0)

    probs = model.predict_proba(df)[:,1]
    for i, c in enumerate(candidates):
        c['ml_score'] = round(float(probs[i]), 4)

    candidates.sort(key=lambda x: -x['ml_score'])
    return candidates


# ── CLI ────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='SeqFinder + miRNAProtPred 2.0: two-stage viral miRNA target scanner'
    )
    parser.add_argument('--mirna',   required=True,
                        help='Mature miRNA sequence (A/T/U/G/C)')
    parser.add_argument('--genome',  required=True,
                        help='Path to viral mRNA/genome FASTA')
    parser.add_argument('--win',     type=int, default=50,
                        help='Sliding window length (default 50 nt)')
    parser.add_argument('--step',    type=int, default=10,
                        help='Step size (default 10 nt)')
    parser.add_argument('--dg',      type=float, default=-5.0,
                        help='delta_G threshold, e.g. -5.0 (default)')
    parser.add_argument('--top',     type=int, default=50,
                        help='Return top N candidates (default 50)')
    parser.add_argument('--seed',    action='store_true',
                        help='Only return windows with 7-mer seed match')
    parser.add_argument('--out',     default='seqfinder_output.csv',
                        help='Output CSV path')
    parser.add_argument('--model',   default=None,
                        help='Path to mirnaprotpred2_best.pkl for Stage 2 scoring')
    parser.add_argument('--traincsv',default=None,
                        help='Path to training CSV (for feature column alignment)')
    args = parser.parse_args()

    print("="*60)
    print("  SeqFinder — Stage 1: Thermodynamic Genome Scan")
    print("="*60)
    print(f"  miRNA   : {args.mirna}")
    print(f"  Genome  : {args.genome}")

    genome_seq = load_fasta(args.genome)
    candidates = scan_genome(
        mirna_seq      = args.mirna,
        genome_seq     = genome_seq,
        win_len        = args.win,
        step           = args.step,
        dg_threshold   = args.dg,
        top_n          = args.top,
        seed_required  = args.seed,
    )

    if not candidates:
        print("No candidates found. Try relaxing --dg or --win.")
        sys.exit(0)

    if args.model and args.traincsv:
        print("\n" + "="*60)
        print("  miRNAProtPred 2.0 — Stage 2: ML Contextual Scoring")
        print("="*60)
        candidates = stage2_score(candidates, args.model, args.traincsv)
        print(f"  Stage 2 complete. ML scores added.")

    # Write output
    out_path = Path(args.out)
    fieldnames = list(candidates[0].keys())
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(candidates)

    print(f"\n  Results saved → {out_path}")
    print(f"  Top 5 candidates:")
    for i, c in enumerate(candidates[:5]):
        score_str = f"  ml_score={c['ml_score']:.3f}" if 'ml_score' in c else ""
        print(f"    [{i+1}] pos={c['start']}-{c['end']}  "
              f"dG={c['delta_G']}  seed={c['seed_match']}"
              f"  gc={c['cts_gc']:.2f}{score_str}")

    print("="*60)

cli = main

if __name__ == '__main__':
    cli()
