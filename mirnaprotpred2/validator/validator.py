#!/usr/bin/env python3
"""
Validator - Targeted miRNA Interaction Verification Engine for miRNAProtPred 2.0.

Validates SPECIFIC user-provided miRNAs against a target sequence.
Internally powered by the SeqFinder2 engine.
"""

import os
import sys
import time
import argparse
from pathlib import Path
import pandas as pd
from Bio import SeqIO

from mirnaprotpred2.SeqFinder.seqfinder import (
    scan_genome_multiple,
    stage2_score,
    load_target_sequence,
    assign_confidence
)

# Resolving default model paths
pkg_data_dir = Path(__file__).resolve().parent.parent / "SeqFinder" / "data"
default_model = pkg_data_dir / "mirnaprotpred2_best.pkl"
default_train = pkg_data_dir / "cts_ml_features_v4_intraviral.csv"
default_xlsx  = pkg_data_dir / "data.xlsx"

def parse_mirna_ids(mirna_input):
    """
    Parse miRNA IDs from comma-separated string, txt file, or FASTA file headers.
    """
    if os.path.isfile(mirna_input):
        ext = os.path.splitext(mirna_input)[1].lower()
        if ext in {".fasta", ".fa", ".fna"}:
            ids = []
            for record in SeqIO.parse(mirna_input, "fasta"):
                ids.append(record.id)
            return list(dict.fromkeys(ids))
        
        # Text/CSV file
        ids = []
        with open(mirna_input, "r") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    for part in line.split(","):
                        part = part.strip()
                        if part:
                            ids.append(part)
        return list(dict.fromkeys(ids))

    # Inline comma-separated
    return [mid.strip() for mid in mirna_input.split(",") if mid.strip()]

def main():
    # Handle backward compatibility and positional argument redirection
    args_list = sys.argv[1:]
    
    # Remove '--mirna' and '--genome' if present
    cleaned_args = []
    skip_next = False
    for arg in args_list:
        if skip_next:
            skip_next = False
            continue
        if arg in ['--mirna', '--genome']:
            skip_next = True
            continue
        cleaned_args.append(arg)

    parser = argparse.ArgumentParser(
        prog='validator2',
        description='validator2: Targeted miRNA Interaction Verification Engine for miRNAProtPred 2.0'
    )
    # The first two positional arguments: mirna and genome
    parser.add_argument('mirna',
                        help='miRNA IDs (comma-separated, .txt file, or .fasta file of miRNAs)')
    parser.add_argument('genome',
                        help='Path to target viral genome or mRNA FASTA file')
    parser.add_argument('--win', type=int, default=50,
                        help='Sliding target window length (default: 50)')
    parser.add_argument('--dg', type=float, default=-5.0,
                        help='delta_G threshold (default: -5.0)')
    parser.add_argument('--mode', choices=['strict', 'relaxed'], default='strict',
                        help='Search mode: strict enforces 7-mer seed matches, relaxed allows non-canonical matches (default: strict)')
    parser.add_argument('--output', choices=['summary', 'raw', 'highconf'], default='summary',
                        help='Output mode (default: summary)')
    parser.add_argument('--out', default=None,
                        help='Save validation results automatically to the specified CSV filename')
    parser.add_argument('--details', action='store_true',
                        help='Alias for --output raw')
    parser.add_argument('--model', default=str(default_model) if default_model.exists() else None,
                        help='Path to custom mirnaprotpred2_best.pkl for Stage 2 scoring')
    parser.add_argument('--traincsv', default=str(default_train) if default_train.exists() else None,
                        help='Path to training CSV for feature column alignment')
    parser.add_argument('--email', default='sudipta@pusan.ac.kr',
                        help='NCBI Entrez email for Protein mode blast searches')

    args = parser.parse_args(cleaned_args)

    if args.details:
        args.output = 'raw'

    print("="*60)
    print("  Validator2 — Targeted miRNA Verification (miRNAProtPred 2.0)")
    print("="*60)
    print(f"  miRNA Input : {args.mirna}")
    print(f"  Genome FASTA: {args.genome}")
    print(f"  Window Size : {args.win}")
    print(f"  delta_G Cut : <= {args.dg}")
    print(f"  Mode        : {args.mode}")
    print(f"  Output      : {args.output}")

    if args.email:
        from Bio import Entrez
        Entrez.email = args.email

    # 1. Parse miRNA input
    requested_ids = parse_mirna_ids(args.mirna)
    if not requested_ids:
        print("Error: No miRNA IDs found in input.", file=sys.stderr)
        sys.exit(1)

    print(f"  Parsed {len(requested_ids)} miRNA ID(s). Searching sequences...")

    # 2. Resolve miRNA sequences
    mirnas = []
    # If miRNA input is a FASTA file, we can read sequences directly from it
    if os.path.isfile(args.mirna) and os.path.splitext(args.mirna)[1].lower() in {".fasta", ".fa", ".fna"}:
        for record in SeqIO.parse(args.mirna, "fasta"):
            mirnas.append((record.id, str(record.seq).strip().upper()))
    else:
        # Load default database to look up sequences
        if not default_xlsx.exists():
            print(f"Error: Packed miRNA database not found at '{default_xlsx}'", file=sys.stderr)
            sys.exit(1)
        
        df_db = pd.read_excel(default_xlsx, engine='openpyxl')
        db_map = {}
        for _, row in df_db.iterrows():
            db_id = str(row.get('Human miRNA ID', '')).strip()
            db_seq = str(row.get('Sequence', '')).strip().upper()
            if db_id and db_seq:
                db_map[db_id.lower()] = (db_id, db_seq)

        for req_id in requested_ids:
            req_clean = req_id.lower().replace("-", "").strip()
            matched = False
            # Try case-insensitive substring matching on hyphen-free IDs
            for k, (db_id, db_seq) in db_map.items():
                k_clean = k.replace("-", "")
                if req_clean in k_clean or k_clean in req_clean:
                    mirnas.append((db_id, db_seq))
                    matched = True
            if not matched:
                print(f"  [Warning] Could not find mature sequence for miRNA: {req_id}")

    if not mirnas:
        print("Error: No valid miRNA sequences could be loaded.", file=sys.stderr)
        sys.exit(1)

    # 3. Load Genome Sequence
    genome_seq = load_target_sequence(args.genome)

    # 4. Scan using SeqFinder2 Multiple Engine
    print(f"\nScanning genome against {len(mirnas)} miRNA sequence(s)...")
    start_time = time.time()
    
    seed_required = (args.mode == 'strict')
    candidates = []
    
    if seed_required:
        candidates = scan_genome_multiple(
            mirnas         = mirnas,
            genome_seq     = genome_seq,
            win_len        = args.win,
            dg_threshold   = args.dg,
            top_n          = 999999, # retrieve all matches
        )
    else:
        # Relaxed mode: run scan_genome for each miRNA sequentially to support wobble/non-canonical scans
        from mirnaprotpred2.SeqFinder.seqfinder import scan_genome
        for m_id, m_seq in mirnas:
            res = scan_genome(
                mirna_id       = m_id,
                mirna_seq      = m_seq,
                genome_seq     = genome_seq,
                win_len        = args.win,
                step           = 1, # check all positions
                dg_threshold   = args.dg,
                top_n          = 999999,
                seed_required  = False
            )
            candidates.extend(res)

    # 5. Run Stage 2 ML Classification
    if candidates and args.model and args.traincsv:
        print("\n" + "="*60)
        print("  miRNAProtPred 2.0 — Stage 2: ML Contextual Scoring")
        print("="*60)
        candidates = stage2_score(candidates, args.model, args.traincsv)

    runtime = time.time() - start_time
    print(f"\nScan complete in {runtime:.2f} seconds.")

    # 6. Build Validation Summary YES/NO Table
    rows = []
    for req_id in requested_ids:
        # Find candidates for this request using hyphen-free matching
        req_clean = req_id.lower().replace("-", "").strip()
        id_candidates = []
        for c in candidates:
            c_clean = c['miRNA_ID'].lower().replace("-", "")
            if req_clean in c_clean or c_clean in req_clean:
                id_candidates.append(c)
        if id_candidates:
            # Sort by ML score descending, fallback to delta_G ascending
            id_candidates.sort(key=lambda x: (-x.get('ml_score', 0.0), x['delta_G']))
            best = id_candidates[0]
            ml_score_val = best.get('ml_score', 'N/A')
            rows.append({
                'miRNA_ID': req_id,
                'Matched_Database_ID': best['miRNA_ID'],
                'Validation_Result': 'YES',
                'Best_ml_score': ml_score_val,
                'Best_delta_G': best['delta_G'],
                'Locus': f"{best['start']}-{best['end']}"
            })
        else:
            rows.append({
                'miRNA_ID': req_id,
                'Matched_Database_ID': 'N/A',
                'Validation_Result': 'NO',
                'Best_ml_score': 0.0,
                'Best_delta_G': 0.0,
                'Locus': 'N/A'
            })

    summary_df = pd.DataFrame(rows)

    print("\n" + "="*60)
    print("  VALIDATION SUMMARY")
    print("="*60)
    print(summary_df.to_string(index=False))
    print("="*60)

    # 7. Print and save detailed output if raw or highconf is requested
    df_detail = None
    if args.output in ['raw', 'highconf']:
        if candidates:
            # Add confidence to all candidates if not present
            for c in candidates:
                if 'confidence' not in c:
                    p = c.get('ml_score', 0.0)
                    c['confidence'] = assign_confidence(p, c['delta_G'])

            display_candidates = list(candidates)
            if args.output == 'highconf':
                display_candidates = [c for c in candidates if c.get('confidence') in ['Very High', 'High', 'Medium']]
            
            if display_candidates:
                print("\n" + "="*60)
                print("  DETAILED INTERACTION ANALYSIS")
                print("="*60)
                cols = ['miRNA_ID', 'start', 'end', 'delta_G', 'seed_match', 'supp_match', 'cts_gc']
                if 'ml_score' in display_candidates[0]:
                    cols.append('ml_score')
                if 'confidence' in display_candidates[0]:
                    cols.append('confidence')
                df_detail = pd.DataFrame(display_candidates)[cols]
                print(df_detail.to_string(index=False))
                print("="*60)
            else:
                print("\nNo high-confidence detailed interactions to show.")
        else:
            print("\nNo detailed interactions to show (no candidates found).")

    # 8. Export to CSV
    if args.out:
        if args.output in ['raw', 'highconf'] and df_detail is not None:
            out_df = df_detail
        else:
            out_df = summary_df
        out_df.to_csv(args.out, index=False)
        print(f"\n  Saved verification results → {args.out}")

cli = main

if __name__ == '__main__':
    cli()
