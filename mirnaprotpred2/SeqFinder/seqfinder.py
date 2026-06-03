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
import shutil
from pathlib import Path


def assign_confidence(prob, dg):
    """
    4-tier confidence grading:
    Very High : prob >= 0.75 AND dG <= -10.0
    High      : prob >= 0.50
    Medium    : prob >= 0.25 AND dG <= -9.0
    Low       : prob < 0.25 OR dG > -9.0
    """
    if prob >= 0.75 and dg <= -10.0:
        return "Very High"
    if prob >= 0.50:
        return "High"
    if prob >= 0.25 and dg <= -9.0:
        return "Medium"
    return "Low"

INTARNA_BIN = shutil.which("IntaRNA") or "/home/somenath/miniconda3/envs/cts_env/bin/IntaRNA"

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

import os
import hashlib
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez, SeqIO

CACHE_DIR = ".seqfinder_cache"
Entrez.email = "sudipta@pusan.ac.kr"

def get_sequence_type(seq):
    if not seq:
        raise ValueError("Empty sequence")
    seq = seq.upper().strip()
    dna_chars = set("ATGCN")
    rna_chars = set("AUGCN")
    protein_chars = set("ACDEFGHIKLMNPQRSTVWY*")
    seq_chars = set(seq)
    if 'U' in seq_chars and seq_chars.issubset(rna_chars):
        return "RNA"
    if seq_chars.issubset(dna_chars):
        return "DNA"
    if seq_chars.issubset(protein_chars):
        return "Protein"
    raise ValueError("Invalid biological sequence")

def run_tblastn(protein_seq):
    seq_hash = hashlib.md5(protein_seq.encode()).hexdigest()[:8]
    os.makedirs(CACHE_DIR, exist_ok=True)
    cache_file = os.path.join(CACHE_DIR, f"blast_cache_{seq_hash}.xml")
    print("[Protein Mode] Initiating tblastn search against NCBI 'nt' database...")
    try:
        if os.path.exists(cache_file):
            print(f"[Protein Mode] Loading cached BLAST results from {cache_file}...")
            with open(cache_file, "r") as f:
                blast_records = NCBIXML.parse(f)
                blast_record = next(blast_records)
        else:
            print("[Protein Mode] Executing remote tblastn (this may take a while)...")
            result_handle = NCBIWWW.qblast("tblastn", "nt", protein_seq)
            xml_data = result_handle.read()
            with open(cache_file, "w") as f:
                f.write(xml_data)
            with open(cache_file, "r") as f:
                blast_records = NCBIXML.parse(f)
                blast_record = next(blast_records)
        if not blast_record.alignments:
            print("[Protein Mode] Error: No tblastn alignments found.")
            return None
        best_alignment = blast_record.alignments[0]
        best_hsp = best_alignment.hsps[0]
        accession = best_alignment.accession
        hit_start = min(best_hsp.sbjct_start, best_hsp.sbjct_end)
        hit_end = max(best_hsp.sbjct_start, best_hsp.sbjct_end)
        print(f"[Protein Mode] Found best hit: {best_alignment.title[:50]}...")
        print(f"[Protein Mode] Alignment coordinates: {hit_start} - {hit_end}")
        return {
            'accession': accession,
            'hit_start': hit_start,
            'hit_end': hit_end,
            'title': best_alignment.title
        }
    except Exception as e:
        print(f"[Protein Mode] Error during tblastn execution: {e}")
        return None

def retrieve_nucleotide_region(accession, start, end):
    fetch_start = start
    fetch_end = end
    print(f"[Protein Mode] Retrieving exact nucleotide alignment region {accession}:{fetch_start}-{fetch_end}...")
    try:
        handle = Entrez.efetch(
            db="nucleotide",
            id=accession,
            seq_start=fetch_start,
            seq_stop=fetch_end,
            rettype="fasta",
            retmode="text"
        )
        records = list(SeqIO.parse(handle, "fasta"))
        if not records:
            print("[Protein Mode] Error: Empty retrieval result.")
            return None
        seq = str(records[0].seq)
        print(f"[Protein Mode] Successfully retrieved {len(seq)} nt sequence.")
        return seq
    except Exception as e:
        print(f"[Protein Mode] Error during sequence retrieval: {e}")
        return None

def load_fasta(path):
    """Returns concatenated sequence from a FASTA file."""
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'): continue
            seq.append(line.upper())
    return ''.join(seq)

def load_target_sequence(genome_input):
    if os.path.isfile(genome_input):
        seq = load_fasta(genome_input)
    else:
        seq = genome_input.strip()
    seq_type = get_sequence_type(seq)
    print(f"\nDetected sequence type: {seq_type}")
    if seq_type == "Protein":
        if len(seq) < 20:
            raise ValueError("Protein sequence too short for reliable tblastn.")
        print("[Protein Mode] Protein sequence detected.")
        print("[Protein Mode] Initiating protein-guided nucleotide region discovery followed by RNA motif interaction analysis.")
        blast_info = run_tblastn(seq)
        if not blast_info:
            raise ValueError("BLAST search failed.")
        retrieved_seq = retrieve_nucleotide_region(
            blast_info['accession'],
            blast_info['hit_start'],
            blast_info['hit_end']
        )
        if not retrieved_seq:
            raise ValueError("Sequence retrieval failed.")
        return retrieved_seq.upper().replace('U', 'T')
    elif seq_type == "RNA":
        return seq.upper().replace('U', 'T')
    else:
        return seq.upper()

# ── main scan ──────────────────────────────────────────────────────

# ── main scan ──────────────────────────────────────────────────────

def scan_genome(mirna_id, mirna_seq, genome_seq, win_len=50, step=10,
                dg_threshold=-5.0, top_n=50, seed_required=False):
    """
    Slide a window across the genome and collect all candidate
    sites where IntaRNA finds binding with delta_G <= dg_threshold.

    Parameters
    ----------
    mirna_id       : str  — identifier for the miRNA
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

    print(f"  miRNA ID       : {mirna_id}")
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
            'miRNA_ID'   : mirna_id,
            'miRNA_seq'  : mirna_seq,
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


def scan_genome_multiple(mirnas, genome_seq, win_len=50, dg_threshold=-5.0, top_n=50):
    """
    Scan a genome against multiple miRNAs using a fast Boyer-Moore seed match pre-filter.
    """
    candidates = []
    genome_len = len(genome_seq)

    print(f"  Genome length  : {genome_len:,} nt")
    print(f"  Window length  : {win_len} nt")
    print(f"  Total miRNAs   : {len(mirnas):,}")
    print(f"  delta_G cutoff : <= {dg_threshold}")
    print(f"  Scanning multiple miRNAs with Boyer-Moore pre-filter...")

    for idx, (mirna_id, mirna_seq) in enumerate(mirnas):
        if (idx + 1) % 100 == 0 or idx + 1 == len(mirnas):
            print(f"    Processed {idx + 1}/{len(mirnas)} miRNAs...")

        if len(mirna_seq) < 8:
            continue

        # Get 7-mer seed (positions 2-8, 1-indexed) and its reverse complement
        seed = mirna_seq[1:8].upper().replace('U', 'T')
        seed_rc = rc(seed)

        # Fast search for the seed reverse complement in the genome
        positions = []
        start_search = 0
        while True:
            pos = genome_seq.find(seed_rc, start_search)
            if pos == -1:
                break
            positions.append(pos)
            start_search = pos + 1

        for pos in positions:
            # Situate the seed match near the 3' end of a 50 nt target window
            start_win = max(0, pos - 35)
            end_win = min(genome_len, start_win + win_len)
            if end_win - start_win < win_len:
                start_win = max(0, end_win - win_len)

            window = genome_seq[start_win:end_win]

            result = run_intarna(mirna_seq, window)
            if result is None:
                continue
            dg, hyb, nbp = result
            if dg > dg_threshold:
                continue

            candidates.append({
                'miRNA_ID'   : mirna_id,
                'miRNA_seq'  : mirna_seq,
                'start'      : start_win,
                'end'        : end_win,
                'window_seq' : window,
                'delta_G'    : round(dg, 3),
                'hybridDP'   : hyb,
                'n_base_pairs': nbp,
                'seed_match' : 1,
                'supp_match' : supp_match(mirna_seq, window),
                'cts_gc'     : round(gc_content(window), 4),
            })

    print(f"  Multi-scan complete. Total candidates found: {len(candidates)}")
    candidates.sort(key=lambda x: x['delta_G'])
    return candidates[:top_n]


def stage2_score(candidates, model_pkl, train_csv, target_name=None):
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

    # Dynamically detect target type and virus family
    if target_name is None:
        for arg in sys.argv:
            if any(arg.lower().endswith(ext) for ext in ['.fasta', '.fa', '.fna', '.txt']):
                target_name = arg
                break
        if target_name is None:
            target_name = sys.argv[-1] if len(sys.argv) > 1 else ""

    target_name_str = str(target_name).lower()
    
    # Detect virus family
    virus_family_col = None
    if "hiv" in target_name_str or "retroviridae" in target_name_str or "lentivirus" in target_name_str or "htlv" in target_name_str or "retrovirus" in target_name_str or "immunodeficiency" in target_name_str:
        virus_family_col = "virusfamily_Retroviridae"
    elif "influenza" in target_name_str or "orthomyxoviridae" in target_name_str or "flu" in target_name_str:
        virus_family_col = "virusfamily_Orthomyxoviridae"
    elif "ebola" in target_name_str or "filoviridae" in target_name_str or "marburg" in target_name_str or "zaire" in target_name_str:
        virus_family_col = "virusfamily_Filoviridae"
    elif "zika" in target_name_str or "dengue" in target_name_str or "flaviviridae" in target_name_str or "hcv" in target_name_str or "hepacivirus" in target_name_str or "hepatitis c" in target_name_str or "west nile" in target_name_str:
        virus_family_col = "virusfamily_Flaviviridae"
    elif "hepatitis b" in target_name_str or "hepadnaviridae" in target_name_str or "hbv" in target_name_str:
        virus_family_col = "virusfamily_Hepadnaviridae"
    elif "herpes" in target_name_str or "herpesviridae" in target_name_str or "ebv" in target_name_str or "hcmv" in target_name_str or "cytomegalovirus" in target_name_str or "epstein-barr" in target_name_str or "kshv" in target_name_str or "hhv" in target_name_str:
        virus_family_col = "virusfamily_Herpesviridae"
    elif "adenovirus" in target_name_str or "adenoviridae" in target_name_str:
        virus_family_col = "virusfamily_Adenoviridae"
    elif "papillomavirus" in target_name_str or "papillomaviridae" in target_name_str or "hpv" in target_name_str:
        virus_family_col = "virusfamily_Papillomaviridae"
    elif "coxsackie" in target_name_str or "picornaviridae" in target_name_str or "poliovirus" in target_name_str or "enterovirus" in target_name_str or "rhinovir" in target_name_str:
        virus_family_col = "virusfamily_Picornaviridae"
    elif "pneumovirus" in target_name_str or "pneumoviridae" in target_name_str or "rsv" in target_name_str or "respiratory syncytial" in target_name_str:
        virus_family_col = "virusfamily_Pneumoviridae"
    elif "rhabdoviridae" in target_name_str or "rabies" in target_name_str or "vesicular" in target_name_str or "vsiv" in target_name_str:
        virus_family_col = "virusfamily_Rhabdoviridae"
    elif "togaviridae" in target_name_str or "sindbis" in target_name_str or "rubella" in target_name_str or "chikungunya" in target_name_str or "eeev" in target_name_str:
        virus_family_col = "virusfamily_Togaviridae"
    elif "bornaviridae" in target_name_str or "borna" in target_name_str or "bdv" in target_name_str:
        virus_family_col = "virusfamily_Bornaviridae"
    elif "arteriviridae" in target_name_str or "arterivirus" in target_name_str:
        virus_family_col = "virusfamily_Arteriviridae"

    # Detect target type
    target_type_col = "targettype_region"
    if any(k in target_name_str for k in ["gene", "nef", "gag", "tat", "rev", "vp40", "gp", "ns5a", "orf", "poly", "ul122", "cds"]):
        target_type_col = "targettype_gene"
    elif any(k in target_name_str for k in ["mrna", "transcript", "utr", "pgrna"]):
        target_type_col = "targettype_transcript"

    def dinuc(seq):
        bases = 'ACGU'
        seq = seq.upper().replace('T', 'U')
        d = {a+b: 0 for a in bases for b in bases}
        total = max(len(seq)-1, 1)
        for i in range(len(seq)-1):
            k = seq[i:i+2]
            if k in d: d[k] += 1
        return {f'dinuc_{k}': v/total for k,v in d.items()}

    def trinuc(seq):
        seq = seq.upper().replace('T', 'U')
        cats = {'AAA':0, 'UUU':0, 'CCC':0, 'GGG':0, 'GCC':0, 'CGC':0, 'GCG':0, 'CGU':0}
        total = max(len(seq)-2, 1)
        for i in range(len(seq)-2):
            t = seq[i:i+3]
            if t in cats: cats[t] += 1
        return {f'trinuc_{k}': v/total for k,v in cats.items()}

    rows = []
    for c in candidates:
        row = {}
        row['cts_len'] = len(c['window_seq'])
        row['mirna_len'] = len(c.get('miRNA_seq', ''))
        row['cts_gc'] = c['cts_gc']
        row['cts_au'] = sum(1 for b in c['window_seq'].upper() if b in 'ATU') / max(len(c['window_seq']),1)
        row['mirna_gc'] = gc_content(c.get('miRNA_seq', ''))
        row['seed_match'] = c['seed_match']
        row['supp_match'] = c['supp_match']
        row['motif_identity'] = c['seed_match'] * c['supp_match']
        
        # Populate target type
        row['targettype_gene'] = 1 if target_type_col == "targettype_gene" else 0
        row['targettype_transcript'] = 1 if target_type_col == "targettype_transcript" else 0
        row['targettype_region'] = 1 if target_type_col == "targettype_region" else 0
        
        # Populate virus family
        for vf in [c for c in feat_cols if c.startswith("virusfamily_")]:
            row[vf] = 1 if vf == virus_family_col else 0

        row.update(dinuc(c['window_seq']))
        row.update(trinuc(c['window_seq']))
        rows.append(row)

    df = pd.DataFrame(rows)
    for col in feat_cols:
        if col not in df.columns: df[col] = 0
    df = df[feat_cols].fillna(0)

    probs = model.predict_proba(df)[:,1]
    for i, c in enumerate(candidates):
        p = round(float(probs[i]), 4)
        c['ml_score'] = p
        c['confidence'] = assign_confidence(p, c['delta_G'])

    candidates.sort(key=lambda x: -x['ml_score'])
    return candidates


# ── CLI ────────────────────────────────────────────────────────────

def main():
    args_list = sys.argv[1:]
    # Remove '--genome' if present to support legacy keyword parameters as positional
    cleaned_args = []
    skip_next = False
    for arg in args_list:
        if skip_next:
            skip_next = False
            continue
        if arg == '--genome':
            skip_next = True
            continue
        cleaned_args.append(arg)

    parser = argparse.ArgumentParser(
        prog='SeqFinder2',
        description='SeqFinder + miRNAProtPred 2.0: two-stage viral miRNA target scanner'
    )
    # The first positional argument: genome (which is the target FASTA or sequence)
    parser.add_argument('genome',
                        help='Path to viral mRNA/genome FASTA or raw sequence')
    parser.add_argument('--mirna',   default='ALL',
                        help='Mature miRNA sequence (A/T/U/G/C), path to a FASTA file of miRNAs, or "ALL" to scan the default database (default: ALL)')
    parser.add_argument('--win',     type=int, default=50,
                        help='Sliding window length (default 50 nt)')
    parser.add_argument('--step',    type=int, default=10,
                        help='Step size (default 10 nt)')
    parser.add_argument('--dg',      type=float, default=-5.0,
                        help='delta_G threshold, e.g. -5.0 (default)')
    parser.add_argument('--top',     type=int, default=50,
                        help='Return top N candidates (default 50)')
    parser.add_argument('--mode',    choices=['strict', 'relaxed'], default='strict',
                        help='Search mode: strict enforces 7-mer seed matches, relaxed allows non-canonical matches (default: strict)')
    parser.add_argument('--output',  choices=['concise', 'raw', 'highconf', 'highconf_raw'], default='concise',
                        help='Output mode (default: concise)')
    parser.add_argument('--out',     default=None,
                        help='Save results automatically to the specified CSV filename')
    pkg_data_dir = Path(__file__).resolve().parent / "data"
    default_model = pkg_data_dir / "mirnaprotpred2.0_best.pkl"
    default_train = pkg_data_dir / "cts_ml_features_v4_intraviral.csv"

    parser.add_argument('--model',   default=str(default_model) if default_model.exists() else None,
                        help='Path to mirnaprotpred2.0_best.pkl for Stage 2 scoring')
    parser.add_argument('--traincsv',default=str(default_train) if default_train.exists() else None,
                        help='Path to training CSV (for feature column alignment)')
    parser.add_argument('--email',   default='sudipta@pusan.ac.kr',
                        help='NCBI Entrez email for Protein mode blast searches')
    
    args = parser.parse_args(cleaned_args)

    print("="*60)
    print("  SeqFinder — Stage 1: Thermodynamic Genome Scan")
    print("="*60)
    print(f"  miRNA Mode: {args.mirna}")
    print(f"  Genome    : {args.genome}")
    print(f"  Mode      : {args.mode}")
    print(f"  Output    : {args.output}")
    # Verify IntaRNA executable exists and is runnable
    if not shutil.which(INTARNA_BIN):
        print(f"Error: IntaRNA executable not found at '{INTARNA_BIN}'", file=sys.stderr)
        print("Please ensure IntaRNA is installed and available in your PATH or current conda environment.", file=sys.stderr)
        sys.exit(1)

    if args.email:
        Entrez.email = args.email

    genome_seq = load_target_sequence(args.genome)

    # Resolve miRNAs to scan
    mirnas = []
    if args.mirna.upper() == "ALL":
        db_path = pkg_data_dir / "data.xlsx"
        if not db_path.exists():
            print(f"Error: Default database not found at '{db_path}'", file=sys.stderr)
            sys.exit(1)
        import pandas as pd
        print(f"Loading miRNA database from '{db_path}'...")
        df_db = pd.read_excel(db_path, engine='openpyxl')
        for _, row in df_db.iterrows():
            mirna_id = str(row.get('Human miRNA ID', 'unknown')).strip()
            mirna_seq = str(row.get('Sequence', '')).strip().upper()
            if mirna_seq:
                mirnas.append((mirna_id, mirna_seq))
    elif Path(args.mirna).exists():
        from Bio import SeqIO
        print(f"Loading miRNAs from FASTA file '{args.mirna}'...")
        for record in SeqIO.parse(args.mirna, "fasta"):
            mirnas.append((record.id, str(record.seq).strip().upper()))
    else:
        mirnas.append(("Query_miRNA", args.mirna.strip().upper()))

    seed_required = (args.mode == 'strict')

    if len(mirnas) > 1:
        if not seed_required:
            print("[Warning] Database-wide scanning is restricted to strict mode for performance.")
        candidates = scan_genome_multiple(
            mirnas         = mirnas,
            genome_seq     = genome_seq,
            win_len        = args.win,
            dg_threshold   = args.dg,
            top_n          = args.top,
        )
    else:
        mirna_id, mirna_seq = mirnas[0]
        candidates = scan_genome(
            mirna_id       = mirna_id,
            mirna_seq      = mirna_seq,
            genome_seq     = genome_seq,
            win_len        = args.win,
            step           = args.step,
            dg_threshold   = args.dg,
            top_n          = args.top,
            seed_required  = seed_required,
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

    # Assign fallback confidence/ml_score if Stage 2 didn't run
    for c in candidates:
        if 'ml_score' not in c:
            c['ml_score'] = 0.0
        if 'confidence' not in c:
            c['confidence'] = assign_confidence(0.0, c['delta_G'])

    # Filter and format results based on --output
    concise_cols = ['miRNA_ID', 'start', 'end', 'delta_G', 'seed_match', 'cts_gc', 'ml_score', 'confidence']
    display_candidates = []
    
    if args.output in ['concise', 'highconf']:
        # Filter to High/Very High for concise, High/Very High/Medium for highconf
        allowed_tiers = ['Very High', 'High'] if args.output == 'concise' else ['Very High', 'High', 'Medium']
        display_candidates = [c for c in candidates if c['confidence'] in allowed_tiers]
    elif args.output == 'highconf_raw':
        display_candidates = [c for c in candidates if c['confidence'] in ['Very High', 'High', 'Medium']]
    else: # raw
        display_candidates = list(candidates)

    if not display_candidates:
        print("No candidates matched the output criteria.")
        sys.exit(0)

    import pandas as pd
    df_display = pd.DataFrame(display_candidates)
    if args.output in ['concise', 'highconf']:
        existing_cols = [col for col in concise_cols if col in df_display.columns]
        df_display = df_display[existing_cols]

    print("\nFinal Results:\n")
    print(df_display.to_string(index=False, max_rows=100))
    print("="*60)

    if args.out:
        out_path = Path(args.out)
        df_display.to_csv(out_path, index=False)
        print(f"  Results saved → {out_path}")
        print("="*60)

cli = main

if __name__ == '__main__':
    cli()
