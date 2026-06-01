#!/usr/bin/env python3
"""
feature_extraction_v2.py
miRNAProtPred 2.0 — Full 58-dimension feature extraction
Runs on BOTH training set and external validation set.

Usage:
  python scripts/feature_extraction_v2.py --mode train
  python scripts/feature_extraction_v2.py --mode external
  python scripts/feature_extraction_v2.py --mode both
"""

import argparse
import pandas as pd
import numpy as np
import re
from pathlib import Path

ROOT         = Path(__file__).resolve().parent.parent
TRAIN_CSV    = ROOT / "virbase_final_dataset/virbase_cts/mirnaprotpred2_training_set.csv"
EXTERNAL_CSV = ROOT / "virmirna_external/virmirna_external_features.csv"
TRAIN_OUT    = ROOT / "virbase_final_dataset/virbase_cts/cts_ml_features_v2.csv"
EXTERNAL_OUT = ROOT / "virmirna_external/virmirna_ml_features_v2.csv"

DINUCLEOTIDES  = ["AA","AC","AG","AU","CA","CC","CG","CU","GA","GC","GG","GU","UA","UC","UG","UU"]
TRINUCLEOTIDES = ["AAA","UUU","CCC","GGG","GCC","CGC","GCG","CGU"]

VIRUS_FAMILY_MAP = {
    "Hepacivirus C":"Flaviviridae","HCV":"Flaviviridae","Hepacivirus C (HCV)":"Flaviviridae",
    "Zika virus":"Flaviviridae","ZIKV":"Flaviviridae",
    "Hepatitis B virus":"Hepadnaviridae","HBV":"Hepadnaviridae","Hepatitis B virus (HBV)":"Hepadnaviridae",
    "Human immunodeficiency virus 1":"Retroviridae","HIV-1":"Retroviridae","HIV":"Retroviridae",
    "Simian immunodeficiency virus":"Retroviridae","SIV":"Retroviridae",
    "Human T-lymphotropic virus 1":"Retroviridae","HTLV-1":"Retroviridae",
    "Influenza A virus (H1N1)":"Orthomyxoviridae","Influenza A virus (H3N2)":"Orthomyxoviridae",
    "Influenza A virus (H5N1)":"Orthomyxoviridae","INFV":"Orthomyxoviridae","IAV":"Orthomyxoviridae",
    "Coxsackievirus B3":"Picornaviridae","CVB3":"Picornaviridae","Coxsackievirus B3 (CVB3)":"Picornaviridae",
    "Enterovirus A71":"Picornaviridae","EV71":"Picornaviridae","EV-A71":"Picornaviridae",
    "Human papillomavirus type 16":"Papillomaviridae","HPV16":"Papillomaviridae",
    "Human papillomavirus type 18":"Papillomaviridae","HPV18":"Papillomaviridae",
    "Human betaherpesvirus 5":"Herpesviridae","HHV-5":"Herpesviridae","HCMV":"Herpesviridae",
    "Human gammaherpesvirus 4":"Herpesviridae","EBV":"Herpesviridae","EPV":"Herpesviridae",
    "Human gammaherpesvirus 8":"Herpesviridae","KSHV":"Herpesviridae",
    "Varicella-zoster virus":"Herpesviridae","VZV":"Herpesviridae",
    "Human orthopneumovirus":"Pneumoviridae","HRSV":"Pneumoviridae","RSV":"Pneumoviridae",
    "Borna disease virus":"Bornaviridae","BDV":"Bornaviridae",
    "Zaire ebolavirus":"Filoviridae","ZEBOV":"Filoviridae",
    "Vesicular stomatitis Indiana virus":"Rhabdoviridae","VSIV":"Rhabdoviridae","VSV":"Rhabdoviridae",
    "Eastern equine encephalitis virus":"Togaviridae","EEEV":"Togaviridae","Sindbis virus":"Togaviridae",
    "Porcine reproductive and respiratory syndrome virus":"Arteriviridae","PRRSV":"Arteriviridae",
    "Human mastadenovirus C":"Adenoviridae","ADENOVIRUS":"Adenoviridae",
}
VIRUS_FAMILIES = sorted(set(VIRUS_FAMILY_MAP.values()))

def to_rna(seq):
    return str(seq).upper().replace("T","U")

def gc_content(seq):
    seq = to_rna(seq)
    return round(sum(1 for c in seq if c in "GC") / max(1,len(seq)), 4)

def au_content(seq):
    seq = to_rna(seq)
    return round(sum(1 for c in seq if c in "AU") / max(1,len(seq)), 4)

def get_dinuc_freq(seq):
    seq = to_rna(seq)
    freqs = {d: 0.0 for d in DINUCLEOTIDES}
    if len(seq) < 2: return freqs
    for i in range(len(seq)-1):
        d = seq[i:i+2]
        if d in freqs: freqs[d] += 1
    total = len(seq)-1
    return {k: round(v/total,4) for k,v in freqs.items()}

def get_trinuc_freq(seq):
    seq = to_rna(seq)
    freqs = {t: 0.0 for t in TRINUCLEOTIDES}
    if len(seq) < 3: return freqs
    for i in range(len(seq)-2):
        t = seq[i:i+3]
        if t in freqs: freqs[t] += 1
    total = len(seq)-2
    return {k: round(v/total,4) for k,v in freqs.items()}

def get_virus_family(virus_str):
    v = str(virus_str).strip()
    if v in VIRUS_FAMILY_MAP: return VIRUS_FAMILY_MAP[v]
    for key,fam in VIRUS_FAMILY_MAP.items():
        if key.lower() in v.lower() or v.lower() in key.lower(): return fam
    return "UnknownFamily"

def get_target_type_onehot(ttype_str):
    t = str(ttype_str).lower().strip()
    if "transcript" in t: return 0.0, 1.0, 0.0
    if any(x in t for x in ["utr","region","coding","3'"]): return 0.0, 0.0, 1.0
    return 1.0, 0.0, 0.0

def parse_hybridDP(hybridDP_str):
    db = str(hybridDP_str) if pd.notna(hybridDP_str) else ""
    if "&" not in db:
        return {"matches":0,"wobbles":0,"mismatches":0,"bulges":0,"loops":0}
    parts    = db.split("&")
    db_mirna = parts[0]
    db_cts   = parts[1] if len(parts)>1 else ""
    total_pairs = db_mirna.count("(")
    est_wobbles = max(0, total_pairs // 7)
    return {
        "matches":    max(0, total_pairs - est_wobbles),
        "wobbles":    est_wobbles,
        "mismatches": db_mirna.count("."),
        "bulges":     db_mirna.count("x") + db_cts.count("x"),
        "loops":      db_mirna.count("o"),
    }

def compute_motif_identity(mirna_seq, cts_seq):
    mirna = to_rna(mirna_seq).replace("U","T")
    cts   = cts_seq.upper().replace("U","T")
    if len(mirna) < 7 or len(cts) < 6: return 0.85
    seed    = mirna[1:7]
    rc_seed = seed.translate(str.maketrans("ATGC","TACG"))[::-1]
    best = 0
    for i in range(len(cts)-5):
        m = sum(1 for a,b in zip(rc_seed, cts[i:i+6]) if a==b)
        if m > best: best = m
    return round(best/6.0, 4)

def extract_features(row):
    mseq  = to_rna(row["miRNA_sequence"])
    cseq  = to_rna(row["CTS_sequence"])
    virus = str(row["virus"])
    ttype = str(row.get("target_type","gene"))

    delta_g      = float(row["delta_G"])     if pd.notna(row.get("delta_G"))     else -5.0
    cts_len      = int(row["cts_len"])        if pd.notna(row.get("cts_len"))     else len(cseq)
    mirna_len    = int(row["mirna_len"])      if pd.notna(row.get("mirna_len"))   else len(mseq)
    delta_g_norm = round(delta_g/max(1,cts_len), 4)
    n_bp         = int(row["n_base_pairs"])   if pd.notna(row.get("n_base_pairs")) else 0
    cts_start    = float(row["CTS_start"])    if pd.notna(row.get("CTS_start"))   else 1.0
    cts_end      = float(row["CTS_end"])      if pd.notna(row.get("CTS_end"))     else float(cts_len)
    site_pos     = round(cts_start/max(1.0,cts_end), 4)

    seed_match   = int(row["seed_match"])     if pd.notna(row.get("seed_match"))  else 0
    supp_match   = int(row["supp_match"])     if pd.notna(row.get("supp_match"))  else 0
    motif_id     = compute_motif_identity(mseq, cseq)

    pf      = parse_hybridDP(row.get("hybridDP",""))
    family  = get_virus_family(virus)
    tt_gene, tt_transcript, tt_region = get_target_type_onehot(ttype)
    dinucs  = get_dinuc_freq(cseq)
    trinucs = get_trinuc_freq(cseq)

    feat = {
        "label":         int(row["label"]),
        "miRNA":         str(row["miRNA"]),
        "virus":         virus,
        "target_symbol": str(row.get("target_symbol","")),
        "delta_G":           delta_g,
        "delta_G_norm":      delta_g_norm,
        "cts_len":           cts_len,
        "mirna_len":         mirna_len,
        "n_base_pairs":      n_bp,
        "site_position_norm": site_pos,
        "mirna_gc":          gc_content(mseq),
        "cts_gc":            gc_content(cseq),
        "cts_au":            au_content(cseq),
        "seed_match":        seed_match,
        "supp_match":        supp_match,
        "motif_identity":    motif_id,
        "struct_matches":    pf["matches"],
        "struct_wobbles":    pf["wobbles"],
        "struct_mismatches": pf["mismatches"],
        "struct_bulges":     pf["bulges"],
        "struct_loops":      pf["loops"],
        "targettype_gene":       tt_gene,
        "targettype_transcript": tt_transcript,
        "targettype_region":     tt_region,
    }
    feat.update({f"virusfamily_{fam}": (1.0 if family==fam else 0.0) for fam in VIRUS_FAMILIES})
    feat.update({f"dinuc_{k}": v for k,v in dinucs.items()})
    feat.update({f"trinuc_{k}": v for k,v in trinucs.items()})
    return feat

def main():
    parser = argparse.ArgumentParser(description="miRNAProtPred 2.0 Feature Extraction")
    parser.add_argument("--mode", choices=["train","external","both"], default="both")
    args = parser.parse_args()

    datasets = []
    if args.mode in ("train","both"):    datasets.append(("TRAINING", TRAIN_CSV, TRAIN_OUT))
    if args.mode in ("external","both"): datasets.append(("EXTERNAL", EXTERNAL_CSV, EXTERNAL_OUT))

    for name, inp, out in datasets:
        print(f"\n{'='*60}\n  {name}: {inp}\n{'='*60}")
        df = pd.read_csv(inp)
        print(f"  Loaded {len(df)} rows")
        features = []
        for idx, row in df.iterrows():
            try:
                features.append(extract_features(row))
            except Exception as e:
                print(f"  WARNING row {idx}: {e}")
        feat_df = pd.DataFrame(features)
        feat_df.to_csv(out, index=False)
        print(f"  Rows:          {len(feat_df)}")
        print(f"  Feature dims:  {feat_df.shape[1]-4}")
        print(f"  Label dist:    {feat_df['label'].value_counts().to_dict()}")
        print(f"  Saved → {out}")

    print(f"\n{'='*60}")
    print("  miRNAProtPred 2.0 — FEATURE EXTRACTION COMPLETE")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()
