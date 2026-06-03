#!/usr/bin/env python3
"""
clean_cts.py
────────────────────────────────────────────────────────────────────────────
Clean and deduplicate the raw CTS dataset:
1. Remove exact duplicate (miRNA_seq, CTS_seq) pairs
2. Remove cross-virus negatives that have the same CTS sequence as positives
   (these would be false negatives / leakage)
3. Ensure no cross-virus negative pairs identical miRNA sequence with the same
   CTS context as positives
4. Enforce minimum CTS length ≥ 10 nt
5. Export cleaned train-ready dataset with feature-ready columns
────────────────────────────────────────────────────────────────────────────
"""
import pandas as pd
import numpy as np

def gc_content(seq):
    seq = str(seq).upper().replace('T', 'U')
    if not seq: return 0.0
    gc = sum(1 for c in seq if c in 'GC')
    return round(gc / len(seq), 4)

def au_content(seq):
    seq = str(seq).upper().replace('T', 'U')
    if not seq: return 0.0
    au = sum(1 for c in seq if c in 'AU')
    return round(au / len(seq), 4)

def seed_match_score(mirna_seq, cts_seq):
    mirna = str(mirna_seq).upper().replace('T', 'U')
    cts_rna = str(cts_seq).upper().replace('T', 'U')
    cts_dna = str(cts_seq).upper().replace('U', 'T')
    seed = mirna[1:7]  # positions 2-7
    comp_rna = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    comp_dna = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rc_seed_rna = ''.join(comp_rna.get(b, 'N') for b in reversed(seed))
    rc_seed_dna = ''.join(comp_dna.get(b, 'N') for b in reversed(seed.replace('U','T')))
    return 1 if (rc_seed_rna in cts_rna or rc_seed_dna in cts_dna) else 0

def supplementary_score(mirna_seq, cts_seq):
    mirna = str(mirna_seq).upper().replace('T', 'U')
    cts_rna = str(cts_seq).upper().replace('T', 'U')
    cts_dna = str(cts_seq).upper().replace('U', 'T')
    if len(mirna) < 16: return 0
    supp = mirna[12:16]
    comp_rna = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    comp_dna = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rc_supp_rna = ''.join(comp_rna.get(b, 'N') for b in reversed(supp))
    rc_supp_dna = ''.join(comp_dna.get(b, 'N') for b in reversed(supp.replace('U','T')))
    return 1 if (rc_supp_rna in cts_rna or rc_supp_dna in cts_dna) else 0

def count_pairs(hybridDP):
    """Count base pairs from hybridDP dot-bracket notation."""
    if pd.isna(hybridDP) or not hybridDP:
        return 0
    return hybridDP.count('(')

def main():
    df = pd.read_csv('virbase_final_dataset/virbase_cts/cts_dataset_raw.csv')
    if 'deltaG' in df.columns and 'delta_G' not in df.columns:
        df.rename(columns={'deltaG': 'delta_G'}, inplace=True)
    print(f"Input: {len(df)} rows  ({(df['label']==1).sum()} pos, {(df['label']==0).sum()} neg)")

    # 1. Enforce minimum CTS length
    df['cts_len'] = df['CTS_sequence'].str.len()
    df = df[df['cts_len'] >= 10].copy()
    print(f"After length filter (≥10 nt): {len(df)} rows")

    # 2. Deduplicate on (miRNA_sequence, CTS_sequence) — keep first occurrence per label
    df = df.drop_duplicates(subset=['label', 'miRNA_sequence', 'virus', 'target_symbol', 'CTS_sequence']).copy()
    print(f"After deduplication: {len(df)} rows  ({(df['label']==1).sum()} pos, {(df['label']==0).sum()} neg)")

    # 3a. Cap hard negatives to max 2 per miRNA+target_symbol combination
    hard_neg_mask = (df['label'] == 0) & (df['target_type'] != 'cross_virus_negative')
    hard_negs = df[hard_neg_mask].copy()
    other_rows = df[~hard_neg_mask].copy()

    hard_negs_capped = hard_negs.groupby(
        ['miRNA', 'target_symbol'], group_keys=False
    ).apply(lambda x: x.head(2))

    df = pd.concat([other_rows, hard_negs_capped], ignore_index=True)
    print(f"After hard-negative capping (max 2 per pair): {len(df)} rows "
          f"({(df['label']==1).sum()} pos, {(df['label']==0).sum()} neg)")

    # 3. Remove cross-virus negatives whose CTS sequence is identical to any positive CTS
    pos_cts_set = set(df[df['label'] == 1]['CTS_sequence'].str.upper())
    neg_mask_leak = (df['label'] == 0) & df['CTS_sequence'].str.upper().isin(pos_cts_set)
    print(f"Leaking negatives removed (CTS = positive CTS): {neg_mask_leak.sum()}")
    df = df[~neg_mask_leak].copy()

    # 4. Feature extraction
    df['mirna_gc']         = df['miRNA_sequence'].apply(gc_content)
    df['cts_gc']           = df['CTS_sequence'].apply(gc_content)
    df['cts_au']           = df['CTS_sequence'].apply(au_content)
    df['cts_len']          = df['CTS_sequence'].str.len()
    df['mirna_len']        = df['miRNA_sequence'].str.len()
    df['seed_match']       = df.apply(lambda r: seed_match_score(r['miRNA_sequence'], r['CTS_sequence']), axis=1)
    df['supp_match']       = df.apply(lambda r: supplementary_score(r['miRNA_sequence'], r['CTS_sequence']), axis=1)
    df['n_base_pairs']     = df['hybridDP'].apply(count_pairs)
    df['delta_G_norm']     = df['delta_G'] / df['cts_len']  # energy per nt
    df['site_position_norm'] = df['CTS_start'] / df.apply(
        lambda r: max(len(str(r['CTS_sequence'])) + r['CTS_start'], 1), axis=1)

    # 5. Save cleaned datasets
    feature_cols = [
        'label', 'miRNA', 'miRNA_sequence', 'virus', 'target_symbol', 'target_type',
        'CTS_sequence', 'CTS_start', 'CTS_end', 'cts_len', 'delta_G', 'delta_G_norm',
        'hybridDP', 'n_base_pairs', 'mirna_gc', 'cts_gc', 'cts_au', 'mirna_len',
        'seed_match', 'supp_match', 'site_position_norm', 'PMID', 'Score'
    ]
    df = df[feature_cols]
    df.to_csv('virbase_final_dataset/virbase_cts/cts_dataset_clean.csv', index=False)

    pos_df = df[df['label'] == 1]
    neg_df = df[df['label'] == 0]
    pos_df.to_csv('virbase_final_dataset/virbase_cts/cts_positives_clean.csv', index=False)
    neg_df.to_csv('virbase_final_dataset/virbase_cts/cts_negatives_clean.csv', index=False)

    print("\n" + "=" * 55)
    print("  CLEANED CTS DATASET SUMMARY")
    print("=" * 55)
    print(f"  Total training rows:   {len(df)}")
    print(f"  Positive CTSs:         {len(pos_df)}")
    print(f"  Negative CTSs:         {len(neg_df)}")
    print(f"  Neg:Pos ratio:         {len(neg_df)/max(len(pos_df),1):.2f}:1")
    print()
    print("  Feature Summary (Positives):")
    print(f"    Mean delta_G:        {pos_df['delta_G'].mean():.2f} kcal/mol")
    print(f"    Mean CTS length:     {pos_df['cts_len'].mean():.1f} nt")
    print(f"    Seed match rate:     {pos_df['seed_match'].mean()*100:.1f}%")
    print(f"    Supp match rate:     {pos_df['supp_match'].mean()*100:.1f}%")
    print(f"    Mean base pairs:     {pos_df['n_base_pairs'].mean():.1f}")
    print()
    print("  Feature Summary (Negatives):")
    print(f"    Mean delta_G:        {neg_df['delta_G'].mean():.2f} kcal/mol")
    print(f"    Mean CTS length:     {neg_df['cts_len'].mean():.1f} nt")
    print(f"    Seed match rate:     {neg_df['seed_match'].mean()*100:.1f}%")
    print(f"    Supp match rate:     {neg_df['supp_match'].mean()*100:.1f}%")
    print(f"    Mean base pairs:     {neg_df['n_base_pairs'].mean():.1f}")
    print("=" * 55)
    print("  Output: virbase_final_dataset/virbase_cts/cts_dataset_clean.csv")

if __name__ == "__main__":
    main()
