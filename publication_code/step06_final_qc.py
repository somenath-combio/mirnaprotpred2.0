import pandas as pd
import numpy as np

def main():
    # Load sequences
    seq_df = pd.read_csv('dataset/intermediate/human_miRNA_viral_targets_sequences.csv')
    seq_df['miRNA_sequence'] = seq_df['miRNA_sequence'].fillna('')
    seq_df['target_sequence'] = seq_df['target_sequence'].fillna('')
    
    # Load scores
    score_df = pd.read_csv('dataset/intermediate/human_miRNA_viral_targets_score70.csv')
    score_sub = score_df[['miRNA', 'Virus_Name', 'Target_Symbol', 'PMID', 'Score']].copy()
    score_sub.rename(columns={
        'Virus_Name': 'virus', 
        'Target_Symbol': 'target_symbol'
    }, inplace=True)
    score_sub.drop_duplicates(subset=['miRNA', 'virus', 'target_symbol', 'PMID'], inplace=True)
    
    # Merge to get Score
    merged_df = pd.merge(seq_df, score_sub, on=['miRNA', 'virus', 'target_symbol', 'PMID'], how='left')
    
    # Filter for complete pairs
    complete_df = merged_df[(merged_df['miRNA_sequence'] != '') & (merged_df['target_sequence'] != '')].copy()
    
    # Export final dataset
    export_cols = ['miRNA', 'miRNA_sequence', 'virus', 'target_symbol', 'target_type', 'target_sequence', 'genome_accession', 'PMID', 'Score']
    complete_df[export_cols].to_csv('dataset/intermediate/final_ml_dataset.csv', index=False)
    
    # 1. Number of complete pairs
    num_complete_pairs = len(complete_df)
    
    # 2. Unique miRNAs
    unique_mirnas = complete_df['miRNA'].nunique()
    
    # 3. Unique viruses
    unique_viruses = complete_df['virus'].nunique()
    
    # 4. Unique target genes
    unique_targets = complete_df['target_symbol'].nunique()
    
    # 5. miRNA length statistics
    complete_df['mirna_len'] = complete_df['miRNA_sequence'].apply(len)
    mirna_len_min = complete_df['mirna_len'].min()
    mirna_len_max = complete_df['mirna_len'].max()
    mirna_len_mean = complete_df['mirna_len'].mean()
    
    # 6. Target sequence length statistics
    complete_df['target_len'] = complete_df['target_sequence'].apply(len)
    target_len_min = complete_df['target_len'].min()
    target_len_max = complete_df['target_len'].max()
    target_len_mean = complete_df['target_len'].mean()
    target_len_median = complete_df['target_len'].median()
    
    # 7. Target type distribution
    type_counts = complete_df['target_type'].value_counts()
    
    # 8. Duplicate sequence pairs
    # Count rows having exactly the same miRNA sequence AND same target sequence
    dup_pairs = complete_df.duplicated(subset=['miRNA_sequence', 'target_sequence'], keep=False)
    dup_count = dup_pairs.sum()
    unique_dup_pairs = len(complete_df[dup_pairs].drop_duplicates(subset=['miRNA_sequence', 'target_sequence']))
    
    print("=== FINAL DATASET QC REPORT ===")
    print(f"1. Number of complete pairs: {num_complete_pairs}")
    print(f"2. Unique miRNAs: {unique_mirnas}")
    print(f"3. Unique viruses: {unique_viruses}")
    print(f"4. Unique target genes: {unique_targets}")
    print(f"\n5. miRNA length statistics:")
    print(f"   - Min:  {mirna_len_min}")
    print(f"   - Max:  {mirna_len_max}")
    print(f"   - Mean: {mirna_len_mean:.2f}")
    print(f"\n6. Target sequence length statistics:")
    print(f"   - Min:    {target_len_min}")
    print(f"   - Max:    {target_len_max}")
    print(f"   - Mean:   {target_len_mean:.2f}")
    print(f"   - Median: {target_len_median}")
    print(f"\n7. Target type distribution:")
    for t_type, count in type_counts.items():
        print(f"   - {t_type}: {count}")
    print(f"\n8. Duplicate sequence pairs:")
    print(f"   - Total rows involved in duplicate pairs: {dup_count}")
    print(f"   - Unique duplicated sequence pairs: {unique_dup_pairs}")
    print("===============================")
    print("Exported: final_ml_dataset.csv")

if __name__ == "__main__":
    main()
