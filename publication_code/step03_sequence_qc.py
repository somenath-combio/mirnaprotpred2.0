import pandas as pd
import numpy as np

# Load data
seq_df = pd.read_csv('dataset/intermediate/human_miRNA_viral_targets_sequences.csv')
score_df = pd.read_csv('dataset/intermediate/human_miRNA_viral_targets_score70.csv')

# Merge to get miRNA_ID and Score back
# score_df has columns: ViRBase_ID,Virus_Name,Virus_Strain_Name,Taxonomy_ID,Virus_Family,Host_Species,miRNA,miRNA_ID,Target_Symbol,Target_ID,PMID,Score,Interactor1_Source,Interactor1_Category,Interactor2_Source,Interactor2_Category
# seq_df has columns: miRNA,miRNA_sequence,virus,target_symbol,target_type,target_sequence,genome_accession,PMID

score_df_subset = score_df[['miRNA', 'miRNA_ID', 'Virus_Name', 'Target_Symbol', 'PMID', 'Score']].copy()
score_df_subset.rename(columns={'Virus_Name': 'virus', 'Target_Symbol': 'target_symbol', 'Score': 'score'}, inplace=True)
score_df_subset.drop_duplicates(subset=['miRNA', 'virus', 'target_symbol', 'PMID'], inplace=True)

df = pd.merge(seq_df, score_df_subset, on=['miRNA', 'virus', 'target_symbol', 'PMID'], how='left')

# Replace NaN with empty strings for easier filtering
df['miRNA_sequence'] = df['miRNA_sequence'].fillna('')
df['target_sequence'] = df['target_sequence'].fillna('')

# 1. missing_miRNA_sequences.csv
missing_mirna = df[df['miRNA_sequence'] == '']
missing_mirna_cols = ['miRNA', 'miRNA_ID', 'virus', 'target_symbol', 'PMID', 'score']
missing_mirna[missing_mirna_cols].to_csv('output/audits/missing_miRNA_sequences.csv', index=False)

# 2. missing_target_sequences.csv
missing_target = df[df['target_sequence'] == '']
missing_target_cols = ['miRNA', 'virus', 'target_symbol', 'target_type', 'genome_accession', 'PMID', 'score']
missing_target[missing_target_cols].to_csv('output/audits/missing_target_sequences.csv', index=False)

# 3. missing_both_sequences.csv
missing_both = df[(df['miRNA_sequence'] == '') & (df['target_sequence'] == '')]
missing_both_cols = ['miRNA', 'miRNA_ID', 'virus', 'target_symbol', 'target_type', 'genome_accession', 'PMID', 'score']
missing_both[missing_both_cols].to_csv('output/audits/missing_both_sequences.csv', index=False)

# 4. recovery_summary_report.csv
total_processed = len(df)
recovered_mirna = len(df[df['miRNA_sequence'] != ''])
missing_mirna_count = len(missing_mirna)
recovered_target = len(df[df['target_sequence'] != ''])
missing_target_count = len(missing_target)
missing_both_count = len(missing_both)

unique_missing_mirnas = missing_mirna['miRNA'].nunique()
unique_missing_targets = missing_target['target_symbol'].nunique()
unique_viruses_affected = missing_target['virus'].nunique()

summary_data = {
    'Metric': [
        'Total interactions processed',
        'Total miRNA sequences recovered',
        'Total miRNA sequences missing',
        'Total target sequences recovered',
        'Total target sequences missing',
        'Total rows missing both',
        'Unique missing miRNAs',
        'Unique missing viral targets',
        'Unique viruses affected'
    ],
    'Count': [
        total_processed,
        recovered_mirna,
        missing_mirna_count,
        recovered_target,
        missing_target_count,
        missing_both_count,
        unique_missing_mirnas,
        unique_missing_targets,
        unique_viruses_affected
    ]
}

summary_df = pd.DataFrame(summary_data)
summary_df.to_csv('output/audits/recovery_summary_report.csv', index=False)

# Additional Analysis Printout
print("=== Additional Analysis ===")
print("\nTop 20 Missing miRNAs:")
print(missing_mirna['miRNA'].value_counts().head(20).to_string())

print("\nTop 20 Missing Viral Targets:")
print(missing_target['target_symbol'].value_counts().head(20).to_string())

print("\nTop 20 Viruses Associated with Failures (Target Sequence):")
print(missing_target['virus'].value_counts().head(20).to_string())
