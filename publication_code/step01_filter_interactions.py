import pandas as pd
import numpy as np

input_file = 'workspace_archive/miRNA_RNA_interactions.txt'
output_file = 'dataset/intermediate/human_miRNA_viral_targets_score70.csv'

# Load dataset
print("Loading dataset...")
df = pd.read_csv(input_file, sep='\t', dtype=str)

# 1. Total rows in original ViRBase
total_rows = len(df)

# 2. Rows with Homo sapiens host
df_human = df[df['Host Species'] == 'Homo sapiens']
rows_homo_sapiens = len(df_human)

# 3. Rows with hsa-miR entries
# Filter by Interactor1 Source = host, Interactor1 Category = miRNA, Interactor1 Symbol starts with hsa-miR
df_hsamir = df_human[
    (df_human['Interactor1 Source'] == 'host') &
    (df_human['Interactor1 Category'] == 'miRNA') &
    (df_human['Interactor1 Symbol'].str.startswith('hsa-miR', na=False))
]
rows_hsamir = len(df_hsamir)

# 4. Rows with Score 0.70–1.00
df_hsamir = df_hsamir.copy()
df_hsamir['Score_num'] = pd.to_numeric(df_hsamir['Score'], errors='coerce')
df_score = df_hsamir[(df_hsamir['Score_num'] >= 0.70) & (df_hsamir['Score_num'] <= 1.00)]
rows_score = len(df_score)

# 5. Final Human miRNA → Viral Target interactions
# Filter by Interactor2 Source = virus
df_final = df_score[df_score['Interactor2 Source'] == 'virus']
final_interactions = len(df_final)

unique_miRNAs = df_final['Interactor1 Symbol'].nunique()
unique_viruses = df_final['Virus Name'].nunique()
unique_targets = df_final['Interactor2 Symbol'].nunique()

# Print Summary Report
print("="*50)
print("Summary Report")
print("="*50)
print(f"Total rows in original ViRBase: {total_rows}")
print(f"Rows with Homo sapiens host: {rows_homo_sapiens}")
print(f"Rows with hsa-miR entries: {rows_hsamir}")
print(f"Rows with Score 0.70-1.00: {rows_score}")
print(f"Final Human miRNA -> Viral Target interactions: {final_interactions}")
print(f"Unique human miRNAs: {unique_miRNAs}")
print(f"Unique viruses: {unique_viruses}")
print(f"Unique viral targets: {unique_targets}")
print("="*50)

# Rename and select columns for output
df_out = pd.DataFrame()
df_out['ViRBase_ID'] = df_final['ViRBase ID']
df_out['Virus_Name'] = df_final['Virus Name']
df_out['Virus_Strain_Name'] = df_final['Virus Strain Name']
df_out['Taxonomy_ID'] = df_final['Taxonomy ID']
df_out['Virus_Family'] = df_final['Virus Family']
df_out['Host_Species'] = df_final['Host Species']
df_out['miRNA'] = df_final['Interactor1 Symbol']
df_out['miRNA_ID'] = df_final['Interactor1 ID']
df_out['Target_Symbol'] = df_final['Interactor2 Symbol']
df_out['Target_ID'] = df_final['Interactor2 ID']
df_out['PMID'] = df_final['PMID']
df_out['Score'] = df_final['Score']
df_out['Interactor1_Source'] = df_final['Interactor1 Source']
df_out['Interactor1_Category'] = df_final['Interactor1 Category']
df_out['Interactor2_Source'] = df_final['Interactor2 Source']
df_out['Interactor2_Category'] = df_final['Interactor2 Category']

df_out.to_csv(output_file, index=False)
print(f"Saved extracted data to {output_file}")

inspect_cols = [
    'Virus Name',
    'Interactor1 Symbol',
    'Interactor2 Source',
    'Interactor2 Category',
    'Interactor2 Symbol',
    'Score'
]

df_score[inspect_cols].head(200).to_csv(
    "virbase_manual_inspection.csv",
    index=False
)
print("Saved manual inspection data to virbase_manual_inspection.csv")
