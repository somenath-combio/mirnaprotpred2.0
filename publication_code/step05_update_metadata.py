import pandas as pd
import time
from Bio import Entrez, SeqIO
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
Entrez.email = "somenath@example.com"

def fetch_updated_sequence(accession, recommended_alias):
    try:
        time.sleep(0.35)
        handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        
        target_lower = recommended_alias.lower()
        if target_lower == 'nonstructural protein 5a':
            target_lower = 'ns5a protein'
        elif target_lower == 'capsid/core protein':
            target_lower = 'core protein'
        elif target_lower == 'viral polyprotein':
            target_lower = 'polyprotein'
        
        for feature in record.features:
            if feature.type in ["CDS", "gene", "mat_peptide", "Protein"]:
                gene = feature.qualifiers.get("gene", [""])[0].lower()
                product = feature.qualifiers.get("product", [""])[0].lower()
                locus_tag = feature.qualifiers.get("locus_tag", [""])[0].lower()
                note = feature.qualifiers.get("note", [""])[0].lower()
                
                # Check for direct or partial match
                if target_lower == gene or target_lower == locus_tag:
                    return str(feature.extract(record.seq)).upper()
                if target_lower in product or target_lower in note:
                    return str(feature.extract(record.seq)).upper()
        return ""
    except Exception as e:
        logging.error(f"Error fetching {accession}: {e}")
        return ""

def get_accession_for_virus(virus_name):
    mapping = {
        'Hepacivirus C (HCV)': 'NC_004102',
        'Borna disease virus (BDV)': 'NC_001607',
        'Enterovirus A71 (EV-A71)': 'NC_001612',
        'Human foamy virus (SFV-11)': 'NC_001736',
        'Semliki Forest virus 4': 'NC_003433',
        'Simian immunodeficiency virus (SIV)': 'NC_001549',
        'Human immunodeficiency virus 1 (HIV-1)': 'NC_001802',
        'Human immunodeficiency virus (HIV)': 'NC_001802',
        'Influenza A virus (H1N1)': 'NC_026437', # PB1 segment for H1N1
        'Influenza A virus (H5N1)': 'NC_007360'  # PB1 segment for H5N1
    }
    if virus_name in mapping:
        return mapping[virus_name]
    return ""

def main():
    # Load mapping
    audit_df = pd.read_csv('output/audits/target_recovery_candidates.csv')
    candidates = audit_df[audit_df['recommended_alias'] != 'UNKNOWN']
    
    # Load main dataset
    main_df = pd.read_csv('dataset/intermediate/human_miRNA_viral_targets_sequences.csv')
    main_df['miRNA_sequence'] = main_df['miRNA_sequence'].fillna('')
    main_df['target_sequence'] = main_df['target_sequence'].fillna('')
    
    # Track metrics
    recovered_count = 0
    
    cache = {}
    
    for idx, row in candidates.iterrows():
        virus = row['virus']
        target_symbol = row['target_symbol']
        genome_acc = row['genome_accession']
        rec_alias = row['recommended_alias']
        
        is_missing_acc = pd.isna(genome_acc) or genome_acc == '' or genome_acc == 'nan'
        
        # Find matching rows in main_df
        if is_missing_acc:
            mask = (main_df['virus'] == virus) & (main_df['target_symbol'] == target_symbol) & (main_df['genome_accession'].isna())
        else:
            mask = (main_df['virus'] == virus) & (main_df['target_symbol'] == target_symbol) & (main_df['genome_accession'] == genome_acc)
        
        if mask.any():
            if (main_df.loc[mask, 'target_sequence'] == '').all():
                actual_acc = get_accession_for_virus(virus) if is_missing_acc else genome_acc
                
                if not actual_acc:
                    logging.warning(f"No accession found for {virus}")
                    continue
                    
                cache_key = (actual_acc, rec_alias)
                if cache_key not in cache:
                    logging.info(f"Fetching {rec_alias} from {actual_acc}...")
                    cache[cache_key] = fetch_updated_sequence(actual_acc, rec_alias)
                
                new_seq = cache[cache_key]
                
                if new_seq:
                    main_df.loc[mask, 'target_sequence'] = new_seq
                    recovered_count += len(main_df[mask])
                    logging.info(f"Successfully recovered sequence for {target_symbol} ({rec_alias}).")
                else:
                    logging.warning(f"Failed to find {rec_alias} in {genome_acc}.")
                    
    # Fill NaN with empty string just to be clean
    main_df['miRNA_sequence'] = main_df['miRNA_sequence'].fillna('')
    main_df['target_sequence'] = main_df['target_sequence'].fillna('')
                    
    main_df.to_csv('dataset/intermediate/human_miRNA_viral_targets_sequences.csv', index=False)
    
    # Final metrics
    total_valid = len(main_df[(main_df['miRNA_sequence'] != '') & (main_df['target_sequence'] != '')])
    
    logging.info("================ SUMMARY ================")
    logging.info(f"Target sequences recovered in this pass: {recovered_count}")
    logging.info(f"Final fully valid interactions (miRNA + Target seq): {total_valid}")
    logging.info("=========================================")

if __name__ == "__main__":
    main()
