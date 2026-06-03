import pandas as pd

def main():
    # Load missing targets
    df_missing_targets = pd.read_csv('output/audits/missing_target_sequences.csv')
    
    # Target recovery candidates mapping
    records = []
    
    for idx, row in df_missing_targets.iterrows():
        virus = row['virus']
        target_symbol = str(row['target_symbol']).strip()
        genome_acc = str(row['genome_accession']).strip()
        
        failure_reason = 'true annotation absence'
        recommended_alias = target_symbol
        
        sym_lower = target_symbol.lower()
        if sym_lower == 'ns5a':
            failure_reason = 'polyprotein component'
            recommended_alias = 'nonstructural protein 5A'
        elif sym_lower == 'nef':
            failure_reason = 'alternative gene name'
            recommended_alias = 'Nef protein'
        elif sym_lower == 'pb1':
            failure_reason = 'alternative gene name'
            recommended_alias = 'polymerase basic protein 1'
        elif sym_lower == 'core':
            failure_reason = 'alternative gene name'
            recommended_alias = 'capsid/core protein'
        elif sym_lower == 'poly':
            failure_reason = 'polyprotein component'
            recommended_alias = 'viral polyprotein'
        elif sym_lower == 'p':
            failure_reason = 'alternative gene name'
            recommended_alias = 'phosphoprotein'
        elif sym_lower == 'n':
            failure_reason = 'alternative gene name'
            recommended_alias = 'nucleoprotein'
        elif sym_lower == 'vp1':
            failure_reason = 'alternative gene name'
            recommended_alias = 'VP1'
        elif sym_lower == 'vp3':
            failure_reason = 'alternative gene name'
            recommended_alias = 'VP3'
        elif sym_lower == 'vp':
            failure_reason = 'alternative gene name'
            recommended_alias = 'VP'
        elif sym_lower == 'tas':
            failure_reason = 'alternative gene name'
            recommended_alias = 'transactivator protein'
        elif sym_lower == 'sfvgp1':
            failure_reason = 'alternative gene name'
            recommended_alias = 'glycoprotein'
        else:
            failure_reason = 'true annotation absence'
            recommended_alias = 'UNKNOWN'
            
        records.append({
            'virus': virus,
            'target_symbol': target_symbol,
            'genome_accession': genome_acc,
            'failure_reason': failure_reason,
            'recommended_alias': recommended_alias
        })
        
    audit_df = pd.DataFrame(records)
    audit_df.to_csv('output/audits/target_recovery_candidates.csv', index=False)
    
    recoverable_count = len(audit_df[audit_df['recommended_alias'] != 'UNKNOWN'])
    non_recoverable_count = len(audit_df[audit_df['recommended_alias'] == 'UNKNOWN'])
    
    print("=== Target Recovery Audit ===")
    print(f"Potentially recoverable targets: {recoverable_count}")
    print(f"Non-recoverable targets: {non_recoverable_count}")
    print(f"Estimated additional interactions recoverable: {recoverable_count}")

if __name__ == "__main__":
    main()
