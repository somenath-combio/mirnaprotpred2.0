import pandas as pd
import urllib.request
import gzip
import os
import time
import logging
import re
import xml.etree.ElementTree as ET
from Bio import Entrez, SeqIO

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
Entrez.email = "somenath@example.com"

INPUT_CSV = "dataset/intermediate/human_miRNA_viral_targets_score70.csv"
OUTPUT_CSV = "dataset/intermediate/human_miRNA_viral_targets_sequences.csv"
MATURE_FA_FILE = "workspace_archive/mature.fa"

ALIAS_MAP = {
    "tat": "Tat",
    "rev": "Rev",
    "nef": "Nef",
    "vpr": "Vpr",
    "vif": "Vif",
    "vpu": "Vpu",
    "gag": "Gag",
    "pol": "Pol",
    "env": "Env",
    "X": "HBx",
    "C": "HBc",
    "P": "Polymerase",
    "S": "HBs",
    "ORF50": "ORF50",
    "BHRF1": "BHRF1",
    "tax": "Tax",
    "UL122": "UL122",
    "PRRSVgp4": "PRRSVgp4",
    "PRRSVgp6": "PRRSVgp6",
    "PRRSVgp8": "PRRSVgp8",
    "POLY": "polyprotein",
    "K2": "K2"
}

def download_mirbase():
    """Download and parse mature.fa from miRBase."""
    if not os.path.exists(MATURE_FA_FILE):
        logging.info("Downloading mature.fa...")
        os.makedirs(os.path.dirname(MATURE_FA_FILE), exist_ok=True)
        try:
            urllib.request.urlretrieve("https://mirbase.org/download/mature.fa", MATURE_FA_FILE)
        except:
            urllib.request.urlretrieve("ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz", MATURE_FA_FILE + ".gz")
            with gzip.open(MATURE_FA_FILE + ".gz", 'rb') as f_in:
                with open(MATURE_FA_FILE, 'wb') as f_out:
                    f_out.write(f_in.read())
            os.remove(MATURE_FA_FILE + ".gz")
            
    mirna_dict = {}
    for record in SeqIO.parse(MATURE_FA_FILE, "fasta"):
        parts = record.description.split()
        for p in parts:
            if p.startswith("MIMAT"):
                mirna_dict[p] = str(record.seq).upper().replace("T", "U")
                break
    return mirna_dict

def get_genome_from_taxonomy(tax_id):
    try:
        time.sleep(0.35)
        handle = Entrez.esearch(db="nuccore", term=f"txid{tax_id}[ORGN] AND refseq[filter]", retmax=10)
        record = Entrez.read(handle)
        handle.close()
        id_list = record.get("IdList", [])
        if id_list:
            time.sleep(0.35)
            handle = Entrez.esummary(db="nuccore", id=",".join(id_list))
            summaries = Entrez.read(handle)
            handle.close()
            
            for summary in summaries:
                if summary["Caption"].startswith("NC_"):
                    return summary["Caption"]
            return summaries[0]["Caption"]
    except:
        pass
    return None

def classify_target(symbol):
    if symbol.upper() == 'POLY':
        return 'gene'
    sym_lower = symbol.lower()
    if 'rna' in sym_lower or 'poly' in sym_lower or 'pgrna' in sym_lower:
        return 'transcript'
    elif 'nt' in sym_lower or 'ntr' in sym_lower:
        return 'region'
    else:
        return 'gene'

def fetch_target_sequence(accession, raw_symbol, target_type):
    try:
        time.sleep(0.35)
        handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        
        full_seq = str(record.seq).upper()
        
        if target_type == 'transcript':
            return full_seq
            
        elif target_type == 'region':
            # e.g., nt1368-nt1383 or CVB3 RNA(nt4343-nt4370)
            matches = re.findall(r'nt(\d+)-nt(\d+)', raw_symbol.lower())
            if matches:
                start = int(matches[0][0])
                end = int(matches[0][1])
                # Convert 1-based inclusive to 0-based exclusive slice
                return full_seq[start-1:end]
            return full_seq # Fallback to full sequence if parsing fails
            
        elif target_type == 'gene':
            target = ALIAS_MAP.get(raw_symbol, raw_symbol)
            clean_target = target.lower().split('(')[0].strip()
            clean_raw = raw_symbol.lower().split('(')[0].strip()
            
            for feature in record.features:
                if feature.type == "CDS" or feature.type == "gene":
                    gene = feature.qualifiers.get("gene", [""])[0].lower()
                    product = feature.qualifiers.get("product", [""])[0].lower()
                    locus_tag = feature.qualifiers.get("locus_tag", [""])[0].lower()
                    
                    match = False
                    if clean_target == gene or clean_target == locus_tag:
                        match = True
                    elif clean_raw == gene or clean_raw == locus_tag:
                        match = True
                    elif clean_target in product or clean_raw in product:
                        match = True
                            
                    if match:
                        return str(feature.extract(record.seq)).upper()
            return full_seq # Fallback
            
    except Exception as e:
        logging.error(f"Error fetching {accession}: {e}")
    return ""

def main():
    logging.info("Downloading miRNA sequences...")
    mirna_dict = download_mirbase()
    
    logging.info("Loading dataset...")
    df = pd.read_csv(INPUT_CSV)
    
    # Filter to unique interactions to save API calls
    df = df.drop_duplicates(subset=['miRNA_ID', 'Target_Symbol', 'Taxonomy_ID']).copy()
    
    records = []
    genome_cache = {}
    seq_cache = {}
    
    for idx, row in df.iterrows():
        mirna = row['miRNA']
        mimat_id = str(row['miRNA_ID']).strip()
        virus_name = row['Virus_Name']
        tax_id = str(row['Taxonomy_ID']).strip()
        target_symbol = str(row['Target_Symbol']).strip()
        pmid = row['PMID']
        
        mirna_seq = mirna_dict.get(mimat_id, "")
        
        target_type = classify_target(target_symbol)
        
        accession = ""
        if tax_id and tax_id != 'nan':
            if tax_id not in genome_cache:
                genome_cache[tax_id] = get_genome_from_taxonomy(tax_id)
            accession = genome_cache[tax_id]
            
        target_seq = ""
        if accession:
            cache_key = (accession, target_symbol, target_type)
            if cache_key not in seq_cache:
                seq_cache[cache_key] = fetch_target_sequence(accession, target_symbol, target_type)
            target_seq = seq_cache[cache_key]
            
        records.append({
            "miRNA": mirna,
            "miRNA_sequence": mirna_seq,
            "virus": virus_name,
            "target_symbol": target_symbol,
            "target_type": target_type,
            "target_sequence": target_seq,
            "genome_accession": accession,
            "PMID": pmid
        })
        
        if (idx + 1) % 10 == 0:
            logging.info(f"Processed {idx + 1}/{len(df)} rows...")
            
    out_df = pd.DataFrame(records)
    out_df.to_csv(OUTPUT_CSV, index=False)
    
    logging.info("================ SUMMARY ================")
    logging.info(f"Total interactions processed: {len(out_df)}")
    logging.info(f"Rows with miRNA sequence:     {len(out_df[out_df['miRNA_sequence'] != ''])}")
    logging.info(f"Rows with Target sequence:    {len(out_df[out_df['target_sequence'] != ''])}")
    logging.info("=========================================")

if __name__ == "__main__":
    main()
