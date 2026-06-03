# 📋 Pipeline Step Guide: miRNAProtPred 2.0 Reproducibility

This folder contains the complete, sequential pipeline scripts used to build, train, benchmark, and visualize the **miRNAProtPred 2.0** framework from scratch.

---

## 🛠️ Pipeline Steps Summary

| Step | Script Name | Purpose | Primary Inputs | Primary Outputs |
| :---: | :--- | :--- | :--- | :--- |
| **01** | `step01_filter_interactions.py` | Filters human miRNA-viral mRNA interactions. | Raw DB files | `virmirna_only.csv` |
| **02** | `step02_retrieve_sequences.py` | Retrieves target sequences from NCBI RefSeq. | Accession numbers | NCBI API downloads |
| **03** | `step03_sequence_qc.py` | Performs quality control on retrieved sequences. | Raw FASTA sequences | Verified sequences |
| **04** | `step04_recovery_audit.py` | Audits sequence recovery rate and missing records. | Retrieval logs | `recovery_summary_report.csv` |
| **05** | `step05_update_metadata.py` | Merges salvaged targets and populates taxonomies. | Verified records | `viral_miRNA_ML_dataset.csv` |
| **06** | `step06_final_qc.py` | Validates dataset integrity and duplicate removal. | Main dataset | Cleaned dataset |
| **07** | `step07_generate_cts.py` | Scans for thermodynamic Candidate Target Sites (CTS). | Cleaned sequences | `cts_positives.csv` |
| **08** | `step08_clean_cts.py` | Filters CTS and generates paired hard negatives. | `cts_positives.csv` | `mirnaprotpred2.0_training_set.csv` |
| **09** | `step09_feature_extraction.py` | Computes the 49-dimensional biophysical features. | Training set CSV | `cts_ml_features_v4_intraviral.csv` |
| **10** | `step10_model_training.py` | Trains RF and XGBoost classifiers via CV. | Features CSV | `mirnaprotpred2.0_xgb.pkl`, `mirnaprotpred2.0_best.pkl` |
| **11** | `step11_prepare_miranda.py` | Prepares the 170-sample independent validation set. | `virmirna_85_normalized.csv` | `mirnas_85.fa`, `targets_85.fa` |
| **12** | `step12_targetscan_benchmark.py`| Evaluates TargetScan seed-match rules on test set. | `virmirna_85_with_targets.csv` | TargetScan AUC, F1 metrics |
| **13** | `step13_miranda_benchmark.py` | Parses and evaluates miRanda scan results. | `miranda_output_85.txt.gz` | miRanda AUC, F1 metrics |
| **14** | `step14_generate_figures.py` | Generates all 6 publication figures at 600 DPI. | Model pkls, Test sets | Fig 1 - Fig 6 (.png files) |

---

## 🔍 Detailed Step Explanations

### **Step 1: Filter Interactions** (`step01_filter_interactions.py`)
Parses raw database exports (from VIRmiRNA/ViRBase), normalizes miRNA names, filters out host gene targets to prevent contamination, and saves the verified human miRNA to viral mRNA interactions list.

### **Step 2: Retrieve Sequences** (`step02_retrieve_sequences.py`)
Automates queries to NCBI's Entrez E-utilities API using viral taxonomy IDs and NC/NM accession numbers to retrieve full-length target viral genome/transcript sequences.

### **Step 3: Sequence QC** (`step03_sequence_qc.py`)
Ensures all retrieved target sequences contain valid IUPAC nucleotide characters, checks sequence lengths, and flags truncated or invalid records.

### **Step 4: Recovery Audit** (`step04_recovery_audit.py`)
Categorizes sequence downloads, separating successfully recovered genes from those with obsolete or missing records, producing a list of targets requiring manual taxonomy-based recovery.

### **Step 5: Update Metadata** (`step05_update_metadata.py`)
Merges manually salvaged and NCBI-downloaded sequences, maps viral family taxonomic groupings, and registers coordinates/target regions.

### **Step 6: Final QC** (`step06_final_qc.py`)
A final sanity check that performs sequence deduplication, validates miRNA-target pair mappings, and ensures data integrity.

### **Step 7: Generate CTS** (`step07_generate_cts.py`)
Locates the 50-nucleotide Candidate Target Site (CTS) representing the most thermodynamically stable binding pocket for each miRNA on its viral target.

### **Step 8: Clean CTS & Generate Negatives** (`step08_clean_cts.py`)
Extracts positive CTS binding sites and shuffles their target sequences using random permutation, filtering for non-specific binding (IntaRNA MFE > -2.0 kcal/mol, seed match < 0.6) to generate biologically rigorous **hard negatives**.

### **Step 9: Feature Extraction** (`step09_feature_extraction.py`)
Computes the 49 biophysical, structural, and sequence-based features (including seed-match parameters, MFE, GC content, and nucleotide frequency patterns).

### **Step 10: Model Training** (`step10_model_training.py`)
Performs hyperparameter tuning and cross-validation to train the final RandomForest (best baseline) and XGBoost (best ensemble classifier) models, saving the frozen binary models.

### **Step 11: Prepare miRanda Inputs** (`step11_prepare_miranda.py`)
Prepares the independent 170-sample validation dataset (85 positive, 85 shuffled negatives) and outputs unique indexed FASTA files for baseline scans.

### **Step 12: TargetScan Benchmark** (`step12_targetscan_benchmark.py`)
Implements TargetScan's canonical seed-matching rules (8mer, 7mer-m8, 7mer-A1, 6mer) directly in Python and evaluates its classification performance.

### **Step 13: miRanda Benchmark** (`step13_miranda_benchmark.py`)
Parses alignment scores and minimum free energy outputs from the miRanda scan on the validation dataset to compute performance metrics.

### **Step 14: Generate Figures** (`step14_generate_figures.py`)
Generates the publication-grade figures (600 DPI) mapping performance comparisons (grouped bar charts, ROC/PR curves, feature importances, and distribution graphs) to visual files in `publication_figures/`.
