# miRNAProtPred 2.0

Two-stage viral miRNA–CTS prediction and validation framework.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

miRNAProtPred 2.0 is a next-generation bioinformatics framework for scanning and validating microRNA (miRNA) interactions against viral genomes. It upgrades the first-generation pipeline (which relied on simple Boyer-Moore string searches and BLAST) to a more robust two-stage architecture: **Stage 1 (IntaRNA thermodynamic scanning)** followed by **Stage 2 (RandomForest/XGBoost Machine Learning context classification)**.

---

## Comparison: miRNAProtPred 1.0 vs 2.0

| Feature | miRNAProtPred 1.0 | miRNAProtPred 2.0 |
| :--- | :--- | :--- |
| **CTS Detection** | Boyer-Moore exact seed matching | IntaRNA thermodynamic genome scanning |
| **Scoring / Verification** | Custom alignment MFE classification | Multi-model Machine Learning classifier (RandomForest/XGBoost) |
| **Seed Requirement** | Strictly required (canonical matches only) | Optional (supports non-canonical/wobble targets) |
| **Novel Virus Support** | Re-training required for species specificity | Zero-shot generalization via 3-tier confidence model |
| **Validation Dataset** | VIRBase curated interactions | Multi-source (VIRBase + external VIRmiRNA benchmarking) |

---

## Key Features

- **High-Performance Whole-Database Scanning**: Scan a viral genome against the entire database of 2,656 human miRNAs in seconds using an optimized Boyer-Moore seed complement pre-filter.
- **Custom miRNA File Input**: Scan genomes using a custom FASTA database of miRNAs or a single target miRNA sequence.
- **Advanced Thermodynamics**: Integrates `IntaRNA` directly to compute minimum free energy (MFE) structural binding dynamics.
- **Premium Machine Learning scoring**: Incorporates a trained Random Forest model predicting interaction validity based on GC content, nucleotide distribution, and thermodynamic scores.
- **Flexible CLI Pipelines**: Custom entry points `SeqFinder2` and `validator2` configured out-of-the-box.

---

## Installation

### Prerequisites

- **Python 3.8 or higher** (Python 3.10 recommended)
- **IntaRNA** (Must be installed via Conda)

### Install via Conda (Recommended)

To ensure all compiled dependencies (especially `IntaRNA` from Bioconda) resolve correctly, run:

```bash
# Create and activate environment
conda create -n mirnaprotpred2_env python=3.10
conda activate mirnaprotpred2_env

# Install IntaRNA via Bioconda
conda install -c bioconda intarna
```

### Install the Package

Ensure your environment is active, then clone and install the repository in editable mode:

```bash
git clone https://github.com/somenath-combio/mirnaprotpred2.0.git
cd mirnaprotpred2.0
pip install -e .
```

---

## Command Line Interface (CLI) Reference

The framework provides two primary entry points:

### 1. `SeqFinder2` (Thermodynamic & ML Target Scanner)

Scans a target viral genome/sequence to discover potential miRNA target sites (CTS).

```bash
SeqFinder2 --genome <viral_genome.fasta> [options]
```

#### Inputs
*   `--genome` (Required): Path to the target viral mRNA or genome FASTA file.
*   `--mirna` (Optional, default: `ALL`):
    *   `ALL`: Scans the genome against the default database of 2,656 human miRNAs (`data.xlsx`).
    *   `path/to/miRNAs.fasta`: Scans the genome against all miRNA sequences in the FASTA file.
    *   `UAGCUUAUCAGACUGAUGUUGA`: Scans the genome for a single inline miRNA sequence.
*   `--win` (Optional, default: `50`): Sliding target window length in nt.
*   `--step` (Optional, default: `10`): Step size for the sliding window in nt.
*   `--dg` (Optional, default: `-5.0`): delta_G thermodynamic threshold (retains interactions with MFE <= threshold).
*   `--top` (Optional, default: `50`): Returns the top N candidates.
*   `--seed` (Optional): If set, enforces a strict 7-mer seed match (positions 2-8 of the miRNA) during the scan.
*   `--out` (Optional, default: `seqfinder_output.csv`): Path to output results CSV.
*   `--model` (Optional): Path to a custom serialized model `.pkl`.
*   `--traincsv` (Optional): Path to training CSV (used for alignment columns).

#### Outputs
Saves a structured CSV file containing candidate interactions. Columns include:
- `miRNA_ID`: The human miRNA identifier (or custom ID).
- `miRNA_seq`: The sequence of the mature miRNA.
- `start`: 0-based start position of the target window in the genome.
- `end`: 0-based end position of the target window in the genome.
- `window_seq`: Sequence of the target site (Complementary Target Sequence).
- `delta_G`: Thermodynamic Minimum Free Energy (MFE) calculated by IntaRNA.
- `hybridDP`: Dot-bracket secondary structure notation of the binding hybrid.
- `n_base_pairs`: Number of base pairs in the hybrid structure.
- `seed_match`: 1 if a 7-mer seed complement is found, 0 otherwise.
- `supp_match`: Complementary match score at positions 13-17 of the miRNA.
- `cts_gc`: GC content of the target site.
- `ml_score`: The probability score (from 0.0 to 1.0) assigned by the ML model.

#### Examples
```bash
# Scan a genome against the entire human miRNA database
SeqFinder2 --genome sars_cov2.fasta --out sars_results.csv

# Scan a genome against a custom miRNA FASTA file
SeqFinder2 --mirna my_mirnas.fa --genome sars_cov2.fasta --dg -7.0

# Scan a single miRNA with strict seed matching
SeqFinder2 --mirna UAGCUUAUCAGACUGAUGUUGA --genome sars_cov2.fasta --seed
```

---

### 2. `validator2` (ML Training & LOVO Cross-Validation Pipeline)

Trains, validates, and serializes the machine learning classifiers on the full dataset, performing hyperparameter tuning and a Leave-One-Virus-Out (LOVO) species generalization check.

```bash
validator2 [options]
```

#### Inputs
*   `-h`, `--help`: Prints the help message and exits.

#### Outputs
Executes grid search training and saves results to the local directory (or `output/`):
- `mirnaprotpred2_best.pkl`: Serialized final production classifier (Random Forest).
- `mirnaprotpred2_xgb.pkl`: Companion XGBoost model for SHAP explainability.
- `cts_cv_predictions_multi.csv`: Out-of-fold cross-validation predictions for all models.
- `virus_generalization_report.csv`: Generalization report (species-level AUC, Accuracy) across all 24 viruses.

---

## Confidence Tiers

Interactions identified by `SeqFinder2` are categorized into three confidence tiers based on the ML scoring output and thermodynamic stability:

| Tier | Condition | Use Case |
| :--- | :--- | :--- |
| **High** | `ml_score` >= 0.50 | Recommended for immediate experimental validation (wet lab assays). |
| **Medium**| `ml_score` >= 0.25 AND `delta_G` <= -9.0 | High thermodynamic stability with supporting contextual features. |
| **Low**  | All other matches | Exploratory computational predictions only. |

---

## How It Works

1.  **Thermodynamic Screening (Stage 1)**:
    - If scanning multiple miRNAs, the engine uses a Boyer-Moore seed search to locate seed complements (positions 2-8) in the genome.
    - `IntaRNA` is invoked on target windows (default: 50 nt) surrounding the seed match to calculate MFE (`delta_G`) and base-pair interactions.
2.  **Machine Learning Classification (Stage 2)**:
    - Candidate structures are parsed into dynamic features including GC/AU content, dinucleotide and trinucleotide frequencies, and seed/supplementary matching flags.
    - The serialized Random Forest classifier uses these features to output the interaction probability score (`ml_score`), filtering out false positive interactions.

---

## Authors

- **Somenath Dutta** (somenath@pusan.ac.kr)
- **Sudipta Sardar** (sudipta@pusan.ac.kr)

---

## License

MIT License
