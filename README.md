# miRNAProtPred 2.0

Two-stage viral miRNA–CTS prediction and validation framework.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

miRNAProtPred 2.0 is a next-generation bioinformatics framework for scanning and validating microRNA (miRNA) interactions against viral genomes. It upgrades the first-generation pipeline to a robust two-stage architecture: **Stage 1 (IntaRNA thermodynamic scanning)** followed by **Stage 2 (RandomForest/XGBoost Machine Learning context classification)**.

---

## Comparison: miRNAProtPred 1.0 vs 2.0

| Feature | miRNAProtPred 1.0 | miRNAProtPred 2.0 |
| :--- | :--- | :--- |
| **CTS Detection** | Boyer-Moore exact seed matching | IntaRNA thermodynamic genome scanning |
| **Scoring / Verification** | Custom alignment MFE classification | Multi-model Machine Learning classifier (RandomForest/XGBoost) |
| **Seed Requirement** | Strictly required (canonical matches only) | Optional (supports non-canonical/wobble targets) |
| **Novel Virus Support** | Re-training required for species specificity | Zero-shot generalization via 4-tier confidence model |
| **Multi-substrate Inputs** | DNA, RNA, Protein (remote tblastn) | DNA, RNA, Protein (remote tblastn - Retained & Optimized) |
| **Tool Division** | `SeqFinder` (scan) & `validator` (target) | `SeqFinder2` (scan), `validator2` (target), `train_model2` (ML training) |

---

## Key Features

- **Multi-Substrate Input Recognition**: Natively accepts DNA, RNA, or Protein sequences (either as direct FASTA files or nucleotide strings). Protein inputs automatically trigger a remote NCBI `tblastn` alignment to resolve and download the source coding sequence region.
- **Strict & Relaxed Search Modes**:
  - `strict`: Enforces canonical 7-mer seed matches (miRNA positions 2-8) for fast filtering.
  - `relaxed`: Allows G·U wobble pairings and non-canonical bindings.
- **High-Performance Whole-Database Scanning**: Scan a viral genome against the entire database of 2,656 human miRNAs in seconds using Boyer-Moore seed complement pre-filtering.
- **Premium Machine Learning Scoring**: Predicts interaction validity based on GC/AU content, dinucleotide/trinucleotide distributions, and thermodynamic parameters.
- **Out-of-the-Box Serialization**: Packaged with pre-trained RandomForest classifiers, enabling instant prediction without local re-training.

---

## Installation

### Prerequisites

- **Python 3.10**
- **IntaRNA** (Must be installed via Conda)

### 1. Set Up Conda Environment

```bash
# Create and activate environment
conda create -n mirnaprotpred2_env python=3.10 -y
conda activate mirnaprotpred2_env

# Install IntaRNA from Bioconda
conda install -c bioconda intarna -y
```

### 2. Install Package

```bash
git clone https://github.com/somenath-combio/mirnaprotpred2.0.git
cd mirnaprotpred2.0
pip install -e .
```

---

## CLI Reference

The framework installs three executable CLI entry points:

### 1. `SeqFinder2` (Genome Target Scanner)

Scans a target viral genome or sequence to discover potential miRNA target sites (CTS).

#### Syntax
```bash
SeqFinder2 <target_sequence_or_fasta> [options]
```

#### Arguments & Options
*   `genome` (Positional): Path to the target viral FASTA file or direct nucleotide/protein sequence string.
*   `--mirna` (Optional, default: `ALL`):
    *   `ALL`: Scans the genome against the entire database of 2,656 human miRNAs.
    *   `path/to/miRNAs.fasta`: Scans the genome against all miRNA sequences in the FASTA file.
    *   `UAGCUUAUCAGACUGAUGUUGA`: Scans the genome for a single mature miRNA sequence.
*   `--mode` (Optional, default: `strict`):
    *   `strict`: Enforces strict 7-mer seed matches.
    *   `relaxed`: Catches wobble pairings and non-canonical matches.
*   `--output` (Optional, default: `concise`):
    *   `concise`: Prints target sites of `High` or `Very High` confidence in a compact layout.
    *   `raw`: Prints all discovered candidates with full biophysical details.
    *   `highconf`: Prints concise details filtered to `Very High`, `High`, and `Medium` confidence.
    *   `highconf_raw`: Prints all columns filtered to `Very High`, `High`, and `Medium` confidence.
*   `--out` (Optional): Saves the display table to a CSV file.
*   `--win` (Optional, default: `50`): Sliding window size.
*   `--step` (Optional, default: `10`): Step size.
*   `--dg` (Optional, default: `-5.0`): delta_G threshold.
*   `--email` (Optional, default: `sudipta@pusan.ac.kr`): Entrez email for remote BLAST reverse-mapping.

#### Examples
```bash
# Scan a genome against the entire database in concise mode
SeqFinder2 sars_cov2.fasta

# Scan a genome in relaxed mode and output raw candidates to CSV
SeqFinder2 sars_cov2.fasta --mode relaxed --output raw --out results.csv
```

---

### 2. `validator2` (Targeted miRNA Verification)

Skip whole-genome databases to verify specific, known miRNA interactions.

#### Syntax
```bash
validator2 <miRNA_IDs> <target_sequence_or_fasta> [options]
```

#### Arguments & Options
*   `mirna` (Positional 1): miRNA IDs (comma-separated list, text file path with one ID per line, or FASTA file path).
*   `genome` (Positional 2): Path to target genome FASTA or raw sequence string.
*   `--mode` (Optional, default: `strict`):
    *   `strict`: Rapid Boyer-Moore seed search validation.
    *   `relaxed`: Sequential sliding-window scan of each target miRNA, allowing wobble and non-canonical binding validation.
*   `--output` (Optional, default: `summary`):
    *   `summary`: Prints a clean YES/NO verification mapping table for the input list.
    *   `raw` (or `--details` alias): Expands to display full candidate parameters for matched sites.
    *   `highconf`: Displays detailed candidate parameters filtered to `Very High`, `High`, and `Medium` confidence.
*   `--out` (Optional): Saves results to CSV.

#### Examples
```bash
# Verify specific miRNA list against HIV genome
validator2 hsa-miR-29a-3p,hsa-miR-3941 hiv.fasta

# Verify using a candidate list text file, outputting raw results to CSV
validator2 candidates.txt hiv.fasta --output raw --out validation.csv
```

---

### 3. `train_model2` (ML Training Engine)

Re-train the model or run local grid search and generalization checks on your datasets.

#### Syntax
```bash
train_model2 [options]
```
Outputs model checkpoint files (`mirnaprotpred2_best.pkl`, `mirnaprotpred2_xgb.pkl`) and validation reports locally.

---

## 4-Tier Confidence Grading

Every candidate interaction is automatically classified into a confidence tier using both biological machine learning scores and thermodynamic stability:

| Tier | Condition | Recommendation |
| :--- | :--- | :--- |
| **Very High** | `ml_score` >= 0.75 AND `delta_G` <= -10.0 | High-confidence binding, optimal for wet lab assay validation. |
| **High** | `ml_score` >= 0.50 | Corresponds to established experimental VIRBase targets. |
| **Medium** | `ml_score` >= 0.25 AND `delta_G` <= -9.0 | Thermodynamically stable with minor contextual mismatching. |
| **Low** | All other matches | Exploratory predictions only. |

---

## Authors

- **Somenath Dutta** (somenath@pusan.ac.kr)
- **Sudipta Sardar** (sudipta@pusan.ac.kr)

---

## License

MIT License
