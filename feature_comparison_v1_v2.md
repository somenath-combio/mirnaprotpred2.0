# Feature Comparison: miRNAProtPred 1.0 vs. miRNAProtPred 2.0

This document provides a detailed technical comparison of the features, computations, and design differences between **miRNAProtPred 1.0** and **miRNAProtPred 2.0**. This summary is formatted for direct integration into academic publication Methods and Results sections.

---

## 1. Feature Representation Philosophy

*   **miRNAProtPred 1.0 (Heuristic/Thermodynamic Ranker)**:
    *   Relied on a hand-crafted heuristic scoring system (`Final_Score`).
    *   Scores were computed by combining ViennaRNA binding energy (`RNAduplex_MFE`) with discrete weights assigned to seed types, flanking AU content, and matching characteristics.
    *   Designed for ranking genomic candidate target sites (CTS) using static weights without context-aware modeling.

*   **miRNAProtPred 2.0 (Context-Aware ML Classifier)**:
    *   Transitioned to a machine learning classification approach using a **Random Forest (RF)** and **XGBoost (XGB)** ensemble.
    *   Eliminated hand-crafted weight formulas.
    *   Uses a **49-dimensional biological feature vector** containing sequence lengths, nucleotide compositions, local nucleotide frequencies (dinucleotides and trinucleotides), target genomic location annotations, and host viral taxonomy context.
    *   Predicts a continuous interaction probability ($P_{ML} \in [0, 1]$) calibrated into a 4-tier confidence grading system.

---

## 2. Detailed Feature Catalog and Computations

### A. Core Thermodynamic and Alignment Features

| Feature Name (v1) | Feature Name (v2) | Computation / Formula | Description |
| :--- | :--- | :--- | :--- |
| `RNAduplex_MFE` | *Dropped from ML* (Used in grading) | $\Delta G \text{ (kcal/mol)} = \text{duplexfold}(miRNA, CTS)$ | Minimum Free Energy (MFE) computed using ViennaRNA package (`RNA.duplexfold` function). Used in v2 for 4-tier confidence classification but excluded from ML training features to prevent learning bias. |
| `motif_identity` | `motif_identity` | $Identity = Seed\_Match \times Supp\_Match$ | In v1, computed as seed match fraction in relaxed mode. In v2, it represents the product of seed and supplementary match binary flags. |
| `motif_extension_score` | `supp_match` | $Matches / 5$ (for flanking nt 7–11/13–17) | Matches in the supplementary flanking region. v1 scanned nt 7–11; v2 targets nt 13–17 of the mature miRNA alignment. |
| `AU_Context` | `cts_au` | $\frac{A + U}{Length_{CTS}}$ | In v1, AU fraction was calculated on 20 nt upstream/downstream flanking regions. In v2, it is computed directly on the target candidate sequence (CTS) window. |
| `seed_col` | `seed_match` | Binary Flag (0 or 1) | Seed interaction matching representation (canonical positions 2–7). v1 mapped to categorical priorities (8mer, 7mer-m8, 7mer-A1), while v2 uses a binary flag. |

---

### B. New Sequence and Composition Features in miRNAProtPred 2.0

These features were introduced in version 2.0 to capture sequence-dependent binding signatures and local target structural accessibility.

*   **Sequence Lengths**:
    *   `cts_len`: Length of the target candidate sequence window (nucleotides).
    *   `mirna_len`: Length of the mature miRNA sequence (nucleotides).
*   **Nucleotide Composition**:
    *   `mirna_gc`: GC content of the mature miRNA:
        $$GC_{miRNA} = \frac{G + C}{Length_{miRNA}}$$
    *   `cts_gc`: GC content of the target candidate sequence window:
        $$GC_{CTS} = \frac{G + C}{Length_{CTS}}$$
*   **Dinucleotide Frequencies (16 Features)**:
    *   `dinuc_AA` to `dinuc_UU`: The relative frequency of adjacent nucleotide pairs in the CTS window (e.g., $AA, AC, \dots, UU$):
        $$Freq(dinuc\_XY) = \frac{\text{Count of } XY}{Length_{CTS} - 1}$$
*   **Trinucleotide Frequencies (8 Features)**:
    *   `trinuc_AAA`, `trinuc_UUU`, `trinuc_CCC`, `trinuc_GGG`, `trinuc_GCC`, `trinuc_CGC`, `trinuc_GCG`, `trinuc_CGU`: Relative frequencies of selected nucleotide triplets in the CTS window:
        $$Freq(trinuc\_XYZ) = \frac{\text{Count of } XYZ}{Length_{CTS} - 2}$$

---

### C. New Contextual and Taxonomic Features in miRNAProtPred 2.0

*   **Target Annotation Type (One-Hot Encoded)**:
    *   `targettype_gene`: Target sequence represents a coding sequence (CDS) or gene region.
    *   `targettype_transcript`: Target sequence represents a transcript/mRNA sequence.
    *   `targettype_region`: Target sequence represents an unannotated genomic region.
*   **Host Viral Taxonomic Family (One-Hot Encoded - 14 Columns)**:
    *   Indicator variables mapping the target genome to its taxonomic family:
        `virusfamily_Adenoviridae`, `virusfamily_Arteriviridae`, `virusfamily_Bornaviridae`, `virusfamily_Filoviridae`, `virusfamily_Flaviviridae`, `virusfamily_Hepadnaviridae`, `virusfamily_Herpesviridae`, `virusfamily_Orthomyxoviridae`, `virusfamily_Papillomaviridae`, `virusfamily_Picornaviridae`, `virusfamily_Pneumoviridae`, `virusfamily_Retroviridae`, `virusfamily_Rhabdoviridae`, `virusfamily_Togaviridae`.

---

## 3. Features Dropped in miRNAProtPred 2.0

The following features and intermediate variables from the v1.0 heuristic scoring engine were removed from the v2.0 ML training pipeline:

1.  **`competing_seeds`**: The density of alternative seed matches in the target sequence window. Dropped to prevent duplicate target site classification bias.
2.  **`match_type`**: Categorical match classification (e.g., 'Exact Match' vs. 'Wobble Pairing'). Replaced with quantitative features (`seed_match`, `motif_identity`).
3.  **`Adjusted_Energy`**: Hardcoded thermodynamic scaling factor. Replaced by a continuous probability calibration.
4.  **`Final_Score`**: Heuristic weighting sum:
    $$Final\_Score = |RNAduplex\_MFE| \times 0.35 + Context\_Score + Identity\_Score + \dots$$
    Dropped in favor of pure ML prediction probabilities.

---

## 4. Paper-Ready Methods Summary Table

The table below outlines the feature set transition from miRNAProtPred 1.0 to 2.0:

| Feature Dimension | Category | miRNAProtPred 1.0 (Heuristics) | miRNAProtPred 2.0 (ML Classifier) | Status |
| :--- | :--- | :--- | :--- | :--- |
| **Thermodynamic** | Hybrid Energy | `RNAduplex_MFE`, `Adjusted_Energy` | None (Decoupled to Confidence Grading) | **Dropped / Moved** |
| **Length** | Sequence Length | None | `cts_len`, `mirna_len` | **New** |
| **Nucleotide Bias** | Global Composition | `AU_Context` (flanking only) | `mirna_gc`, `cts_gc`, `cts_au` | **Updated / New** |
| **Seed Match** | Structural Match | `seed_col` (Categorical), `match_type` | `seed_match`, `supp_match`, `motif_identity` | **Updated** |
| **K-Mer Signatures** | Dinucleotides | None | 16 features (`dinuc_AA` to `dinuc_UU`) | **New** |
| **K-Mer Signatures** | Trinucleotides | None | 8 features (`trinuc_AAA` to `trinuc_CGU`) | **New** |
| **Annotation Context**| Substrate Location | None | `targettype_gene`, `targettype_transcript`, `targettype_region` | **New** |
| **Taxonomy Context** | Viral Host | None | 14 taxonomic family indicator columns | **New** |
| **Multi-Site** | Target Density | `competing_seeds` | None | **Dropped** |
| **Total Features** | - | **6 features** | **49 features** | **Expansion** |
