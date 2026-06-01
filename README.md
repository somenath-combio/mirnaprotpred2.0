# miRNAProtPred 2.0

**Two-stage viral miRNA–CTS prediction pipeline.**
Upgrade from v1 (Boyer-Moore + BLAST) → v2 (IntaRNA thermodynamics + XGBoost ML).

## vs miRNAProtPred 1.0
| Feature              | v1                      | v2                          |
|----------------------|-------------------------|-----------------------------|
| CTS detection        | Boyer-Moore seed match  | IntaRNA thermodynamic scan  |
| Scoring              | BLAST alignment         | XGBoost ML classifier       |
| Seed required        | Yes (strict)            | No (non-canonical supported)|
| Novel virus support  | Retrain needed          | 3-tier confidence (no retrain)|
| Validated pairs      | VIRBase curated         | VIRBase + VIRmiRNA external |

## Confidence Tiers
| Tier   | Condition                       | Use                          |
|--------|---------------------------------|------------------------------|
| High   | prob ≥ 0.50                     | Experimental-grade evidence  |
| Medium | prob ≥ 0.25 AND dG ≤ −9.0      | Prioritize for wet lab       |
| Low    | everything else                 | Computational only           |

## Key Results
- Training F1        : 1.000 (VIRBase, 462 pairs)
- External AUC       : 0.83  (VIRmiRNA 85-pair, zero-shot)
- Shuffle test       : 86.7% significant (p<0.05)
- SARS-CoV-2 recovery: 16/16 miRNAs (SeqFinder Stage 1)
- Novel virus        : hsa-miR-3941 → 3'UTR recovered as Medium confidence

## Author
Somenath Dutta (somenath@pusan.ac.kr)
