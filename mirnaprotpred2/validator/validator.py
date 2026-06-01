#!/usr/bin/env python3
"""
train_miRNAProtPred2.py
────────────────────────────────────────────────────────────────────────────
Machine Learning Training and Validation Pipeline for miRNAProtPred 2.0.
1. Loads the 67-dimension feature vector from main_work/cts_ml_features.csv.
2. Conducts grid searches and Stratified 5-Fold CV for 5 models:
   - Logistic Regression
   - Random Forest
   - XGBoost
   - LightGBM
   - CatBoost
3. Identifies the best-performing model based on cross-validation ROC-AUC.
4. Performs Leave-One-Virus-Out (LOVO) cross-validation on all 24 viruses.
5. Saves all out-of-fold prediction probabilities and generalization reports.
6. Serializes the final optimal classifier object on the full dataset.
────────────────────────────────────────────────────────────────────────────
"""

import pandas as pd
import numpy as np
import pickle
import os
from pathlib import Path
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    matthews_corrcoef,
    f1_score,
    accuracy_score
)

# Import models
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
from catboost import CatBoostClassifier

# Suppress warnings
import warnings
warnings.filterwarnings("ignore")

ROOT = Path("/home/somenath/Pictures/Publication_somenath")
INPUT_FEATURES_CSV = ROOT / "virbase_final_dataset/virbase_cts/cts_ml_features_v4_intraviral.csv"
OUTPUT_PREDS_CSV      = ROOT / "output/cts_cv_predictions_multi.csv"
OUTPUT_BEST_MODEL_PKL = ROOT / "output/mirnaprotpred2_best.pkl"
OUTPUT_XGB_MODEL_PKL  = ROOT / "output/mirnaprotpred2_xgb.pkl"
OUTPUT_LOVO_REPORT_CSV= ROOT / "output/virus_generalization_report.csv"

# Taxonomic Virus Family lookup for LOVO reporting
VIRUS_TO_FAMILY = {
    'Human mastadenovirus C': 'Adenoviridae',
    'Porcine reproductive and respiratory syndrome virus (PRRSV)': 'Arteriviridae',
    'Borna disease virus (BDV)': 'Bornaviridae',
    'Zaire ebolavirus (ZEBOV)': 'Filoviridae',
    'Hepacivirus C (HCV)': 'Flaviviridae',
    'Zika virus (ZIKV)': 'Flaviviridae',
    'Hepatitis B virus (HBV)': 'Hepadnaviridae',
    'Human betaherpesvirus 5 (HHV-5)': 'Herpesviridae',
    'Human gammaherpesvirus 4 (Epstein-Barr virus, EPV)': 'Herpesviridae',
    "Human gammaherpesvirus 8 (Kaposi's sarcoma herpesvirus, KSHV)": 'Herpesviridae',
    'Influenza A virus (H1N1)': 'Orthomyxoviridae',
    'Influenza A virus (H3N2)': 'Orthomyxoviridae',
    'Influenza A virus (H5N1)': 'Orthomyxoviridae',
    'Human papillomavirus type 16 (HPV16)': 'Papillomaviridae',
    'Human papillomavirus type 18': 'Papillomaviridae',
    'Coxsackievirus B3 (CVB3)': 'Picornaviridae',
    'Enterovirus A71 (EV-A71)': 'Picornaviridae',
    'Human orthopneumovirus (HRSV)': 'Pneumoviridae',
    'Human immunodeficiency virus 1 (HIV-1)': 'Retroviridae',
    'Human T-lymphotropic virus 1 (HTLV-1)': 'Retroviridae',
    'Simian immunodeficiency virus (SIV)': 'Retroviridae',
    'Vesicular stomatitis Indiana virus (VSIV)': 'Rhabdoviridae',
    'Eastern equine encephalitis virus (EEEV)': 'Togaviridae',
    'Sindbis virus': 'Togaviridae'
}

def main():
    print("Loading extracted feature dataset...")
    df = pd.read_csv(INPUT_FEATURES_CSV)
    
    # ─── DATA PREPARATION ───
    meta_cols = ["label", "miRNA", "virus", "target_symbol"]
    leaky_cols = [
        "delta_G", "delta_G_norm", "n_base_pairs",
        "struct_matches", "struct_wobbles", "struct_mismatches",
        "struct_bulges", "struct_loops", "site_position_norm"
    ]
    feature_cols = [c for c in df.columns if c not in meta_cols and c not in leaky_cols]
    
    print(f"Total feature dimensions: {len(feature_cols)}")
    
    X = df[feature_cols].copy()
    y = df["label"].copy()
    
    # Stratified 5-Fold CV Setup
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    # Define models and their grids
    models_config = {
        "LogisticRegression": {
            "model": LogisticRegression(random_state=42, solver='liblinear'),
            "grid": {
                "C": [0.01, 0.1, 1.0, 10.0],
                "penalty": ["l1", "l2"]
            }
        },
        "RandomForest": {
            "model": RandomForestClassifier(random_state=42),
            "grid": {
                "n_estimators": [50, 100, 200],
                "max_depth": [3, 5, 8, None],
                "min_samples_split": [2, 5]
            }
        },
        "XGBoost": {
            "model": XGBClassifier(random_state=42, eval_metric='logloss', use_label_encoder=False),
            "grid": {
                "n_estimators": [50, 100, 200],
                "learning_rate": [0.01, 0.05, 0.1],
                "max_depth": [3, 4, 5]
            }
        },
        "LightGBM": {
            "model": LGBMClassifier(random_state=42, verbose=-1),
            "grid": {
                "n_estimators": [50, 100, 200],
                "learning_rate": [0.01, 0.05, 0.1],
                "max_depth": [3, 4, 5]
            }
        },
        "CatBoost": {
            "model": CatBoostClassifier(random_state=42, verbose=0),
            "grid": {
                "iterations": [50, 100, 200],
                "learning_rate": [0.01, 0.05, 0.1],
                "depth": [3, 4, 5]
            }
        }
    }
    
    print("\n" + "=" * 65)
    print("  EXECUTING MULTI-MODEL GRID SEARCH & BENCHMARK")
    print("=" * 65)
    
    best_cv_models = {}
    best_cv_scores = {}
    oof_predictions_dict = {}
    oof_probabilities_dict = {}
    
    for name, config in models_config.items():
        print(f"\n[Training {name}] Running hyperparameter Grid Search...")
        grid_search = GridSearchCV(
            estimator=config["model"],
            param_grid=config["grid"],
            scoring='roc_auc',
            cv=skf,
            n_jobs=1,
            verbose=0
        )
        grid_search.fit(X, y)
        
        best_params = grid_search.best_params_
        best_score = grid_search.best_score_
        best_cv_scores[name] = best_score
        
        print(f"  - Best Params: {best_params}")
        print(f"  - Best CV ROC-AUC: {best_score:.4f}")
        
        # Instantiate best estimator with optimal parameters
        if name == "LogisticRegression":
            best_clf = LogisticRegression(**best_params, solver='liblinear', random_state=42)
        elif name == "RandomForest":
            best_clf = RandomForestClassifier(**best_params, random_state=42)
        elif name == "XGBoost":
            best_clf = XGBClassifier(**best_params, random_state=42, eval_metric='logloss', use_label_encoder=False)
        elif name == "LightGBM":
            best_clf = LGBMClassifier(**best_params, random_state=42, verbose=-1)
        elif name == "CatBoost":
            best_clf = CatBoostClassifier(**best_params, random_state=42, verbose=0)
            
        best_cv_models[name] = (best_clf, best_params)
        
        # Extract OOF Predictions
        oof_preds = np.zeros(len(df))
        oof_probs = np.zeros(len(df))
        
        for fold, (train_idx, val_idx) in enumerate(skf.split(X, y)):
            X_train, y_train = X.iloc[train_idx], y.iloc[train_idx]
            X_val, y_val = X.iloc[val_idx], y.iloc[val_idx]
            
            if name == "LogisticRegression":
                fold_clf = LogisticRegression(**best_params, solver='liblinear', random_state=42)
            elif name == "RandomForest":
                fold_clf = RandomForestClassifier(**best_params, random_state=42)
            elif name == "XGBoost":
                fold_clf = XGBClassifier(**best_params, random_state=42, eval_metric='logloss', use_label_encoder=False)
            elif name == "LightGBM":
                fold_clf = LGBMClassifier(**best_params, random_state=42, verbose=-1)
            elif name == "CatBoost":
                fold_clf = CatBoostClassifier(**best_params, random_state=42, verbose=0)
                
            fold_clf.fit(X_train, y_train)
            oof_preds[val_idx] = fold_clf.predict(X_val)
            oof_probs[val_idx] = fold_clf.predict_proba(X_val)[:, 1]
            
        oof_predictions_dict[name] = oof_preds
        oof_probabilities_dict[name] = oof_probs
        
        # Compute summary metrics for this model
        auc = roc_auc_score(y, oof_probs)
        pr_auc = average_precision_score(y, oof_probs)
        mcc = matthews_corrcoef(y, oof_preds)
        f1 = f1_score(y, oof_preds)
        acc = accuracy_score(y, oof_preds)
        
        print(f"  - OOF Metrics: AUC={auc:.4f} | PR-AUC={pr_auc:.4f} | MCC={mcc:.4f} | F1={f1:.4f} | ACC={acc:.4f}")
        
    # Determine overall best model based on CV ROC-AUC
    best_model_name = max(best_cv_scores, key=best_cv_scores.get)
    print("\n" + "=" * 65)
    print(f"  OVERALL BEST PERFORMING MODEL: {best_model_name}")
    print("=" * 65)
    
    # Save multi-model out-of-fold predictions to CSV
    pred_df = df[meta_cols].copy()
    for name in models_config.keys():
        pred_df[f"{name}_predicted_label"] = oof_predictions_dict[name].astype(int)
        pred_df[f"{name}_prediction_probability"] = oof_probabilities_dict[name]
        
    # For backwards-compatibility with single-model evaluation scripts
    pred_df["cv_predicted_label"] = oof_predictions_dict["XGBoost"].astype(int)
    pred_df["cv_prediction_probability"] = oof_probabilities_dict["XGBoost"]
    
    pred_df.to_csv(OUTPUT_PREDS_CSV, index=False)
    print(f"Saved all out-of-fold model predictions to: {OUTPUT_PREDS_CSV}")
    
    # ─── LEAVE-ONE-VIRUS-OUT (LOVO) CROSS-VALIDATION ───
    print("\n" + "=" * 65)
    print("  RUNNING LEAVE-ONE-VIRUS-OUT (LOVO) GENERALIZATION BENCHMARK")
    print("=" * 65)
    
    unique_viruses = df["virus"].unique()
    best_clf_template, best_params_template = best_cv_models[best_model_name]
    
    lovo_records = []
    
    # Global LOVO predictions for micro-average metrics
    global_lovo_preds = np.zeros(len(df))
    global_lovo_probs = np.zeros(len(df))
    
    n_viruses = len(unique_viruses)
    for v_idx, virus in enumerate(unique_viruses):
        test_mask = df["virus"] == virus
        train_mask = ~test_mask
        
        X_train, y_train = X[train_mask], y[train_mask]
        X_test, y_test = X[test_mask], y[test_mask]
        
        pos_test = (y_test == 1).sum()
        neg_test = (y_test == 0).sum()
        total_test = len(y_test)
        
        # Instantiate and train model
        if best_model_name == "LogisticRegression":
            lovo_clf = LogisticRegression(**best_params_template, solver='liblinear', random_state=42)
        elif best_model_name == "RandomForest":
            lovo_clf = RandomForestClassifier(**best_params_template, random_state=42)
        elif best_model_name == "XGBoost":
            lovo_clf = XGBClassifier(**best_params_template, random_state=42, eval_metric='logloss', use_label_encoder=False)
        elif best_model_name == "LightGBM":
            lovo_clf = LGBMClassifier(**best_params_template, random_state=42, verbose=-1)
        elif best_model_name == "CatBoost":
            lovo_clf = CatBoostClassifier(**best_params_template, random_state=42, verbose=0)
            
        lovo_clf.fit(X_train, y_train)
        
        # Predict on left-out virus
        test_preds = lovo_clf.predict(X_test)
        test_probs = lovo_clf.predict_proba(X_test)[:, 1]
        
        # Store in global array
        global_lovo_preds[test_mask] = test_preds
        global_lovo_probs[test_mask] = test_probs
        
        # Calculate virus-specific metrics
        # Standardize single-class metrics
        if pos_test > 0 and neg_test > 0:
            v_auc = roc_auc_score(y_test, test_probs)
            v_pr_auc = average_precision_score(y_test, test_probs)
            v_mcc = matthews_corrcoef(y_test, test_preds)
            v_f1 = f1_score(y_test, test_preds)
        else:
            # Single class present - AUC, PR-AUC, F1, MCC are mathematically undefined
            v_auc = np.nan
            v_pr_auc = np.nan
            v_mcc = np.nan
            v_f1 = np.nan
            
        v_acc = accuracy_score(y_test, test_preds)
        family = VIRUS_TO_FAMILY.get(virus, "Unknown_Family")
        
        print(f"  [{v_idx+1:2d}/{n_viruses}] {virus:<40} (Pos: {pos_test:2d}, Neg: {neg_test:2d}) | AUC: {v_auc:.4f} | ACC: {v_acc:.4f}")
        
        lovo_records.append({
            "Virus": virus,
            "Virus_family": family,
            "Positives": pos_test,
            "Negatives": neg_test,
            "ROC-AUC": round(float(v_auc), 5) if not np.isnan(v_auc) else "N/A",
            "PR-AUC": round(float(v_pr_auc), 5) if not np.isnan(v_pr_auc) else "N/A",
            "MCC": round(float(v_mcc), 5) if not np.isnan(v_mcc) else "N/A",
            "F1": round(float(v_f1), 5) if not np.isnan(v_f1) else "N/A",
            "Accuracy": round(float(v_acc), 5)
        })
        
    # Create LOVO Generalization DataFrame
    lovo_df = pd.DataFrame(lovo_records)
    lovo_df.to_csv(OUTPUT_LOVO_REPORT_CSV, index=False)
    print(f"\nSaved Leave-One-Virus-Out Generalization Report to: {OUTPUT_LOVO_REPORT_CSV}")
    
    # Calculate LOVO CV averages
    # Macro averages (averaging defined values across valid viruses)
    valid_aucs = [r["ROC-AUC"] for r in lovo_records if r["ROC-AUC"] != "N/A"]
    valid_pr_aucs = [r["PR-AUC"] for r in lovo_records if r["PR-AUC"] != "N/A"]
    valid_mccs = [r["MCC"] for r in lovo_records if r["MCC"] != "N/A"]
    valid_f1s = [r["F1"] for r in lovo_records if r["F1"] != "N/A"]
    valid_accs = [r["Accuracy"] for r in lovo_records]
    
    macro_auc = np.mean(valid_aucs)
    macro_pr_auc = np.mean(valid_pr_aucs)
    macro_mcc = np.mean(valid_mccs)
    macro_f1 = np.mean(valid_f1s)
    macro_acc = np.mean(valid_accs)
    
    # Micro averages (globally pooled)
    micro_auc = roc_auc_score(y, global_lovo_probs)
    micro_pr_auc = average_precision_score(y, global_lovo_probs)
    micro_mcc = matthews_corrcoef(y, global_lovo_preds)
    micro_f1 = f1_score(y, global_lovo_preds)
    micro_acc = accuracy_score(y, global_lovo_preds)
    
    print("\n" + "=" * 65)
    print("  LEAVE-ONE-VIRUS-OUT GLOBAL PERFORMANCE METRICS")
    print("=" * 65)
    print(f"  Macro-Averages (Valid species-level average):")
    print(f"    - ROC-AUC:  {macro_auc:.4f}")
    print(f"    - PR-AUC:   {macro_pr_auc:.4f}")
    print(f"    - MCC:      {macro_mcc:.4f}")
    print(f"    - F1-Score: {macro_f1:.4f}")
    print(f"    - Accuracy: {macro_acc:.4f}")
    print(f"  Micro-Averages (Globally pooled interactions):")
    print(f"    - ROC-AUC:  {micro_auc:.4f}")
    print(f"    - PR-AUC:   {micro_pr_auc:.4f}")
    print(f"    - MCC:      {micro_mcc:.4f}")
    print(f"    - F1-Score: {micro_f1:.4f}")
    print(f"    - Accuracy: {micro_acc:.4f}")
    print("=" * 65)
    
    # Train the best classifier model on the FULL dataset
    print(f"\nTraining final {best_model_name} production model on the FULL dataset...")
    if best_model_name == "LogisticRegression":
        final_clf = LogisticRegression(**best_params_template, solver='liblinear', random_state=42)
    elif best_model_name == "RandomForest":
        final_clf = RandomForestClassifier(**best_params_template, random_state=42)
    elif best_model_name == "XGBoost":
        final_clf = XGBClassifier(**best_params_template, random_state=42, eval_metric='logloss', use_label_encoder=False)
    elif best_model_name == "LightGBM":
        final_clf = LGBMClassifier(**best_params_template, random_state=42, verbose=-1)
    elif best_model_name == "CatBoost":
        final_clf = CatBoostClassifier(**best_params_template, random_state=42, verbose=0)
        
    final_clf.fit(X, y)
    
    # Save best model
    with open(OUTPUT_BEST_MODEL_PKL, "wb") as f:
        pickle.dump(final_clf, f)
    print(f"Successfully serialized final optimal model to: {OUTPUT_BEST_MODEL_PKL}")
    
    # If the best model is NOT XGBoost, we still need to serialize the best XGBoost model
    # to OUTPUT_XGB_MODEL_PKL to ensure SHAP Explainer and other tools have an XGBoost model.
    if best_model_name != "XGBoost":
        print("\nTraining companion XGBoost model for SHAP/backward-compatibility...")
        xgb_clf, xgb_params = best_cv_models["XGBoost"]
        final_xgb = XGBClassifier(**xgb_params, random_state=42, eval_metric='logloss', use_label_encoder=False)
        final_xgb.fit(X, y)
        with open(OUTPUT_XGB_MODEL_PKL, "wb") as f:
            pickle.dump(final_xgb, f)
        print(f"Successfully serialized XGBoost model to: {OUTPUT_XGB_MODEL_PKL}")
    else:
        # If XGBoost is the best model, just copy it to OUTPUT_XGB_MODEL_PKL
        with open(OUTPUT_XGB_MODEL_PKL, "wb") as f:
            pickle.dump(final_clf, f)
        print(f"Successfully cloned XGBoost model to: {OUTPUT_XGB_MODEL_PKL}")
        
    print("\n" + "=" * 65)
    print("  miRNAProtPred 2.0 TRAINING & BENCHMARK PIPELINE COMPLETE")
    print("=" * 65)

cli = main

if __name__ == "__main__":
    cli()
