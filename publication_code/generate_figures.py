import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pickle, pandas as pd
import os
from sklearn.metrics import roc_curve, auc, precision_recall_curve

DPI = 600
FIGDIR = "/home/somenath/Pictures/Publication_somenath/mirnaprotpred2/publication_figures"

os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(os.path.join(FIGDIR, "raw_data"), exist_ok=True)

# ── Colors ────────────────────────────────────────────
C = {
    "targetscan": "#b0bec5",
    "miranda":    "#78909c",
    "rf":         "#4f98a3",
    "xgb":        "#01696f",
    "pos":        "#01696f",
    "neg":        "#e07b39",
}

# ── DATA ─────────────────────────────────────────────
tools  = ["TargetScan", "miRanda", "miRNAProtPred2\n(RF)", "miRNAProtPred2\n(XGB)"]
aucs   = [0.5176, 0.5497, 0.5325, 0.8255]
precs  = [0.5625, 0.5030, 0.7500, 1.0000]
recs   = [0.2118, 0.9765, 0.1412, 0.5059]
f1s    = [0.3077, 0.6640, 0.2376, 0.6719]
tcols  = [C["targetscan"], C["miranda"], C["rf"], C["xgb"]]

# ════════════════════════════════════════════════════
# FIGURE 1 — Grouped bar: AUC / Precision / Recall / F1
# ════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(10, 6))
x = np.arange(len(tools))
w = 0.2
metrics_data = [aucs, precs, recs, f1s]
metric_labels = ["AUC", "Precision", "Recall", "F1"]
metric_colors = ["#4f98a3", "#e07b39", "#7a5cb8", "#2d9e5e"]
for i, (vals, label, col) in enumerate(zip(metrics_data, metric_labels, metric_colors)):
    bars = ax.bar(x + i*w - 1.5*w, vals, width=w, label=label, color=col, alpha=0.9)
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.015,
                f"{v:.3f}", ha="center", va="bottom", fontsize=7.5, fontweight="bold")
ax.set_xticks(x)
ax.set_xticklabels(tools, fontsize=10)
ax.set_ylim(0, 1.25)
ax.set_ylabel("Score", fontsize=12)
ax.set_xlabel("Tool", fontsize=12)
ax.set_title("Benchmark: Viral miRNA Target Prediction Tools (n=170)", fontsize=13, fontweight="bold")
ax.legend(loc="upper left", fontsize=10, framealpha=0.9)
ax.axhline(0.5, color="red", linestyle="--", linewidth=0.8, alpha=0.5, label="Random baseline")
ax.spines[["top","right"]].set_visible(False)
plt.tight_layout()
plt.savefig(f"{FIGDIR}/Fig1_benchmark_grouped_bar.png", dpi=DPI, bbox_inches="tight")
plt.close()
print("✅ Fig1 saved")

# ════════════════════════════════════════════════════
# FIGURE 2 — ROC curves for miRNAProtPred2 XGB + RF
# ════════════════════════════════════════════════════
df_ext = pd.read_csv("output/external_test_features_final.csv")
y_true = df_ext["label"].values

fig, ax = plt.subplots(figsize=(6, 6))
model_files = [
    ("output/mirnaprotpred2_xgb.pkl",  "miRNAProtPred2 (XGB)", C["xgb"],  "-"),
    ("output/mirnaprotpred2_best.pkl", "miRNAProtPred2 (RF)",  C["rf"],   "--"),
]
for pkl, label, col, ls in model_files:
    m = pickle.load(open(pkl,"rb"))
    X = df_ext[list(m.feature_names_in_)].fillna(0)
    probs = m.predict_proba(X)[:,1]
    fpr, tpr, _ = roc_curve(y_true, probs)
    roc_auc = auc(fpr, tpr)
    ax.plot(fpr, tpr, color=col, lw=2, ls=ls, label=f"{label} (AUC={roc_auc:.3f})")
ax.plot([0,1],[0,1],"k--",lw=1,alpha=0.5,label="Random (AUC=0.500)")
ax.set_xlim([0,1]); ax.set_ylim([0,1.02])
ax.set_xlabel("False Positive Rate", fontsize=12)
ax.set_ylabel("True Positive Rate", fontsize=12)
ax.set_title("ROC Curve — miRNAProtPred 2.0 (External Test Set, n=170)", fontsize=12, fontweight="bold")
ax.legend(loc="lower right", fontsize=10)
ax.spines[["top","right"]].set_visible(False)
plt.tight_layout()
plt.savefig(f"{FIGDIR}/Fig2_ROC_curves.png", dpi=DPI, bbox_inches="tight")
plt.close()
print("✅ Fig2 saved")

# ════════════════════════════════════════════════════
# FIGURE 3 — Precision-Recall curve (XGB)
# ════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(6, 6))
for pkl, label, col, ls in model_files:
    m = pickle.load(open(pkl,"rb"))
    X = df_ext[list(m.feature_names_in_)].fillna(0)
    probs = m.predict_proba(X)[:,1]
    prec_arr, rec_arr, _ = precision_recall_curve(y_true, probs)
    pr_auc = auc(rec_arr, prec_arr)
    ax.plot(rec_arr, prec_arr, color=col, lw=2, ls=ls, label=f"{label} (AUC={pr_auc:.3f})")
baseline = y_true.mean()
ax.axhline(baseline, color="red", lw=1, ls="--", alpha=0.6, label=f"Random baseline ({baseline:.2f})")
ax.set_xlim([0,1]); ax.set_ylim([0,1.05])
ax.set_xlabel("Recall", fontsize=12)
ax.set_ylabel("Precision", fontsize=12)
ax.set_title("Precision-Recall Curve — miRNAProtPred 2.0 (n=170)", fontsize=12, fontweight="bold")
ax.legend(loc="upper right", fontsize=10)
ax.spines[["top","right"]].set_visible(False)
plt.tight_layout()
plt.savefig(f"{FIGDIR}/Fig3_precision_recall.png", dpi=DPI, bbox_inches="tight")
plt.close()
print("✅ Fig3 saved")

# ════════════════════════════════════════════════════
# FIGURE 4 — Feature importance (XGB top 15)
# ════════════════════════════════════════════════════
m = pickle.load(open("output/mirnaprotpred2_xgb.pkl","rb"))
importances = m.feature_importances_
feat_names  = list(m.feature_names_in_)
idx = np.argsort(importances)[-15:]
fig, ax = plt.subplots(figsize=(8, 6))
colors_bar = [C["xgb"] if importances[i] > np.median(importances) else C["rf"] for i in idx]
ax.barh([feat_names[i] for i in idx], importances[idx], color=colors_bar, alpha=0.9)
ax.set_xlabel("Feature Importance (Gain)", fontsize=12)
ax.set_title("Top 15 Features — miRNAProtPred 2.0 XGBoost", fontsize=13, fontweight="bold")
ax.spines[["top","right"]].set_visible(False)
plt.tight_layout()
plt.savefig(f"{FIGDIR}/Fig4_feature_importance.png", dpi=DPI, bbox_inches="tight")
plt.close()
print("✅ Fig4 saved")

# ════════════════════════════════════════════════════
# FIGURE 5 — dG distribution: positives vs negatives
# ════════════════════════════════════════════════════
train = pd.read_csv("virbase_final_dataset/virbase_cts/mirnaprotpred2_training_set.csv")
pos_dg = train[train["label"]==1]["delta_G"]
neg_dg = train[train["label"]==0]["delta_G"]
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(pos_dg, bins=30, color=C["pos"], alpha=0.7, label=f"Positive (n={len(pos_dg)})", density=True)
ax.hist(neg_dg, bins=30, color=C["neg"], alpha=0.7, label=f"Negative (n={len(neg_dg)})", density=True)
ax.axvline(pos_dg.mean(), color=C["pos"], lw=2, ls="--", label=f"Pos mean={pos_dg.mean():.1f}")
ax.axvline(neg_dg.mean(), color=C["neg"], lw=2, ls="--", label=f"Neg mean={neg_dg.mean():.1f}")
ax.set_xlabel("ΔG (kcal/mol)", fontsize=12)
ax.set_ylabel("Density", fontsize=12)
ax.set_title("ΔG Distribution: True Viral Targets vs Negatives\n(Mann-Whitney p = 2.76×10⁻⁹⁷)", fontsize=12, fontweight="bold")
ax.legend(fontsize=10)
ax.spines[["top","right"]].set_visible(False)
plt.tight_layout()
plt.savefig(f"{FIGDIR}/Fig5_dG_distribution.png", dpi=DPI, bbox_inches="tight")
plt.close()
print("✅ Fig5 saved")

# ════════════════════════════════════════════════════
# FIGURE 6 — Seed match pie: positives
# ════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(6, 5))
sizes = [19.5, 80.5]
labels = ["Seed Match\n(19.5%)", "No Seed Match\n(80.5%)"]
ax.pie(sizes, labels=labels, colors=[C["neg"], C["rf"]],
       autopct="%1.1f%%", startangle=90,
       textprops={"fontsize":12},
       wedgeprops={"edgecolor":"white","linewidth":2})
ax.set_title("Seed Match Distribution\nin Experimentally Validated Viral Targets (n=154)",
             fontsize=12, fontweight="bold")
plt.tight_layout()
plt.savefig(f"{FIGDIR}/Fig6_seed_match_pie.png", dpi=DPI, bbox_inches="tight")
plt.close()
print("✅ Fig6 saved")

print("\n🎉 All 6 publication figures saved to publication_figures/")
