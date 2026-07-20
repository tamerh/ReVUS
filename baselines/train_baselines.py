#!/usr/bin/env python
"""Trained reference baselines for ReVUS: are the tasks learnable, and how much of the signal is
gene-level (ascertainment + gene-specific biology) rather than variant-level? Compares a random
split against a GENE-DISJOINT split (GroupKFold on gene); the gap between them is the amount of
apparent signal attributable to gene identity, and gene-disjoint is the benchmark's primary split.

Requires scikit-learn / pandas / numpy (not needed for the pure-Python baselines/evaluate.py).

    python baselines/train_baselines.py --data-dir data/
"""
import argparse
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, GroupKFold
from sklearn.metrics import average_precision_score, roc_auc_score

META = ["type", "review_2022", "n_submitters_2022", "chrom", "gene_vus_count"]
T1 = ["saprot", "phylop", "gnomad_af", "type", "review_2022", "gene_vus_count"]


def _design(df, cols):
    X = df[cols].copy()
    for c in cols:
        if str(X[c].dtype) in ("category", "object"):
            X[c] = X[c].astype("category").cat.codes
    return X.astype(float).values


def _oof(df, ycol, cols, split, seed=0):
    y = df[ycol].values
    oof = np.zeros(len(df))
    if split == "random":
        it = StratifiedKFold(5, shuffle=True, random_state=seed).split(df, y)
    else:
        it = GroupKFold(5).split(df, y, df.gene.values)
    X = _design(df, cols)
    for tr, te in it:
        m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.1, random_state=seed)
        m.fit(X[tr], y[tr])
        oof[te] = m.predict_proba(X[te])[:, 1]
    return y, oof


def _ap_ci(y, s, b=1000, seed=0):
    rng = np.random.default_rng(seed)
    n = len(y)
    aps = []
    for _ in range(b):
        idx = rng.integers(0, n, n)
        if y[idx].sum() == 0:
            continue
        aps.append(average_precision_score(y[idx], s[idx]))
    lo, hi = np.percentile(aps, [2.5, 97.5])
    return lo, hi


def row(tag, y, s):
    base = y.mean()
    ap = average_precision_score(y, s)
    lo, hi = _ap_ci(y, s)
    print(f"  {tag:<32} n={len(y):>7,} base={base:.4f} AUPRC={ap:.3f} "
          f"[{lo:.3f},{hi:.3f}] lift={ap/base:.2f}x ROC={roc_auc_score(y, s):.3f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data/")
    args = ap.parse_args()
    R = args.data_dir.rstrip("/") + "/"

    main = pd.read_csv(R + "revus_2022-06_to_2026-07.tsv", sep="\t", low_memory=False)
    main["y2"] = main.outcome.isin(["to_pathogenic", "to_benign"]).astype(int)
    main["y3"] = (main.outcome == "to_pathogenic").astype(int)
    main = main.join(main.groupby("gene").size().rename("gene_vus_count"), on="gene")
    for c in ("type", "review_2022", "chrom"):
        main[c] = main[c].astype("category")

    print("=== Task 2: will a freeze-time VUS resolve at all? (metadata only) ===")
    row("random split", *_oof(main, "y2", META, "random"))
    row("GENE-DISJOINT (primary)", *_oof(main, "y2", META, "group"))

    print("\n=== Task 3 (headline): resolve AND toward pathogenic? (metadata only) ===")
    row("random split", *_oof(main, "y3", META, "random"))
    row("GENE-DISJOINT (primary)", *_oof(main, "y3", META, "group"))

    feat = pd.read_csv(R + "revus_resolved_features.tsv", sep="\t", low_memory=False)
    feat = feat.join(main.groupby("gene").size().rename("gene_vus_count"), on="gene")
    feat["y1"] = (feat.outcome == "to_pathogenic").astype(int)
    sub = feat.dropna(subset=["saprot"]).copy()
    for c in ("type", "review_2022"):
        sub[c] = sub[c].astype("category")
    print(f"\n=== Task 1: direction, trained on license-clean features (n={len(sub):,}) ===")
    y1, s1 = sub.y1.values, -sub.saprot.values
    row("SaProt alone (rank)", y1, s1)
    row("trained, GENE-DISJOINT", *_oof(sub, "y1", T1, "group"))
    row("trained, random split", *_oof(sub, "y1", T1, "random"))


if __name__ == "__main__":
    main()
