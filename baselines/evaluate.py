#!/usr/bin/env python
"""ReVUS reference baselines and metrics — pure Python, no heavy dependencies.

Task 2 (resolution likelihood) is scored here directly from the released label table using only
freeze-time columns (review status, submitter count). Task 1 (resolution direction) needs
per-variant predictor scores; pass them with --features (variation_id <tab> score) to score a
predictor, and declare whether it is ClinVar-trained.

    python baselines/evaluate.py --data data/revus_2022-06_to_2026-07.tsv
    python baselines/evaluate.py --data data/... --features revel.tsv --clinvar-trained
"""
import argparse, csv, collections

RESOLVED = {"to_pathogenic", "to_benign"}
# freeze review status ordinal (higher = stronger review)
REVIEW_RANK = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting classifications": 1,
    "no assertion criteria provided": 0,
}


def load(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def roc_auc(scores, labels):
    """Rank-based ROC-AUC (Mann-Whitney). labels in {0,1}."""
    pos = [s for s, y in zip(scores, labels) if y]
    neg = [s for s, y in zip(scores, labels) if not y]
    if not pos or not neg:
        return float("nan")
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    rsum = sum(r for r, y in zip(ranks, labels) if y)
    return (rsum - len(pos) * (len(pos) + 1) / 2.0) / (len(pos) * len(neg))


def average_precision(scores, labels):
    """Area under the precision-recall curve (AUPRC), via the step-wise AP estimator."""
    order = sorted(range(len(scores)), key=lambda i: -scores[i])
    P = sum(labels)
    if P == 0:
        return float("nan")
    # Sum over DISTINCT-score thresholds of (recall increment) x precision, collapsing tied
    # scores into a single PR point (matches sklearn.average_precision_score). Walking tied
    # items in file order would let the arbitrary ordering of a tie block bias the score — which
    # matters here because absent-from-gnomAD variants share one score at the top of the ranking.
    ap = 0.0
    tp = fp = 0
    prev_recall = 0.0
    n = len(order)
    i = 0
    while i < n:
        j = i
        while j < n and scores[order[j]] == scores[order[i]]:
            if labels[order[j]]:
                tp += 1
            else:
                fp += 1
            j += 1
        recall = tp / P
        ap += (recall - prev_recall) * (tp / (tp + fp))
        prev_recall = recall
        i = j
    return ap


def balanced_accuracy(scores, labels, thresh):
    tp = fp = tn = fn = 0
    for s, y in zip(scores, labels):
        pred = s >= thresh
        if y and pred: tp += 1
        elif y and not pred: fn += 1
        elif not y and pred: fp += 1
        else: tn += 1
    tpr = tp / (tp + fn) if (tp + fn) else 0.0
    tnr = tn / (tn + fp) if (tn + fp) else 0.0
    return (tpr + tnr) / 2.0


def task2_review_baseline(rows):
    """Resolution likelihood from freeze review status + submitter count (no external data)."""
    labels = [1 if r["outcome"] in RESOLVED else 0 for r in rows]
    print(f"\n== Task 2: resolution likelihood (n={len(rows):,}, "
          f"{sum(labels):,} resolved = {100*sum(labels)/len(rows):.1f}%) ==")
    # resolution rate by review status (the descriptive, strong non-predictor signal)
    by = collections.defaultdict(lambda: [0, 0])
    for r, y in zip(rows, labels):
        b = by[r["review_2022"]]; b[0] += y; b[1] += 1
    print("  resolution rate by freeze review status:")
    for rv, (res, tot) in sorted(by.items(), key=lambda x: -x[1][1]):
        print(f"    {rv[:44]:44} {res:>6,}/{tot:>7,} = {100*res/tot:4.1f}%")
    # as rankers
    for name, score in (("review-status rank", [REVIEW_RANK.get(r["review_2022"], 0) for r in rows]),
                        ("n_submitters", [float(r["n_submitters_2022"] or 0) for r in rows])):
        print(f"  {name:20} ROC-AUC={roc_auc(score, labels):.3f}  AUPRC={average_precision(score, labels):.3f}")


def task1_direction(rows, features, clinvar_trained, name):
    """Direction: among RESOLVED variants, P/LP (1) vs B/LB (0). Needs a predictor score."""
    idx = {r["variation_id"]: r for r in rows if r["outcome"] in RESOLVED}
    scores, labels = [], []
    for vid, sc in features.items():
        if vid in idx:
            scores.append(sc); labels.append(1 if idx[vid]["outcome"] == "to_pathogenic" else 0)
    if not scores:
        print("\n== Task 1: no overlap between features and resolved variants =="); return
    tag = "ClinVar-TRAINED (circularity-suspect)" if clinvar_trained else "ClinVar-independent"
    print(f"\n== Task 1: resolution direction — {name} [{tag}] ==")
    print(f"  n={len(scores):,} resolved ({sum(labels):,} path : {len(labels)-sum(labels):,} benign)")
    print(f"  ROC-AUC={roc_auc(scores, labels):.3f}  AUPRC={average_precision(scores, labels):.3f}")


# Task-1 predictors: (column, display, orient-to-"higher=path", clinvar_trained, threshold).
# SaProt LLR is more negative = more damaging, so it is negated. gnomAD AF: rarer = more likely
# pathogenic, so it is negated. Thresholds are conventions, not universal standards (see paper
# footnotes): AM 0.564, REVEL 0.5 (ClinGen SVI uses 0.644/0.773), SaProt -10.0 is an in-house
# ClinVar-calibrated cut, phyloP 2.0, gnomAD AF 1e-4.
SEQ = ["am", "revel", "saprot"]                 # the paired-intersection sequence predictors
PREDICTORS = [
    ("am",     "AlphaMissense",        lambda x: x,  False, 0.564),
    ("revel",  "REVEL",                lambda x: x,  True,  0.5),
    ("saprot", "SaProt",               lambda x: -x, False, 10.0),
    ("phylop", "Conservation (phyloP)", lambda x: x, False, 2.0),
    ("gnomad_af", "gnomAD rarity",     lambda x: -x, False, -1e-4),
]


def _boot_ci_ap(scores, labels, b=1000, seed=0):
    import random
    rng = random.Random(seed)
    n = len(scores)
    aps = []
    for _ in range(b):
        idx = [rng.randrange(n) for _ in range(n)]
        aps.append(average_precision([scores[i] for i in idx], [labels[i] for i in idx]))
    aps = sorted(a for a in aps if a == a)
    if not aps:
        return (float("nan"), float("nan"))
    return (aps[int(0.025 * len(aps))], aps[int(0.975 * len(aps))])


def _score_one(rows, col, orient, thr):
    # For gnomAD rarity, absence from gnomAD is the STRONGEST rarity signal, so a missing AF is
    # scored as absent (af=0, rarest) rather than dropped; other predictors skip missing values.
    fill_absent = (col == "gnomad_af")
    scores, labels = [], []
    for r in rows:
        v = r.get(col, "")
        if v not in ("", None):
            scores.append(orient(float(v)))
            labels.append(1 if r["outcome"] == "to_pathogenic" else 0)
        elif fill_absent:
            scores.append(orient(0.0))
            labels.append(1 if r["outcome"] == "to_pathogenic" else 0)
    if not scores:
        return None
    base = sum(labels) / len(labels)
    ap = average_precision(scores, labels)
    lo, hi = _boot_ci_ap(scores, labels)
    # Sanity check: a correctly-computed point estimate must lie inside its bootstrap CI. When it
    # does not (as with the earlier tie-ordering bug on absence-imputed scores), the estimator is
    # order-dependent and the AUPRC is unreliable.
    if not (lo - 1e-9 <= ap <= hi + 1e-9):
        print(f"WARNING [{col}]: point estimate {ap:.4f} lies outside its bootstrap CI "
              f"[{lo:.4f},{hi:.4f}] -- likely a tie/ordering artifact in average_precision().")
    return dict(n=len(scores), base=base, ap=ap, lo=lo, hi=hi,
                auc=roc_auc(scores, labels), bacc=balanced_accuracy(scores, labels, thr))


def _row(name, tag, m):
    # AUPRC has no fixed chance level; report base rate and lift (AUPRC/base) so figures are
    # comparable across predictors evaluated on different-coverage sets.
    print(f"{name:24} {tag:8} {m['n']:>6,} {m['base']:>5.3f} {m['ap']:>6.3f} "
          f"[{m['lo']:.3f},{m['hi']:.3f}] {m['ap']/m['base']:>5.2f} {m['auc']:>7.3f} {m['bacc']:>7.3f}")


def leaderboard(path):
    """Task-1 direction leaderboard. Sequence predictors are scored on the SHARED intersection
    (all of am/revel/saprot present) so they are compared on identical data; predictors with
    broader coverage (phyloP, gnomAD rarity) are additionally scored on their own coverage, with
    the per-set base rate and lift reported. Path=1, benign=0."""
    rows = load(path)
    # sequence predictors actually present in this file (AlphaMissense/REVEL are not
    # redistributed; a user who adds them from dbNSFP recovers the full paired comparison).
    seq_present = [c for c in SEQ if any(r.get(c) not in ("", None) for r in rows)]
    paired = [r for r in rows if all(r.get(c) not in ("", None) for c in seq_present)]
    hdr = f"{'Predictor':24} {'trained':8} {'n':>6} {'base':>5} {'AUPRC':>6} {'95% CI':>15} {'lift':>5} {'ROC-AUC':>7} {'bal.acc':>7}"

    print(f"\n== Task 1: matched comparison on the sequence-predictor intersection "
          f"(n={len(paired):,}) ==")
    print(hdr); print("-" * len(hdr))
    for col, name, orient, trained, thr in PREDICTORS:
        m = _score_one(paired, col, orient, thr)   # every predictor on the SAME variant set
        if m:
            _row(name, "yes" if trained else "no", m)
    # majority-class floor on the same set
    base = sum(1 for r in paired if r["outcome"] == "to_pathogenic") / max(1, len(paired))
    print(f"{'majority (always benign)':24} {'--':8} {len(paired):>6,} {base:>5.3f} {base:>6.3f} "
          f"{'--':>15} {'1.00':>5} {'0.500':>7} {'0.500':>7}")

    print(f"\n== broader-coverage predictors on their own coverage ==")
    print(hdr); print("-" * len(hdr))
    for col, name, orient, trained, thr in PREDICTORS:
        if col in SEQ:
            continue
        m = _score_one(rows, col, orient, thr)
        if m:
            _row(name, "yes" if trained else "no", m)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--features", default="", help="TSV: variation_id <tab> score (single predictor)")
    ap.add_argument("--leaderboard", default="", help="multi-column features TSV -> Task-1 leaderboard")
    ap.add_argument("--clinvar-trained", action="store_true")
    ap.add_argument("--name", default="predictor")
    args = ap.parse_args()

    rows = load(args.data)
    tr = collections.Counter(r["outcome"] for r in rows)
    print(f"ReVUS: {len(rows):,} variants; outcomes {dict(tr)}")
    task2_review_baseline(rows)

    if args.features:
        feats = {}
        with open(args.features) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    try: feats[p[0]] = float(p[1])
                    except ValueError: pass
        task1_direction(rows, feats, args.clinvar_trained, args.name)
    elif args.leaderboard:
        leaderboard(args.leaderboard)
    else:
        print("\n(Task 1: pass --leaderboard <features.tsv> for the predictor leaderboard.)")


if __name__ == "__main__":
    main()
