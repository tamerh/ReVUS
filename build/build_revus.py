#!/usr/bin/env python
"""Build the ReVUS benchmark: every ClinVar VUS at a freeze release, labelled with its OBSERVED
germline reclassification by a strictly later release. A predict-the-transition benchmark with a
temporal holdout — contamination-resistant by construction. Deterministic; no external services.

    python build/build_revus.py \
        --old archive/variant_summary_2022-06.txt.gz \
        --cur archive/variant_summary_2026-07.txt.gz \
        --out data/revus_2022-06_to_2026-07.tsv
"""
import argparse, gzip, collections, os


def norm(s):
    s = (s or "").lower()
    if "conflict" in s: return "CONF"
    if "pathogenic" in s and "likely" in s: return "LP"
    if "pathogenic" in s: return "P"
    if "likely benign" in s: return "LB"
    if "benign" in s: return "B"
    if "uncertain" in s: return "VUS"
    return "OTHER"


def transition(cur):
    n = norm(cur)
    if n in ("P", "LP"): return "to_pathogenic"
    if n in ("B", "LB"): return "to_benign"
    if n == "CONF": return "to_conflicting"
    if n == "VUS": return "still_vus"
    return "other"


def load(path, cols, assembly="GRCh38"):
    """Load a ClinVar variant_summary.txt.gz, keyed by VariationID, one row per variant on the
    requested assembly (the file has one row per variant per assembly)."""
    with gzip.open(path, "rt") as f:
        hdr = f.readline().lstrip("#").rstrip("\n").split("\t")
        ci = {c: hdr.index(c) for c in cols}
        ia = hdr.index("Assembly")
        out = {}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(list(ci.values()) + [ia]) or p[ia] != assembly:
                continue
            out[p[ci["VariationID"]]] = {c: p[ci[c]] for c in cols}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True, help="freeze-date variant_summary.txt.gz")
    ap.add_argument("--cur", required=True, help="horizon-date variant_summary.txt.gz")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    # Use the VCF-normalised allele/position columns; ClinVar leaves ReferenceAllele /
    # AlternateAllele as "na" for nearly all rows, so those cannot be joined to a VCF source.
    keep = ["VariationID", "GeneSymbol", "Type", "ClinicalSignificance", "ReviewStatus",
            "NumberSubmitters", "Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF"]
    print("loading freeze + horizon snapshots...", flush=True)
    old = load(args.old, keep)
    cur = load(args.cur, ["VariationID", "ClinicalSignificance"])

    vus = {v: d for v, d in old.items() if norm(d["ClinicalSignificance"]) == "VUS"}
    rows = []
    for v, d in vus.items():
        rows.append({
            "variation_id": v, "gene": d["GeneSymbol"], "type": d["Type"],
            "chrom": d["Chromosome"], "pos": d["PositionVCF"],
            "ref": d["ReferenceAlleleVCF"], "alt": d["AlternateAlleleVCF"],
            "review_2022": d["ReviewStatus"], "n_submitters_2022": d["NumberSubmitters"],
            "outcome": transition(cur[v]["ClinicalSignificance"]) if v in cur else "removed",
        })

    n = len(rows)
    tr = collections.Counter(r["outcome"] for r in rows)
    print(f"\nReVUS: {n:,} VUS at freeze (GRCh38)\noutcome distribution:")
    for k, c in tr.most_common():
        print(f"   {k:16} {c:>7,}  ({100*c/n:.1f}%)")
    resolved = tr["to_pathogenic"] + tr["to_benign"]
    print(f"\n  resolved to a definite call: {resolved:,}  (path {tr['to_pathogenic']:,} :"
          f" benign {tr['to_benign']:,} = 1 : {tr['to_benign']/max(1,tr['to_pathogenic']):.1f})")

    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        cols = ["variation_id", "gene", "type", "chrom", "pos", "ref", "alt",
                "review_2022", "n_submitters_2022", "outcome"]
        with open(args.out, "w") as f:
            f.write("\t".join(cols) + "\n")
            for r in rows:
                f.write("\t".join(str(r[c]) for c in cols) + "\n")
        print(f"\n-> wrote {args.out} ({n:,} rows)")


if __name__ == "__main__":
    main()
