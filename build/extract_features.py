#!/usr/bin/env python
"""Extract freeze-time predictor features for the ReVUS resolved variants, directly
from BioBTree — no dependency on the Sugi Variant application.

Each shipped predictor column reduces to a BioBTree `entry()` lookup keyed by the
GRCh38 coordinate already carried in the label TSV (chrom/pos/ref/alt):

    gnomad_af : entry("chrom:pos:ref:alt", "gnomad_variant")  -> af (fallback af_grpmax)
    phylop    : entry("chrom:pos",         "conservation")    -> phylop
    saprot    : entry("chrom:pos:ref:alt", "alphamissense")   -> uniprot_id, protein_variant
                entry("uniprot:protein_variant", "saprot")    -> saprot_llr

This mirrors, field-for-field, how the companion Sugi Variant resource surfaces the
same values (sugivariant/collect.py:166-171, enrich._coord_entry) — only the variant
resolution is skipped, because ReVUS already ships the coordinate. am/revel are also in
BioBTree and can be emitted with --with-restricted to regenerate the full leaderboard
locally, but are NOT written to the shipped table (their licenses forbid redistributing
the scores).

IMPORTANT: the values are a frozen snapshot of BioBTree at extraction time. BioBTree
re-ingests its sources on its own cycle, so this is transparent and re-runnable but not
bit-for-bit reproducible without the same snapshot. The BioBTree per-dataset build dates
are written to <out>.provenance.txt.

Requires a running BioBTree (BIOBTREE_WS, default http://localhost:9291) and the
`sugibiobtree` client (pip install it, or set SUGIBIOBTREE_DIR to its checkout).

    python build/extract_features.py \
        --labels data/revus_2022-06_to_2026-07.tsv \
        --out data/revus_resolved_features.tsv --workers 12
"""
import argparse, csv, json, os, sys, urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

try:
    from sugibiobtree import entry
except ImportError:  # companion client not installed — fall back to a checkout
    sys.path.insert(0, os.environ.get("SUGIBIOBTREE_DIR", "the sugi-biobtree repo"))
    from sugibiobtree import entry

RESOLVED = {"to_pathogenic", "to_benign"}


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _attrs(ident, source):
    """Attributes dict for a coordinate/protein-keyed dataset via entry() — the only
    working access for gnomad_variant / alphamissense / conservation / saprot (they are
    not map-chainable). Mirrors sugivariant/enrich._coord_entry."""
    if not ident:
        return None
    try:
        a = (entry(ident, source) or {}).get("Attributes") or {}
    except Exception:
        return None
    if not a:
        return None
    v = next(iter(a.values())) if len(a) == 1 else a   # unwrap single-key wrapper
    return v if isinstance(v, dict) else None


def features(row):
    chrom, pos, ref, alt = row.get("chrom"), row.get("pos"), row.get("ref"), row.get("alt")
    if not (chrom and pos and ref and alt):
        return {"saprot": None, "phylop": None, "gnomad_af": None}
    coord = f"{chrom}:{pos}:{ref}:{alt}"

    # gnomAD allele frequency: global af, falling back to grpmax popmax (as Sugi surfaces it).
    g = _attrs(coord, "gnomad_variant") or {}
    af, popmax = _f(g.get("af")), _f(g.get("af_grpmax"))
    gnomad_af = af if af else popmax

    # phyloP: conservation is keyed by chr:pos (ref/alt-agnostic).
    c = _attrs(f"{chrom}:{pos}", "conservation") or {}
    phylop = _f(c.get("phylop"))

    # SaProt LLR: keyed by uniprot:protein_variant, both taken from the AlphaMissense entry
    # (only when AM actually resolved a pathogenicity, mirroring enrich.alphamissense_for).
    saprot = None
    am = _attrs(coord, "alphamissense")
    if am and am.get("am_pathogenicity") is not None \
            and am.get("uniprot_id") and am.get("protein_variant"):
        s = _attrs(f"{am['uniprot_id']}:{am['protein_variant']}", "saprot") or {}
        saprot = _f(s.get("saprot_llr"))

    # AlphaMissense + REVEL scores are also in BioBTree. Extractable locally to
    # regenerate the full leaderboard, but NOT written to the shipped table:
    # their licenses (AlphaMissense CC BY-NC-SA, REVEL non-commercial) forbid
    # redistributing the scores. Only emitted under --with-restricted.
    am_score = _f(am.get("am_pathogenicity")) if am else None
    revel = _f((_attrs(coord, "revel") or {}).get("revel"))

    return {"saprot": saprot, "phylop": phylop, "gnomad_af": gnomad_af,
            "am": am_score, "revel": revel}


def _biobtree_provenance():
    """Pin the BioBTree snapshot from /ws/meta: the app version/commit/build_date
    (appparams) plus each used dataset's last_built date, status and source URL.
    Returns (base_url, appparams_dict, {group: {last_built, status, source_url}})."""
    base = os.environ.get("BIOBTREE_WS", "http://localhost:9291")
    appp, datasets = {}, {}
    try:
        with urllib.request.urlopen(base + "/ws/meta", timeout=5) as r:
            meta = json.load(r)
        appp = meta.get("appparams") or {}
        for v in (meta.get("datasets") or {}).values():
            grp = v.get("group")
            if grp and grp not in datasets:
                datasets[grp] = {"last_built": (v.get("last_built") or "")[:10],
                                 "status": v.get("status"),
                                 "source_url": v.get("source_url")}
    except Exception:
        pass
    return base, appp, datasets


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", required=True, help="ReVUS label TSV (has chrom/pos/ref/alt)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--limit", type=int, default=0, help="only the first N resolved rows (testing)")
    ap.add_argument("--drop-empty", action="store_true",
                    help="omit variants with no feature at all (saprot, phylop and gnomad_af all missing)")
    ap.add_argument("--with-restricted", action="store_true",
                    help="also emit am + revel columns (from BioBTree) for LOCAL leaderboard "
                         "regeneration only; do NOT redistribute (AlphaMissense/REVEL license)")
    args = ap.parse_args()

    rows = [r for r in csv.DictReader(open(args.labels, newline=""), delimiter="\t")
            if r["outcome"] in RESOLVED]
    if args.limit:
        rows = rows[:args.limit]
    print(f"resolved variants: {len(rows):,} — extracting from BioBTree...", flush=True)

    cols = ["variation_id", "outcome", "gene", "type", "review_2022"]
    if args.with_restricted:
        cols += ["am", "revel"]
    cols += ["saprot", "phylop", "gnomad_af"]
    out = open(args.out, "w")
    out.write("\t".join(cols) + "\n")
    done = written = 0
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(features, r): r for r in rows}
        for fut in as_completed(futs):
            r = futs[fut]
            f = fut.result()
            done += 1
            if args.drop_empty and f["saprot"] is None and f["phylop"] is None and f["gnomad_af"] is None:
                pass
            else:
                vals = {"variation_id": r["variation_id"], "outcome": r["outcome"],
                        "gene": r.get("gene", ""), "type": r.get("type", ""),
                        "review_2022": r.get("review_2022", ""), **f}
                out.write("\t".join(str(vals[c]) if vals.get(c) is not None else "" for c in cols) + "\n")
                written += 1
            if done % 2000 == 0:
                out.flush()
                print(f"  {done:,}/{len(rows):,}", flush=True)
    out.close()

    base, appp, dsets = _biobtree_provenance()
    prov = args.out + ".provenance.txt"
    with open(prov, "w") as p:
        p.write("ReVUS predictor features — provenance\n")
        p.write(f"labels: {os.path.basename(args.labels)}\n")
        p.write(f"resolved variants read: {len(rows):,}   feature rows written: {written:,}\n")
        p.write("columns saprot/phylop/gnomad_af extracted via sugibiobtree.entry() (no Sugi Variant app).\n")
        p.write("gnomad_af = global af, falling back to af_grpmax (popmax) when af is absent/zero.\n")
        p.write("am/revel omitted (dbNSFP redistribution license); obtain from dbNSFP and join.\n\n")
        if appp:
            p.write(f"BioBTree snapshot: {appp.get('biobtree_version','?')} "
                    f"(commit {appp.get('biobtree_commit','?')}, build {appp.get('biobtree_build_date','?')}) "
                    f"at {base}\n")
            p.write(f"  /ws/meta freshness_checked_at: {appp.get('freshness_checked_at','?')}\n")
        else:
            p.write(f"BioBTree snapshot: /ws/meta unreachable at {base} at run time\n")
        if dsets:
            p.write("BioBTree per-dataset snapshot (group: last_built [status] source):\n")
            for grp in ("gnomad_variant", "conservation", "alphamissense", "saprot"):
                d = dsets.get(grp)
                if d:
                    p.write(f"  {grp}: {d['last_built']} [{d['status']}] {d['source_url']}\n")
    print(f"-> wrote {args.out} ({written:,} rows); provenance -> {prov}")


if __name__ == "__main__":
    main()
