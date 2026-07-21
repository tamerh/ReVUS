# ReVUS — a temporal-holdout benchmark for VUS reclassification

**ReVUS** (REclassification of Variants of Uncertain Significance) is a benchmark for a
question the field has no standard evaluation for: *given a variant classified as a Variant
of Uncertain Significance (VUS) at one point in time, will it later be reclassified — toward
pathogenic, toward benign, or remain uncertain — and can a method predict which?*

Existing longitudinal variant benchmarks score **static** pathogenic-vs-benign discrimination
on newly observed variants. None of them target the **transition itself**. ReVUS labels each
VUS by its *observed* reclassification over a temporal holdout, so it is:

- **A prediction target, not a re-encoding of static ClinVar labels** — the label is the
  realised transition of a variant that was uncertain at freeze time.
- **Contamination-resistant by construction** — a method may use only information available at
  the freeze date; the outcomes are drawn from a strictly later ClinVar release.

## The dataset

Every variant classified VUS in the **ClinVar release of 2022-06** (GRCh38), labelled by its
germline classification in the **ClinVar release of 2026-07**.

| Outcome | n | % |
|---|---:|---:|
| still VUS | 487,967 | 84.1 |
| → conflicting | 60,116 | 10.4 |
| **→ benign / likely-benign** | 18,161 | 3.1 |
| removed / retired | 8,559 | 1.5 |
| **→ pathogenic / likely-pathogenic** | 5,516 | 1.0 |
| **Total** | **580,386** | |

**Resolved to a definite call: 23,677** (pathogenic : benign = 1 : 3.3 — imbalanced; report
balanced metrics).

Schema (`data/revus_2022-06_to_2026-07.tsv`, tab-separated):

| column | description |
|---|---|
| `variation_id` | ClinVar VariationID (stable join key) |
| `gene`, `type` | gene symbol; variant type (SNV, deletion, …) |
| `chrom`, `pos`, `ref`, `alt` | GRCh38 coordinate |
| `review_2022` | ClinVar review status at freeze |
| `n_submitters_2022` | number of submitters at freeze |
| `outcome` | `still_vus` / `to_pathogenic` / `to_benign` / `to_conflicting` / `removed` / `other` |

Only freeze-time columns (`review_2022`, `n_submitters_2022`, coordinate) may be used as
inputs; `outcome` is the label.

## Tasks

1. **Resolution direction** — among variants that resolved to a definite call, predict
   pathogenic/likely-pathogenic vs benign/likely-benign. n = 23,677 (a missense subset is
   provided for sequence-based predictors). The primary clinical task.
2. **Resolution likelihood** — will a 2022 VUS resolve at all by the horizon? n = 580,386,
   ≈ 4 % positive.

## Splits and metrics

- **Split:** temporal holdout (freeze = ClinVar 2022-06; outcomes = ClinVar 2026-07). A
  *forward-living* variant — register predictions on today's VUS and auto-score them as they
  resolve — is planned as a v1 upgrade.
- **Metrics:** both tasks are imbalanced, so report **AUPRC** and **balanced accuracy**, not
  raw accuracy/AUC, and report per-stratum (variant type, freeze review status, gene).

## Baselines

A leaderboard entry must declare whether its scores were trained on ClinVar labels, because a
ClinVar-trained predictor is circularity-suspect on a ClinVar-derived target. ReVUS *exposes*
this rather than hiding it.

- **Freeze review status** is a strong non-predictor baseline for Task 2: expert-panel VUS
  resolve at ≈ 13 % vs ≈ 5 % for single-submitter VUS.
- **Sequence predictors** on Task 1 (direction): reported in `baselines/`, split by whether the
  predictor is ClinVar-trained (e.g. REVEL) or ClinVar-independent (e.g. AlphaMissense, SaProt).

See `baselines/evaluate.py` for the single-feature baselines and the current leaderboard, and
`baselines/train_baselines.py` (requires scikit-learn / pandas) for the trained baselines that
establish learnability and quantify the gene-level component of the signal via a random-versus-
gene-disjoint comparison.

### Baseline features and third-party predictor scores

`data/revus_resolved_features.tsv` ships freeze-time predictor features for the resolved
variants — **SaProt** LLR, **phyloP** conservation, and **gnomAD** allele frequency — plus
`gene`, `type`, and `review_2022` for stratification. Those three predictor columns are extracted
from a pinned **BioBTree** snapshot by `build/extract_features.py`, keyed on the GRCh38 coordinate
the label file already carries; they are *not* read from the primary sources directly (see
[Reproducing the features](#reproducing-the-features)). Only redistributable predictors are
shipped: SaProt (MIT — computed in-house from the MIT-licensed weights and ingested into BioBTree),
phyloP (public), and gnomAD allele frequency (CC0). `gnomad_af` is the global allele frequency,
falling back to the grpmax (popmax) value when the global figure is absent. Running
`baselines/evaluate.py --leaderboard` on this file reproduces the SaProt, conservation,
gnomAD-rarity, and majority-floor rows of the leaderboard directly.

The feature table covers all 23,677 resolved variants; a cell is left blank where the snapshot has
no value for that position.

**AlphaMissense** (CC BY-NC-SA) and **REVEL** (academic, non-commercial) are *not* redistributed
here. To reproduce their leaderboard rows, add `am` and `revel` columns to a local copy of the
features file — either regenerate them from BioBTree with
`python build/extract_features.py --with-restricted ...` (which emits `am` and `revel` alongside
the others; the licenses forbid *redistributing* those scores, not computing them locally), or
obtain them from dbNSFP (the paper used dbNSFP v4.x) and join by `variation_id` or
`chrom:pos:ref:alt`. `evaluate.py` then scores them automatically on the shared intersection.
Because these predictors apply only to missense variants and cover slightly different sets, the
leaderboard reports per-predictor `n` and base rate.

## Reproducing the dataset

```bash
bash build/download_clinvar.sh          # fetch the two ClinVar variant_summary archives
python build/build_revus.py \
    --old archive/variant_summary_2022-06.txt.gz \
    --cur archive/variant_summary_2026-07.txt.gz \
    --out data/revus_2022-06_to_2026-07.tsv
```

The label build is a deterministic diff of two public ClinVar releases and uses no external
services or data. (The predictor *features* are separate — a dated extract from a BioBTree
snapshot, not covered by this determinism; see below.)

## Reproducing the features

The predictor columns (`saprot`, `phylop`, `gnomad_af`) are extracted from BioBTree, keyed on the
GRCh38 coordinate the label file already carries:

```bash
# needs a running BioBTree (BIOBTREE_WS) and the `sugibiobtree` client importable
python build/extract_features.py \
    --labels data/revus_2022-06_to_2026-07.tsv \
    --out data/revus_resolved_features.tsv --workers 12
```

Unlike the labels, the features are **not** a deterministic function of the two ClinVar releases:
they are a frozen snapshot of BioBTree at extraction time, which re-ingests its sources on its own
cycle. The extractor therefore pins the snapshot in
`data/revus_resolved_features.tsv.provenance.txt` — the BioBTree version, commit and build date,
and each dataset's build date and source URL. The shipped table was extracted from BioBTree
v2.10.0 (commit `e08dd52`); reproducing it bit-for-bit needs the same snapshot, but every number in
the paper reproduces from the shipped table regardless.

## Licensing and provenance

- **Code:** MIT (`LICENSE`).
- **Labels/coordinates:** derived from ClinVar, which is public domain. Redistributed here as a
  convenience; the authoritative source is NCBI ClinVar.
- **Shipped features:** SaProt (MIT — self-computed), phyloP (public), and gnomAD allele frequency
  (CC0) only, surfaced through a pinned BioBTree snapshot (see *Reproducing the features*).
- **Not redistributed:** AlphaMissense (CC BY-NC-SA) and REVEL (non-commercial) scores; obtain
  from dbNSFP under their own licenses (see above).

## Citation

Preprint in preparation. Until then, please cite this repository.
