# Adding a dataset

The pipeline is dataset-agnostic by construction. Adding a dataset means
writing a manifest, not editing code.

**The rule that keeps it that way:** no `if dataset_id == ...` anywhere in
`src/`. Branching on *file format* is allowed (that is a generic capability);
branching on *dataset identity* is not. A dataset's peculiarities go in its
manifest, or into a recorded exclusion decision in `docs/DECISIONS.md`.

---

## 1. Write the manifest

Copy `manifests/VT2018.yaml` and edit. Required keys: `dataset_id`,
`local_path`, `file_format`, `species`, `gene_id`, `obs_columns`,
`celltype_map`, `counts_layer`.

`obs_columns.donor` and `obs_columns.celltype` are mandatory.
`compartment` must come from either `obs_columns` or `fixed_values` — the
loader refuses a manifest where it is absent from both, because a silent
missing compartment is how cross-compartment contrasts get mixed.

### Gene identity

`gene_id.var_index_type` is one of:

| value | meaning | extra keys |
|---|---|---|
| `symbol` | `var_names` are HGNC symbols | — |
| `ensembl` | `var_names` are Ensembl IDs | `symbol_column` (e.g. `feature_name`) |
| `symbol_ensembl_concat` | joined, e.g. `CCL5_ENSG00000271503` | `id_separator`, `symbol_field`, `ensembl_field` |

CELLxGENE-derived h5ad files are almost always `ensembl` with symbols in
`var['feature_name']`. Getting this wrong is not a cosmetic error: it is how
a run silently matches ~0 genes and prints a full results table anyway. The
`symbol_ensembl_concat` reader anchors the Ensembl ID to the **last** field,
so symbols containing the separator (`HLA-G`, `HLA_DRA`) survive.

Adding a fourth identity scheme means extending `VAR_INDEX_TYPES` and
`split_var_names` in `src/dnkchem/manifest.py` — still generic, still no
dataset branch.

### celltype_map

Maps published labels to this project's vocabulary. Copy labels **exactly**,
including whitespace oddities (VT2018 publishes `dNK p`, with a space).
Anything unmapped is excluded and printed by `validate_manifest.py` — read
that list, since an unmapped label is a decision, not an accident.

Reference lineages are not optional. Soup calibration needs `Myeloid`
(ruler B) and benefits from `Stromal`; the R8 noise check needs `T`. `cDC1`
is kept separate from `Myeloid` wherever the annotation supports it, because
it is the only decidual population that truly expresses `XCR1` — folding it
into `Myeloid` would corrupt the negative control.

---

## 2. Ingest

```bash
python src/ingest.py manifests/<id>.yaml
```

Reads `source.path` in `source.format` and writes the canonical h5ad at
`local_path`: cells x genes, sparse CSR, **raw integer counts**, unified obs
columns (`donor`, `compartment`, `subset`, `library`, `celltype_raw`).

Aborts if the matrix is not integer counts. Normalised or log-transformed
input is a manifest error, not something to work around.

Supported `source.format`: `h5ad`, `mtx_dir`, `dense_tsv_genes_by_cells`.

Then set `raw_counts_verified: true` in the manifest.

---

## 3. Validate — the admission gate

```bash
python src/validate_manifest.py manifests/<id>.yaml
```

Checks, in order:

1. manifest schema
2. matrix is raw integer counts
3. **panel match rate >= 90% — raises `PanelMatchError` and exits non-zero
   below that (R11).** There is no flag to downgrade this to a warning
4. every `celltype_map` target label, with its cell count
5. cells per `donor x compartment x subset`, written to
   `out/<id>/unit_counts_pre_qc.tsv`
6. UMI depth distribution, overall and for dNK cells

---

## 4. QC report

```bash
python src/qc_report.py manifests/<id>.yaml
```

Produces T5 (`qc_report.tsv`), the first part of T6 (`excluded_units.tsv`),
`annotation_check.tsv` and `nk_gate_purity.tsv`.

Exits **3** if triple-positive T contamination exceeds 15% in any dNK subset.
That is the pre-registered stop: fix the annotation question before running
the main analysis, do not push past it.

Annotation is *verified*, never repaired. If subset markers contradict the
published labels, report it and stop — re-clustering and re-labelling are
barred.

---

## 5. Only then, the main analysis

```bash
python src/detection_rates.py  manifests/<id>.yaml   # T1 (+ retention, R6)
python src/donor_tests.py      manifests/<id>.yaml   # T2
python src/soup_calibration.py manifests/<id>.yaml   # T3
python src/null_distribution.py manifests/<id>.yaml  # T4
```

`detection_rates.py` must run first: it fixes the depth floor and writes
`detection_meta.json`, which the others read so that every table refers to the
same depth matching.

---

## 6. Combining datasets

Never at the matrix level (**R9**). Each dataset runs the pipeline
independently to donor-level statistics; those statistics are combined
afterwards (random effects, or Fisher-z via `stats.fisher_z_combine`).
Integration re-imports the batch/platform/tissue confounds the design exists
to avoid.

---

## Checklist

- [ ] manifest written, `notes` records the dataset's known defects
- [ ] `ingest.py` clean, `raw_counts_verified: true`
- [ ] `validate_manifest.py` passes, unmapped labels reviewed
- [ ] `qc_report.py` passes, purity stop does not fire
- [ ] donor counts checked against `docs/PREREGISTRATION.md` §2 —
      **if n is small, `min_achievable_p` decides what may be claimed before
      any p value is looked at**
- [ ] no `if dataset_id` was added to `src/`
