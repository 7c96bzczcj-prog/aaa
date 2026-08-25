#!/usr/bin/env python3
"""Emit results/PROVENANCE.md: every delivered number back to its script,
its inputs and the exact data version it came from (spec section 6)."""
import hashlib, json, os, subprocess, sys, datetime

ROOT = "/home/user/aaa"
PRE = os.path.join(ROOT, "preflight")


def sha(p, n=1 << 20):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while True:
            b = f.read(n)
            if not b:
                break
            h.update(b)
    return h.hexdigest()[:16]


CHAIN = [
    ("manifest.yaml", "preflight/scripts/geo_scan.py + gsm_meta.py + manual "
     "verification (barcode line counts, GEO SOFT, HCA Azul API)",
     "GEO esearch/esummary/filelist.txt; HCA Azul index/projects and index/files"),
    ("results/soup_by_compartment.tsv", "preflight/scripts/task_a.py -> task_a_stats.py",
     "HCA ICA *_raw_feature_bc_matrix.h5 (dcp60); GSE233304 *_matrix.mtx.gz; "
     "GSE181543 *_raw_feature_bc_matrix.h5; cell labels from preflight/scripts/annotate.py"),
    ("results/soup_by_donor.tsv", "preflight/scripts/task_a_stats.py", "as above"),
    ("results/target_gene_detection.tsv", "preflight/scripts/task_a_stats.py",
     "admission criterion C2, per dataset x compartment"),
    ("results/early_nk_stages.tsv", "preflight/scripts/task_c.py -> task_c_stats.py "
     "-> task_c_supp.py -> assemble_task_c.py",
     "HCA ICA MantonBM1-8 (71 libraries) and GSE233304 marrow samples, unfiltered"),
    ("results/ltnk_age.tsv", "preflight/scripts/task_b.py -> task_b_stats.py",
     "GSE120221 (20 donors, ages 24-84) and HCA ICA MantonBM1-8 (ages 26-52)"),
    ("results/rho_by_lineage.tsv", "preflight/scripts/rho_by_lineage.py",
     "cross-lineage rho, to check whether the NK-specific soup excess reproduces here"),
    ("VERDICT.md", "hand-written from the tsv files; "
     "preflight/scripts/verify_verdict.py re-checks every number in it", "-"),
]

lines = ["# Provenance", "",
         f"Generated {datetime.date.today().isoformat()} by "
         "`preflight/scripts/make_provenance.py`.", "",
         "Frozen constants for this run live in `preflight/lib/pf_core.py` "
         "(`SPEC_VERSION`); the verdict lookup tables are dumped next to the "
         "results as `task_[abc]_verdict_table.json` so they can be diffed "
         "against the spec.", "",
         "| deliverable | produced by | inputs |", "|---|---|---|"]
for f, s, i in CHAIN:
    lines.append(f"| `{f}` | `{s}` | {i} |")

lines += ["", "## Data versions", ""]
try:
    hca = json.load(open(os.path.join(PRE, "scan", "hca_ica_files.json")))
    keep = [r for r in hca if r["name"].startswith(("MantonBM", "MantonBL"))
            and r["name"][8:9].isdigit()]
    lines.append(f"* HCA Census of Immune Cells, project "
                 f"`cc95ff89-2e68-4a08-a234-480eca21ce79`, catalog `dcp60`; "
                 f"{len(keep)} raw h5 files admitted, "
                 f"{sum(r['size'] for r in keep)/1e9:.2f} GB. "
                 f"Per-file sha256 in `preflight/scan/hca_ica_files.json`.")
except Exception as e:
    lines.append(f"* HCA file list unavailable ({e})")
lines.append("* GEO series downloaded from `ftp.ncbi.nlm.nih.gov/geo`; the "
             "`filelist.txt` of every scanned series is summarised in "
             "`preflight/scan/bm_scan.json` and `preflight/scan/pb_scan.json`.")

lines += ["", "## Environment", ""]
try:
    out = subprocess.run([os.path.join(ROOT, ".venv/bin/python"), "-c",
                          "import sys,numpy,scipy,pandas,scanpy,anndata,h5py,statsmodels,skimage;"
                          "print(sys.version.split()[0]);"
                          "print('numpy',numpy.__version__);print('scipy',scipy.__version__);"
                          "print('pandas',pandas.__version__);"
                          "import importlib.metadata as m;print('scanpy',m.version('scanpy'));"
                          "print('anndata',anndata.__version__);print('h5py',h5py.__version__);"
                          "print('statsmodels',statsmodels.__version__);"
                          "print('scikit-image',skimage.__version__)"],
                         capture_output=True, text=True)
    lines += ["```", "python " + out.stdout.strip().replace("\n", "\n"), "```"]
except Exception as e:
    lines.append(f"(environment capture failed: {e})")

lines += ["", "## Output checksums", "", "| file | sha256 (first 16) | bytes |", "|---|---|---|"]
for f, _, _ in CHAIN:
    p = os.path.join(ROOT, f)
    if os.path.exists(p):
        lines.append(f"| `{f}` | `{sha(p)}` | {os.path.getsize(p)} |")

open(os.path.join(ROOT, "results", "PROVENANCE.md"), "w").write("\n".join(lines) + "\n")
print("wrote results/PROVENANCE.md")
