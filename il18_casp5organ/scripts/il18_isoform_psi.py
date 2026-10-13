"""Is the SHORT IL-18 transcript isoform actually used, and where?

Why this exists: the main screen counts IL18 at GENE level (CELLxGENE Census stores one value per
Ensembl gene), so it cannot distinguish the full-length transcript from the shorter ones. Ensembl
annotates a shorter human isoform — IL18-211 (ENST00000945849, 148 aa) and IL18-210 (ENST00000945848,
144 aa) against the canonical IL18-201 (ENST00000280357, 193 aa). The difference is one in-frame
cassette exon.

Structure (GRCh38, chr11, MINUS strand; exons listed in transcription order):
    E1 112163906-112164094   E2 112154975-112155061   E3 112153592-112153603 (12 bp micro-exon)
    E4 112150072-112150206 (135 bp CASSETTE)          E5 112148603-112148736   E6 112143260-112143817
Skipping E4 removes exactly 45 codons -> 193 aa becomes 148 aa. Skipping E4 and the micro-exon gives
144 aa. So "short IL-18" at the transcript level is measurable as exon-E4 skipping, and every junction
that answers it shares the same acceptor at 112148737 (the E5 side). That makes PSI a clean ratio of
junctions competing for one splice site, not a comparison across different sites.

Measured here as:
    PSI = inclusion / (inclusion + exclusion)
    inclusion = reads on 112148737-112150071   (E4 -> E5, cassette present)
    exclusion = reads on 112148737-112153591   (E3 -> E5, E4 skipped      -> 148 aa)
              + reads on 112148737-112154974   (E2 -> E5, E4+micro skipped -> 144 aa)
A fourth junction sharing the same acceptor, 112148737-112163905 (E1 -> E5), skips E2/E3/E4 together
and is reported separately rather than folded into either side: it is not the annotated short isoform.

Data source: Snaptron (snaptron.cs.jhu.edu), which serves per-sample exon-exon junction counts over
recount/SRA and GTEx. Junction coordinates were not taken from memory — they were read off the Ensembl
transcript structures and then confirmed against Snaptron's own annotated junctions at this locus.

Limits, stated up front: this is TRANSCRIPT evidence, not protein. It says nothing about whether a
short IL-18 PROTEIN is made, nor about caspase-1 processing. GTEx is bulk tissue, so a tissue-level PSI
is a mixture over the cell types in that tissue and cannot be attributed to one of them.
"""
import argparse
import csv
import io
import json
import sys
import subprocess
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = "https://snaptron.cs.jhu.edu"
csv.field_size_limit(10 ** 9)

ACCEPTOR = 112148737           # shared E5-side splice site; every junction below competes for it
INCLUSION = (112148737, 112150071)
EXCLUSION = {(112148737, 112153591): "skip E4 only -> 148 aa (IL18-211)",
             (112148737, 112154974): "skip E4 + micro-exon -> 144 aa (IL18-210)"}
OTHER = {(112148737, 112163905): "E1->E5, skips E2/E3/E4 — not the annotated short isoform"}
REGION = "chr11:112148500-112164500"
MIN_READS = 20                 # a PSI from fewer than this many junction reads is noise


CA = "/root/.ccr/ca-bundle.crt"      # outbound HTTPS goes through the agent proxy; urllib does not


def get(url, out=None):              # pick up this CA bundle on its own, so shell out to curl
    cmd = ["curl", "-sS", "-m", "300", "--retry", "3", "--retry-delay", "2"]
    if Path(CA).exists():
        cmd += ["--cacert", CA]
    cmd += [url]
    r = subprocess.run(cmd, capture_output=True)
    if r.returncode != 0:
        raise RuntimeError(f"curl failed ({r.returncode}): {r.stderr.decode()[:300]}")
    if out:
        Path(out).write_bytes(r.stdout)
    return r.stdout.decode("utf-8", "replace")


def junction_rows(instance, region, cache):
    """Snaptron junction rows for one locus. Cached on disk so re-runs do not re-query."""
    p = Path(cache)
    if not p.exists():
        get(f"{BASE}/{instance}/snaptron?regions={region}", out=p)
    rows = list(csv.DictReader(io.StringIO(p.read_text()), delimiter="\t"))
    return [r for r in rows if r.get("chromosome") == "chr11"]


def per_sample(row):
    """Snaptron packs per-sample coverage as ',rail_id:count,rail_id:count,...'."""
    out = {}
    for tok in (row.get("samples") or "").strip(",").split(","):
        if ":" in tok:
            sid, cov = tok.split(":")[:2]
            out[sid] = out.get(sid, 0) + int(cov)
    return out


def sample_meta(instance, ids, cache, key_cols):
    """rail_id -> label, fetched in batches from Snaptron's samples endpoint."""
    # One ?all=1 request, cached. Batched ?ids= requests get reset by the server once the id list grows,
    # and the full table is only tens of MB, so fetching it whole is both simpler and more reliable.
    p = Path(cache)
    raw = p.with_suffix(".all.tsv")
    if not raw.exists():
        get(f"{BASE}/{instance}/samples?all=1", out=raw)
    meta = {}
    want = set(ids)
    for r in csv.DictReader(io.StringIO(raw.read_text(errors="replace")), delimiter="\t"):
        rid = r.get("rail_id")
        if rid in want:
            meta[rid] = {k: (r.get(k) or "").strip() for k in key_cols}
    p.write_text(json.dumps(meta))
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--instance", default="gtexv2", choices=["gtexv2", "srav3h"])
    args = ap.parse_args()
    work = ROOT / "data" / "isoform"
    work.mkdir(parents=True, exist_ok=True)

    rows = junction_rows(args.instance, REGION, work / f"{args.instance}_junctions.tsv")
    idx = {(int(r["start"]), int(r["end"])): r for r in rows}
    missing = [j for j in [INCLUSION, *EXCLUSION] if j not in idx]
    if missing:
        sys.exit(f"junction(s) absent from {args.instance}: {missing}")

    inc = per_sample(idx[INCLUSION])
    exc = defaultdict(int)
    for j in EXCLUSION:
        for sid, c in per_sample(idx[j]).items():
            exc[sid] += c
    oth = defaultdict(int)
    for j in OTHER:
        if j in idx:
            for sid, c in per_sample(idx[j]).items():
                oth[sid] += c

    print(f"== {args.instance} locus totals ==")
    print(f"  inclusion  {INCLUSION}  samples={idx[INCLUSION]['samples_count']:>7} "
          f"cov={idx[INCLUSION]['coverage_sum']:>10}")
    for j, what in EXCLUSION.items():
        r = idx[j]
        print(f"  exclusion  {j}  samples={r['samples_count']:>7} cov={r['coverage_sum']:>10}  {what}")
    for j, what in OTHER.items():
        if j in idx:
            r = idx[j]
            print(f"  other      {j}  samples={r['samples_count']:>7} cov={r['coverage_sum']:>10}  {what}")

    key_cols = (["SMTSD", "SMTS"] if args.instance == "gtexv2"
                else ["study", "study_title"])
    ids = set(inc) | set(exc)
    meta = sample_meta(args.instance, ids, work / f"{args.instance}_samples.json", key_cols)

    groups = defaultdict(lambda: {"n": 0, "inc": 0, "exc": 0, "oth": 0, "n_exc_pos": 0})
    per_sample_rows = []
    for sid in ids:
        i, e, o = inc.get(sid, 0), exc.get(sid, 0), oth.get(sid, 0)
        if i + e < MIN_READS:
            continue
        m = meta.get(sid, {})
        if args.instance == "gtexv2":
            label = m.get("SMTSD") or m.get("SMTS") or "(unlabelled)"
        else:   # SRA: the study accession alone says nothing — carry the title, which is the finding
            label = f"{m.get('study') or '(?)'} — {(m.get('study_title') or '')[:90]}"
        g = groups[label]
        g["n"] += 1
        g["inc"] += i
        g["exc"] += e
        g["oth"] += o
        g["n_exc_pos"] += 1 if e > 0 else 0
        per_sample_rows.append({"rail_id": sid, "label": label, "inclusion": i, "exclusion": e,
                                "other_E1_E5": o, "psi": round(i / (i + e), 5)})

    out = []
    for label, g in groups.items():
        tot = g["inc"] + g["exc"]
        out.append({"label": label, "n_samples": g["n"], "inclusion_reads": g["inc"],
                    "exclusion_reads": g["exc"], "other_E1_E5_reads": g["oth"],
                    "psi_pooled": round(g["inc"] / tot, 5) if tot else None,
                    "pct_short": round(100 * g["exc"] / tot, 4) if tot else None,
                    "n_samples_with_short": g["n_exc_pos"]})
    out.sort(key=lambda r: -(r["pct_short"] or 0))

    res = ROOT / "results"
    res.mkdir(exist_ok=True)
    for name, data, cols in [
        (f"il18_isoform_psi_{args.instance}.tsv", out,
         ["label", "n_samples", "inclusion_reads", "exclusion_reads", "other_E1_E5_reads",
          "psi_pooled", "pct_short", "n_samples_with_short"]),
        (f"il18_isoform_psi_{args.instance}_per_sample.tsv", per_sample_rows,
         ["rail_id", "label", "inclusion", "exclusion", "other_E1_E5", "psi"])]:
        with open(res / name, "w", newline="") as fh:
            fh.write(f"# IL18 cassette-exon E4 (chr11:112150072-112150206, GRCh38) PSI from Snaptron "
                     f"{args.instance}. PSI = inclusion/(inclusion+exclusion); 短 isoform = exclusion.\n"
                     f"# 仅转录本层面证据，与蛋白形式无关；bulk 组织的 PSI 是该组织内各细胞类型的混合值。\n"
                     f"# 每个样本需 inclusion+exclusion >= {MIN_READS} 条跨接读段才纳入。\n")
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(data)

    print(f"\n== top groups by short-isoform share ({args.instance}, {len(per_sample_rows)} samples pass depth) ==")
    for r in out[:20]:
        print(f"  {r['label'][:46]:46} n={r['n_samples']:>5} 短 isoform {r['pct_short']:>7.3f}%  "
              f"({r['exclusion_reads']}/{r['inclusion_reads'] + r['exclusion_reads']} 读段, "
              f"{r['n_samples_with_short']} 个样本检出)")


if __name__ == "__main__":
    main()
