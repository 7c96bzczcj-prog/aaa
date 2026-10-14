"""Resolve each reference database/dataset to its PubMed record and print a formatted citation.

Nothing here is written from memory: every field (journal, year, volume, issue, pages) comes back from
NCBI E-utilities. A query that returns no hit, or more than one, is printed as such and left for review
rather than guessed at.
"""
import json
import subprocess
import sys
import time
import urllib.parse

CA = "/root/.ccr/ca-bundle.crt"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def get(url):
    r = subprocess.run(["curl", "-sS", "--cacert", CA, "--max-time", "60", url],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[:200])
    return r.stdout


def esearch(term, retmax=3):
    u = f"{EUTILS}/esearch.fcgi?db=pubmed&retmode=json&retmax={retmax}&term={urllib.parse.quote(term)}"
    return json.loads(get(u))["esearchresult"].get("idlist", [])


def esummary(pmids):
    u = f"{EUTILS}/esummary.fcgi?db=pubmed&retmode=json&id={','.join(pmids)}"
    d = json.loads(get(u))["result"]
    return [d[p] for p in d.get("uids", [])]


def fmt(rec):
    """Journal. Year Mon Day;Vol(Issue):Pages — exactly the fields PubMed returns."""
    j = rec.get("source", "")
    dp = rec.get("pubdate", "")
    vol = rec.get("volume", "")
    iss = rec.get("issue", "")
    pg = rec.get("pages", "")
    s = f"{j}. {dp}"
    if vol:
        s += f";{vol}"
    if iss:
        s += f"({iss})"
    if pg:
        s += f":{pg}"
    return s + "."


# label -> (query, expected-journal-substring or None). DOIs come from datasets.yaml / the project's
# own verified records where available; otherwise an exact-title search.
QUERIES = [
    ("CZ CELLxGENE Discover", '"CZ CELLxGENE Discover"[Title] AND 2025[DP]', "Nucleic Acids Res"),
    ("Allen Human Variation Study (Jorstad 2023)", "10.1126/science.adf2359[AID]", "Science"),
    ("SEA-AD", '"Seattle Alzheimer\'s Disease Brain Cell Atlas"[Title]', None),
    ("PsychAD", '"PsychAD"[Title]', None),
    ("Tabula Muris Senis", '"Tabula Muris Senis"[Title]', "Nature"),
    ("Tabula Sapiens", '"Tabula Sapiens"[Title]', "Science"),
    ("Human Lung Cell Atlas (Sikkema 2023)", '"integrated cell atlas of the human lung"[Title]', "Nat Med"),
    ("LungMAP CellRef", '"LungMAP"[Title] AND "CellRef"[Title]', None),
    ("Sathyamurthy 2018 spinal cord", "29466745[UID]", "Cell Rep"),
    ("Siletti 2023 human brain atlas", '"Transcriptomic diversity of cell types across the adult human brain"[Title]', "Science"),
    ("Yao 2023 mouse CNS taxonomy", '"high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain"[Title]', "Nature"),
    ("Guilliams 2022 liver", '"Spatial proteogenomics reveals distinct and evolutionarily conserved hepatic macrophage niches"[Title]', "Cell"),
    ("Human pancreas islet multiomics", "10.1038/s42255-022-00531-x[AID]", None),
    ("Mouse islet cross-species map", "10.1016/j.molmet.2022.101595[AID]", None),
    ("CNS white matter (human spinal cord WM)", "10.1186/s40478-023-01568-z[AID]", None),
    ("ImmGen GSE109125 (Yoshida 2019)", '"The cis-Regulatory Atlas of the Mouse Immune System"[Title]', "Cell"),
    ("brainrnaseq mouse (Zhang 2014)", '"An RNA-sequencing transcriptome and splicing database of glia, neurons, and vascular cells of the cerebral cortex"[Title]', "J Neurosci"),
    ("brainrnaseq human (Zhang 2016)", '"Purification and Characterization of Progenitor and Mature Human Astrocytes"[Title]', "Neuron"),
    ("Pressl sFANS", '"Pressl C"[Author] AND striat*[Title]', "Neuron"),
    ("Allen human cortex SMART-seq (Hodge 2019)", '"Conserved cell types with divergent features in human versus mouse cortex"[Title]', "Nature"),
    ("GTEx v8", '"The GTEx Consortium atlas of genetic regulatory effects across human tissues"[Title]', "Science"),
    ("Human Protein Atlas (Uhlen 2015)", '"Tissue-based map of the human proteome"[Title]', "Science"),
    ("Monaco 2019 immune cells", '"RNA-Seq Signatures Normalized by mRNA Abundance Allow Absolute Deconvolution of Human Immune Cell Types"[Title]', "Cell Rep"),
    ("Haemopedia", '"Haemopedia RNA-seq"[Title]', "Nucleic Acids Res"),
    ("Snaptron", '"Snaptron: querying splicing patterns across tens of thousands of RNA-seq samples"[Title]', "Bioinformatics"),
    ("recount3", '"recount3: summaries and queries for large-scale RNA-seq expression and splicing"[Title]', "Genome Biol"),
    ("Ensembl 2025", '"Ensembl 2025"[Title]', "Nucleic Acids Res"),
]

for label, term, expect in QUERIES:
    try:
        ids = esearch(term)
    except Exception as e:
        print(f"!! {label}: query failed — {e}")
        continue
    if not ids:
        print(f"?? {label}: NO HIT  ({term})")
        time.sleep(0.4)
        continue
    recs = esummary(ids)
    for r in recs:
        flag = ""
        if expect and expect.lower() not in r.get("source", "").lower():
            flag = "  <-- journal differs from expected"
        print(f"{label}\n    PMID {r['uid']}  {fmt(r)}{flag}\n    {r.get('title','')[:110]}")
    if len(recs) > 1:
        print(f"    ^^ {len(recs)} hits for {label} — pick one")
    time.sleep(0.4)
