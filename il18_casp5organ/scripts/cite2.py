"""Resolve each reference source to its PubMed record. Author/keyword queries, retried, with the
returned title printed so every hit can be eyeballed before it goes on a slide."""
import json
import subprocess
import time
import urllib.parse

CA = "/root/.ccr/ca-bundle.crt"
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def get(url, tries=4):
    for i in range(tries):
        r = subprocess.run(["curl", "-sS", "--cacert", CA, "--max-time", "60", url],
                           capture_output=True, text=True)
        if r.returncode == 0 and r.stdout.strip().startswith(("{", "<")):
            return r.stdout
        time.sleep(1.5 * (i + 1))
    return ""


def lookup(label, term, n=4):
    s = get(f"{E}/esearch.fcgi?db=pubmed&retmode=json&retmax={n}&term={urllib.parse.quote(term)}")
    try:
        ids = json.loads(s)["esearchresult"].get("idlist", [])
    except Exception:
        print(f"!! {label}: esearch failed"); return
    if not ids:
        print(f"?? {label}: NO HIT  <{term}>"); return
    s = get(f"{E}/esummary.fcgi?db=pubmed&retmode=json&id={','.join(ids)}")
    try:
        d = json.loads(s)["result"]
    except Exception:
        print(f"!! {label}: esummary failed for {ids}"); return
    print(f"== {label}")
    for p in d.get("uids", []):
        r = d[p]
        c = f"{r.get('source','')}. {r.get('pubdate','')}"
        if r.get("volume"):
            c += f";{r['volume']}"
        if r.get("issue"):
            c += f"({r['issue']})"
        if r.get("pages"):
            c += f":{r['pages']}"
        print(f"   PMID {p}  {c}.")
        print(f"      {r.get('title','')[:120]}")


Q = [
    ("CELLxGENE Discover", 'CELLxGENE[Title] AND ("Nucleic Acids Res"[TA])'),
    ("SEA-AD", 'Gabitto[Author] AND Alzheimer[Title] AND 2024[DP]'),
    ("PsychAD", 'PsychAD[All Fields]'),
    ("Tabula Muris Senis", '"Tabula Muris Senis"[All Fields]'),
    ("Tabula Sapiens (original)", '"Tabula Sapiens"[All Fields] AND Science[TA] AND 2022[DP]'),
    ("HLCA Sikkema", 'Sikkema[Author] AND lung[Title] AND 2023[DP]'),
    ("LungMAP CellRef", 'CellRef[Title]'),
    ("Siletti human brain atlas", 'Siletti[Author] AND Science[TA] AND 2023[DP]'),
    ("Yao mouse whole-brain atlas", 'Yao[Author] AND Nature[TA] AND 2023[DP] AND "mouse brain"[Title]'),
    ("Guilliams liver niches", 'Guilliams[Author] AND Cell[TA] AND 2022[DP] AND macrophage[Title]'),
    ("ImmGen open-chromatin atlas", 'Yoshida[Author] AND Cell[TA] AND 2019[DP] AND immune[Title]'),
    ("brainrnaseq mouse (Zhang 2014)", 'Zhang[Author] AND "J Neurosci"[TA] AND 2014[DP] AND splicing[Title]'),
    ("brainrnaseq human (Zhang 2016)", 'Zhang[Author] AND Neuron[TA] AND 2016[DP] AND astrocyte*[Title]'),
    ("Pressl sFANS striatum", 'Pressl[Author] AND Neuron[TA]'),
    ("Allen human cortex SMART-seq (Hodge)", 'Hodge[Author] AND Nature[TA] AND 2019[DP] AND cortex[Title]'),
    ("GTEx v8", '"GTEx Consortium"[Author] AND Science[TA] AND 2020[DP]'),
    ("Human Protein Atlas", 'Uhlen[Author] AND Science[TA] AND 2015[DP] AND proteome[Title]'),
    ("Monaco immune deconvolution", 'Monaco[Author] AND "Cell Rep"[TA] AND 2019[DP]'),
    ("Haemopedia", 'Haemopedia[All Fields]'),
    ("Snaptron", 'Snaptron[All Fields]'),
    ("recount3", 'recount3[Title]'),
    ("Ensembl 2025", '"Ensembl 2025"[Title]'),
    ("Human lung CellRef / LungMAP", 'LungMAP[All Fields] AND atlas[Title]'),
]
for lab, term in Q:
    lookup(lab, term)
    time.sleep(0.5)
