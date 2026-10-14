import json
import subprocess
import time
import urllib.parse

CA = "/root/.ccr/ca-bundle.crt"
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def get(url, tries=6):
    for i in range(tries):
        r = subprocess.run(["curl", "-sS", "--cacert", CA, "--max-time", "60", url],
                           capture_output=True, text=True)
        if r.returncode == 0 and r.stdout.strip().startswith("{"):
            return r.stdout
        time.sleep(2.0 * (i + 1))
    return ""


def show(label, pmids):
    s = get(f"{E}/esummary.fcgi?db=pubmed&retmode=json&id={','.join(pmids)}")
    try:
        d = json.loads(s)["result"]
    except Exception:
        print(f"!! {label}: esummary failed"); return
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
        print(f"   PMID {p}  {c}.\n      {r.get('title','')[:120]}")
    time.sleep(1.0)


def find(label, term, n=5):
    s = get(f"{E}/esearch.fcgi?db=pubmed&retmode=json&retmax={n}&term={urllib.parse.quote(term)}")
    try:
        ids = json.loads(s)["esearchresult"].get("idlist", [])
    except Exception:
        print(f"!! {label}: esearch failed <{term}>"); return
    if not ids:
        print(f"?? {label}: NO HIT <{term}>"); return
    time.sleep(1.0)
    show(label, ids)


show("Tabula Sapiens (original)", ["35549404"])
show("Siletti candidates", ["37824674", "37824663", "37824655", "37824643"])
show("Pressl", ["38237588"])
show("Snaptron candidates", ["28968689", "35188075"])
find("Tabula Muris Senis (Nature 2020)", '"Tabula Muris"[Title] AND Nature[TA] AND 2020[DP]')
find("HLCA Sikkema", 'Sikkema[Author] AND "Nat Med"[TA]')
find("Yao mouse whole brain", 'Yao[Author] AND Nature[TA] AND atlas[Title] AND "2023/01/01"[PDAT] : "2024/06/30"[PDAT]')
find("recount3", 'recount3[All Fields] AND "Genome Biol"[TA]')
find("GTEx v8", 'GTEx[Title] AND Science[TA]')
find("Ensembl latest", 'Ensembl[Title] AND "Nucleic Acids Res"[TA] AND "2024/09/01"[PDAT] : "2025/12/31"[PDAT]')
find("CellRef lung", 'Guo[Author] AND "Nat Commun"[TA] AND lung[Title] AND reference[Title]')
find("PsychAD prefrontal atlas", 'PsychAD[All Fields] AND prefrontal[Title]')
