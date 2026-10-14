import json, subprocess, time, urllib.parse
CA="/root/.ccr/ca-bundle.crt"; E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
def get(u,t=6):
    for i in range(t):
        r=subprocess.run(["curl","-sS","--cacert",CA,"--max-time","60",u],capture_output=True,text=True)
        if r.returncode==0 and r.stdout.strip().startswith("{"): return r.stdout
        time.sleep(2.5*(i+1))
    return ""
def find(lab,term,n=5):
    s=get(f"{E}/esearch.fcgi?db=pubmed&retmode=json&retmax={n}&term={urllib.parse.quote(term)}")
    try: ids=json.loads(s)["esearchresult"].get("idlist",[])
    except Exception: print(f"!! {lab}: esearch failed"); return
    if not ids: print(f"?? {lab}: NO HIT <{term}>"); return
    time.sleep(1.2)
    s=get(f"{E}/esummary.fcgi?db=pubmed&retmode=json&id={','.join(ids)}")
    try: d=json.loads(s)["result"]
    except Exception: print(f"!! {lab}: esummary failed"); return
    print(f"== {lab}")
    for p in d.get("uids",[]):
        r=d[p]; c=f"{r.get('source','')}. {r.get('pubdate','')}"
        if r.get("volume"): c+=f";{r['volume']}"
        if r.get("issue"): c+=f"({r['issue']})"
        if r.get("pages"): c+=f":{r['pages']}"
        print(f"   PMID {p}  {c}.\n      {r.get('title','')[:120]}")
    time.sleep(1.2)
find("Tabula Muris Senis", '"ageing tissues"[Title] AND Nature[TA]')
find("CellRef lung reference", '"single cell reference"[Title] AND lung[Title]')
find("PsychAD prefrontal", 'PsychAD[All Fields] AND ("prefrontal"[Title] OR "lifespan"[Title])')
find("SEA-AD confirm", '39402379[UID]')
