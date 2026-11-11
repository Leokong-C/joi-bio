"""Step 1 - map the CASF-2016 core set to UniProt and pull BindingDB affinities.

Writes:  work/casf_uniprot.csv, work/bindingdb/<accession>.json

Needs network access to data.rcsb.org and www.bindingdb.org.
"""
import argparse, json, os, time, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor
import pandas as pd

UA = {"User-Agent": "headroom-pipeline"}
RCSB = "https://data.rcsb.org/graphql?query="
BDB = ("https://www.bindingdb.org/rest/getLigandsByUniprots"
       "?uniprot={acc}&cutoff=1000000&response=application/json")


def get(url, timeout=600):
    return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout).read()


def pdb_to_uniprot(codes, batch=40):
    """Resolve PDB entry ids to UniProt accessions via the RCSB GraphQL endpoint."""
    out = {}
    tmpl = ("{ entries(entry_ids: %s) { rcsb_id polymer_entities "
            "{ rcsb_polymer_entity_container_identifiers "
            "{ reference_sequence_identifiers { database_accession database_name } } } } }")
    for i in range(0, len(codes), batch):
        chunk = [c.upper() for c in codes[i:i + batch]]
        d = json.loads(get(RCSB + urllib.parse.quote(tmpl % json.dumps(chunk)), timeout=120))
        for e in d["data"]["entries"]:
            accs = []
            for pe in e.get("polymer_entities") or []:
                ids = pe["rcsb_polymer_entity_container_identifiers"]
                for r in ids.get("reference_sequence_identifiers") or []:
                    if r["database_name"] == "UniProt":
                        accs.append(r["database_accession"])
            out[e["rcsb_id"].lower()] = sorted(set(accs))
        time.sleep(0.2)
    return out


def pull_bindingdb(accessions, outdir, workers=6):
    """One JSON per accession. Already-downloaded files are kept."""
    os.makedirs(outdir, exist_ok=True)

    def one(acc):
        fp = os.path.join(outdir, acc + ".json")
        if os.path.exists(fp) and os.path.getsize(fp) > 100:
            return acc, "cached"
        try:
            open(fp, "wb").write(get(BDB.format(acc=acc)))
            return acc, "ok"
        except Exception as e:                       # noqa: BLE001
            return acc, f"FAILED {type(e).__name__}"

    with ThreadPoolExecutor(workers) as ex:
        return list(ex.map(one, accessions))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True,
                    help="CSV with a 'code' column of CASF-2016 PDB ids")
    ap.add_argument("--work", default="work")
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)

    codes = pd.read_csv(a.index)["code"].astype(str).tolist()
    acc = pdb_to_uniprot(codes)
    pd.DataFrame([{"code": c, "uniprot": ";".join(v)} for c, v in acc.items()]).to_csv(
        os.path.join(a.work, "casf_uniprot.csv"), index=False)
    accessions = sorted({x for v in acc.values() for x in v})
    print(f"[01] {len(acc)} entries -> {len(accessions)} UniProt accessions", flush=True)

    res = pull_bindingdb(accessions, os.path.join(a.work, "bindingdb"))
    bad = [r for r in res if r[1].startswith("FAILED")]
    print(f"[01] BindingDB: {len(res) - len(bad)}/{len(res)} ok; failures: {bad}", flush=True)


if __name__ == "__main__":
    main()
