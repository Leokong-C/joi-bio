"""Step 2 - label-noise ceiling for CASF-2016 from independent BindingDB measurements.

Match CASF ligands to BindingDB molecules, convert nM to pAffinity, and estimate
how much of the label variance is measurement error.

Estimator (Kalliokoski et al. 2013, doi:10.1371/journal.pone.0061007):
    sigma_error = SD(differences) / sqrt(2)
    ceiling_R   = sqrt(1 - sigma_error^2 / SD(label)^2)

Writes: work/matched_records.csv, work/ceiling.json
"""
import argparse, json, os
import numpy as np, pandas as pd
from scipy import stats
from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
RDLogger.DisableLog("rdApp.*")


def load_bindingdb(outdir):
    rows = []
    for fn in sorted(os.listdir(outdir)):
        d = json.load(open(os.path.join(outdir, fn), encoding="utf-8"))
        for r in list(d.values())[0].get("affinities", []):   # key is misspelled upstream
            rows.append((fn[:-5], str(r.get("monomerid")), r.get("smile"),
                         r.get("affinity_type"), str(r.get("affinity")),
                         r.get("doi"), r.get("pmid")))
    return pd.DataFrame(rows, columns=["uniprot", "monomerid", "smiles",
                                       "affinity_type", "affinity", "doi", "pmid"])


def match(casf_smiles, bdb):
    """Formula pre-filter, then InChIKey. Returns {monomerid: (code, level)}."""
    ours, formulas = {}, set()
    for code, smi in casf_smiles.items():
        m = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
        if m is None:
            continue
        ik = Chem.MolToInchiKey(m)
        ours[ik] = code
        formulas.add(CalcMolFormula(m))
    skel = {}
    for ik, code in ours.items():
        skel.setdefault(ik.split("-")[0], []).append(code)

    uniq = bdb.drop_duplicates("monomerid")[["monomerid", "smiles"]]
    out = {}
    for mid, smi in zip(uniq.monomerid, uniq.smiles):
        m = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
        if m is None or CalcMolFormula(m) not in formulas:
            continue
        ik = Chem.MolToInchiKey(m)
        if ik in ours:
            out[mid] = (ours[ik], "exact")
        else:
            s = ik.split("-")[0]
            if s in skel:
                out[mid] = (skel[s][0], "skeleton" if len(skel[s]) == 1 else "skeleton_ambiguous")
    return out


def ceiling(label_sd, diffs):
    sd_diff = float(np.std(diffs, ddof=1))
    sigma_e = sd_diff / np.sqrt(2)
    return sd_diff, sigma_e, float(np.sqrt(max(0.0, 1 - sigma_e ** 2 / label_sd ** 2)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True, help="CSV with code, logKa")
    ap.add_argument("--smiles", required=True, help="CSV with code, smiles")
    ap.add_argument("--work", default="work")
    a = ap.parse_args()

    idx = pd.read_csv(a.index)[["code", "logKa"]]
    smi = pd.read_csv(a.smiles).set_index("code")["smiles"].to_dict()
    bdb = load_bindingdb(os.path.join(a.work, "bindingdb"))
    hits = match(smi, bdb)
    print(f"[02] matched {len(hits)} BindingDB molecules", flush=True)

    m = bdb[bdb.monomerid.isin(hits)].copy()
    m["code"] = [hits[x][0] for x in m.monomerid]
    m["match_level"] = [hits[x][1] for x in m.monomerid]
    m.to_csv(os.path.join(a.work, "matched_records.csv"), index=False)

    raw = m.affinity.astype(str)
    m["val_nM"] = pd.to_numeric(raw.str.replace(r"[<>~\s]", "", regex=True), errors="coerce")
    m = m[(~raw.str.contains(r"[<>]")) & m.val_nM.notna() & (m.val_nM > 0)
          & (m.match_level != "skeleton_ambiguous")]
    m = m[m.affinity_type.isin(["Ki", "Kd"])].copy()
    m["pAff"] = 9 - np.log10(m.val_nM)

    g = m.groupby("code")["pAff"].median().rename("bdb").reset_index().merge(idx, on="code")
    sd_label = float(g.logKa.std(ddof=1))
    sd_diff, sigma_e, cross = ceiling(sd_label, g.bdb - g.logKa)

    m["ref"] = m.doi.fillna("").astype(str) + "|" + m.pmid.fillna("").astype(str)
    rep = m.groupby(["code", "ref"], as_index=False)["pAff"].median()
    rep = rep.groupby("code")["pAff"].agg(["count", "std"])
    rep = rep[(rep["count"] >= 3) & rep["std"].notna()]
    pooled = float(np.sqrt(np.average(rep["std"] ** 2, weights=rep["count"] - 1)))
    replicate = float(np.sqrt(max(0.0, 1 - pooled ** 2 / sd_label ** 2)))

    res = {"n_complexes": int(len(g)), "n_measurements": int(len(m)),
           "pearson_R_sources": float(stats.pearsonr(g.logKa, g.bdb)[0]),
           "bias": float((g.bdb - g.logKa).mean()), "sd_diff": sd_diff,
           "sigma_error_cross_source": sigma_e, "ceiling_cross_source": cross,
           "pooled_within_sd": pooled, "ceiling_replicate_only": replicate,
           "shared_measurement_fraction": float(((g.bdb - g.logKa).abs() < 0.02).mean()),
           "headline_ceiling": [round(replicate, 3), 0.90],
           "note": ("cross-source ceiling is optimistic when BindingDB re-reports the same "
                    "primary measurement the label was curated from; the replicate-only "
                    "estimate never touches the label. Kramer et al. 2012 (10.1021/jm300131x) "
                    "put the independent maximum at R^2 = 0.81, i.e. R = 0.90.")}
    json.dump(res, open(os.path.join(a.work, "ceiling.json"), "w"), indent=1)
    print(json.dumps(res, indent=1), flush=True)


if __name__ == "__main__":
    main()
