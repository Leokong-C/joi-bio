"""Step 4 - the external-training floor, using BindingDB instead of PDBbind General.

PDBbind General is not freely downloadable from a scripted pipeline, but the
external-training regime it represents can be reproduced from BindingDB: a
ligand-only model needs SMILES and an affinity, no structures.

Two regimes are reported, both trained on BindingDB and tested on CASF-2016:
  target-overlapping  - train on every BindingDB ligand except the CASF ones
  target-disjoint     - leave-one-target-out: for each CASF target, train only
                        on BindingDB ligands of the OTHER targets

Writes: work/external_floor.csv
"""
import argparse, json, os
import numpy as np, pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestRegressor


def build_training_table(bdb_dir, exclude_monomers):
    rows = []
    for fn in sorted(os.listdir(bdb_dir)):
        d = json.load(open(os.path.join(bdb_dir, fn), encoding="utf-8"))
        for r in list(d.values())[0].get("affinities", []):
            if r.get("affinity_type") not in ("Ki", "Kd"):
                continue
            v = str(r.get("affinity", "")).strip()
            if any(c in v for c in "<>"):
                continue
            try:
                val = float(v)
            except ValueError:
                continue
            if val <= 0:
                continue
            rows.append((fn[:-5], str(r.get("monomerid")), r.get("smile"), 9 - np.log10(val)))
    T = pd.DataFrame(rows, columns=["uniprot", "monomerid", "smiles", "pAff"]).dropna()
    T = T.groupby(["uniprot", "monomerid", "smiles"], as_index=False)["pAff"].median()
    return T[~T.monomerid.isin(exclude_monomers)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True)
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--casf-uniprot", required=True, help="work/casf_uniprot.csv from step 1")
    ap.add_argument("--matched", required=True, help="work/matched_records.csv from step 2")
    ap.add_argument("--train-features", required=True, help="featurize.py output for BindingDB")
    ap.add_argument("--test-features", required=True, help="featurize.py output for CASF ligands")
    ap.add_argument("--work", default="work")
    ap.add_argument("--trees", type=int, default=400)
    a = ap.parse_args()

    idx = pd.read_csv(a.index)[["code", "logKa"]]
    smi = pd.read_csv(a.smiles).set_index("code")["smiles"]
    excl = set(pd.read_csv(a.matched)["monomerid"].astype(str))
    T = build_training_table(os.path.join(a.work, "bindingdb"), excl)

    ftr = pd.read_csv(a.train_features).set_index("smiles")
    fte = pd.read_csv(a.test_features).set_index("smiles")
    cols = [c for c in ftr.columns if c in fte.columns]
    ftr, fte = ftr[cols], fte[cols]
    ok = ftr.notna().all() & (ftr.std() > 0)
    cols = [c for c in cols if ok.get(c, False)]

    T = T[T.smiles.isin(ftr.index)]
    Xtr_all = ftr.loc[T.smiles.values, cols].values
    ytr_all = T.pAff.values

    test = idx[idx.code.isin(smi.index)].copy()
    test["smiles"] = smi.loc[test.code].values
    test = test[test.smiles.isin(fte.index)]
    Xte = fte.loc[test.smiles.values, cols].values
    yte = test.logKa.values

    rows = []
    rf = RandomForestRegressor(n_estimators=a.trees, random_state=0, n_jobs=-1)
    rf.fit(Xtr_all, ytr_all)
    r_over = float(stats.pearsonr(yte, rf.predict(Xte))[0])
    rows.append({"regime": "target_overlapping", "n_train": len(T),
                 "n_test": len(test), "R": round(r_over, 3)})
    print(f"[04] target-overlapping R = {r_over:.3f} (train {len(T):,})", flush=True)

    cu = pd.read_csv(a.casf_uniprot)
    code2acc = {r.code: set(str(r.uniprot).split(";")) for r in cu.itertuples()}
    preds = np.full(len(test), np.nan)
    for i, code in enumerate(test.code.values):
        accs = code2acc.get(code, set())
        mask = ~T.uniprot.isin(accs).values
        if mask.sum() < 100:
            continue
        rf = RandomForestRegressor(n_estimators=a.trees, random_state=0, n_jobs=-1)
        rf.fit(Xtr_all[mask], ytr_all[mask])
        preds[i] = rf.predict(Xte[i:i + 1])[0]
    m = ~np.isnan(preds)
    r_dis = float(stats.pearsonr(yte[m], preds[m])[0])
    rows.append({"regime": "target_disjoint_LOTO", "n_train": int(len(T)),
                 "n_test": int(m.sum()), "R": round(r_dis, 3)})
    print(f"[04] target-disjoint (leave-one-target-out) R = {r_dis:.3f}", flush=True)

    pd.DataFrame(rows).to_csv(os.path.join(a.work, "external_floor.csv"), index=False)


if __name__ == "__main__":
    main()
