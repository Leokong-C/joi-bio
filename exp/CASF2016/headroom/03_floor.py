"""Step 3 - the floor: how well do models that never see the interaction do?

Three bias-only baselines (after Durant et al., Bioinformatics 41(2) btaf040):
  LigandBias   - RDKit 2D descriptors of the ligand only
  ProteinBias  - amino-acid composition of the binding-site pocket only
  BothBias     - the two concatenated

Each is scored under two cross-validation splits on the SAME data:
  target-disjoint (GroupKFold by target)  - the target is never in training
  random (KFold, shuffled)                - other ligands of the target are in training

The difference between the two is the size of target-identity memorisation.

Writes: work/floor.csv, work/split_effect.csv
"""
import argparse, collections, os
import numpy as np, pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold, KFold, cross_val_predict

AA = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
      "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]


def pocket_composition(pdb_dir):
    """Residue counts per entry from trimmed binding-site PDB files."""
    rows = []
    for fn in sorted(os.listdir(pdb_dir)):
        if not fn.endswith(".pdb"):
            continue
        seen, cnt = set(), collections.Counter()
        for line in open(os.path.join(pdb_dir, fn), errors="replace"):
            if line.startswith("ATOM"):
                key = (line[21], line[22:27])
                if key not in seen:
                    seen.add(key)
                    cnt[line[17:20].strip()] += 1
        row = {a: cnt.get(a, 0) for a in AA}
        row["n_res"] = sum(row.values())
        row["code"] = fn[:-4]
        rows.append(row)
    return pd.DataFrame(rows).set_index("code")


def clean(df):
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.loc[:, df.notna().all()]
    return df.loc[:, df.std() > 0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True, help="CSV with code, logKa, target")
    ap.add_argument("--ligand-features", required=True,
                    help="CSV from featurize.py keyed by the ligand SMILES")
    ap.add_argument("--smiles", required=True, help="CSV with code, smiles")
    ap.add_argument("--pockets", required=True, help="directory of trimmed pocket PDBs")
    ap.add_argument("--work", default="work")
    ap.add_argument("--trees", type=int, default=500)
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)

    idx = pd.read_csv(a.index)[["code", "logKa", "target"]].set_index("code")
    smi = pd.read_csv(a.smiles).set_index("code")["smiles"]
    feats = pd.read_csv(a.ligand_features).set_index("smiles")
    L = clean(feats.reindex(smi.loc[idx.index].values).set_index(idx.index))
    K = clean(pocket_composition(a.pockets).loc[idx.index])
    B = pd.concat([L, K], axis=1)

    y = idx.logKa.values
    g = idx.target.values
    rows = []
    for name, X in [("LigandBias", L), ("ProteinBias", K), ("BothBias", B)]:
        r = {}
        for split, cv, kw in [("target_disjoint", GroupKFold(10), {"groups": g}),
                              ("random", KFold(10, shuffle=True, random_state=0), {})]:
            rf = RandomForestRegressor(n_estimators=a.trees, random_state=0, n_jobs=-1)
            p = cross_val_predict(rf, X.values, y, cv=cv, **kw)
            r[split] = round(float(stats.pearsonr(y, p)[0]), 3)
        rows.append({"model": name, "n_features": X.shape[1],
                     "R_target_disjoint": r["target_disjoint"],
                     "R_random_split": r["random"],
                     "gain_from_seeing_target": round(r["random"] - r["target_disjoint"], 3)})
        print(f"[03] {name:12} disjoint {r['target_disjoint']:.3f}  random {r['random']:.3f}", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(a.work, "floor.csv"), index=False)
    print(out.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
