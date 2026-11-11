"""
Trivial multivariate null model for CASF-2016 binding affinity.

The point of this script is adversarial: it predicts logKa from 14 descriptors
that contain NO model of binding — ligand bulk/physicochemistry plus raw counts
of how many protein atoms sit near the ligand.  If a physics-based scoring
function cannot beat this, its correlation with affinity is not evidence that
it models binding.

Cross-validation is GroupKFold(10) grouped by CASF target cluster, so the four
other ligands of a target are never in the training fold with it.  Without
grouping the scores are inflated by target-level leakage.

Protein coordinates come from the prepared receptor PDBQT files; rows whose
AutoDock type is H or HD are skipped so the counts are heavy-atom counts.
"""
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from scipy import stats
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

FEATURES = ["n_heavy", "mw", "logp", "tpsa", "hbd", "hba", "n_rot", "n_arom",
            "frac_csp3", "rgyr", "c40", "c45", "c60", "c45_per_heavy"]


def receptor_heavy_coords(path):
    xyz = []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")) and line[77:79].strip() not in ("H", "HD"):
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(xyz)


def features_for(code, sdf_path, rec_path):
    mol = Chem.MolFromMolFile(sdf_path)
    lig = mol.GetConformer().GetPositions()
    heavy = np.array([a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1])
    lig = lig[heavy]
    prot = receptor_heavy_coords(rec_path)
    d = np.sqrt(((prot[:, None, :] - lig[None, :, :]) ** 2).sum(-1))
    n_heavy = mol.GetNumHeavyAtoms()
    return {
        "code": code,
        "n_heavy": n_heavy,
        "mw": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "tpsa": Descriptors.TPSA(mol),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "n_rot": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "n_arom": rdMolDescriptors.CalcNumAromaticRings(mol),
        "frac_csp3": rdMolDescriptors.CalcFractionCSP3(mol),
        "rgyr": float(np.sqrt(((lig - lig.mean(0)) ** 2).sum(1).mean())),
        "c40": int((d <= 4.0).any(1).sum()),
        "c45": int((d <= 4.5).any(1).sum()),
        "c60": int((d <= 6.0).any(1).sum()),
        "c45_per_heavy": int((d <= 4.5).any(1).sum()) / n_heavy,
    }


def within_target_rho(df, pred_col):
    out = []
    for _, g in df.groupby("target"):
        if g[pred_col].nunique() > 1 and len(g) > 2:
            out.append(stats.spearmanr(g[pred_col], g.logKa)[0])
    return float(np.nanmean(out)), len(out), float(np.mean(np.array(out) >= 0.9) * 100)


def main(feat_csv="casf2016_null_features.csv"):
    df = pd.read_csv(feat_csv)
    X, y, groups = df[FEATURES].values, df.logKa.values, df.target.values
    cv = GroupKFold(n_splits=10)
    preds = {
        "null_ridge": cross_val_predict(make_pipeline(StandardScaler(), Ridge(alpha=1.0)),
                                        X, y, cv=cv, groups=groups),
        "null_rf": cross_val_predict(
            RandomForestRegressor(n_estimators=400, random_state=42, n_jobs=-1),
            X, y, cv=cv, groups=groups),
    }
    for name, p in preds.items():
        df[name] = p
        r = stats.pearsonr(p, y)[0]
        rho_mean, n_t, pct = within_target_rho(df, name)
        print(f"{name:12s} R={r:+.3f}  rho={rho_mean:+.3f} over {n_t} targets  "
              f"rho>=0.9 in {pct:.1f}%  resid SD={np.std(y - p):.2f}")

    rf = RandomForestRegressor(n_estimators=400, random_state=42, n_jobs=-1).fit(X, y)
    imp = pd.Series(rf.feature_importances_, index=FEATURES).sort_values(ascending=False)
    imp.to_csv("casf2016_null_feature_importance.csv")
    df.to_csv("casf2016_null_models_predictions.csv", index=False)
    print("\n", imp.head(5).to_string())


if __name__ == "__main__":
    main()
