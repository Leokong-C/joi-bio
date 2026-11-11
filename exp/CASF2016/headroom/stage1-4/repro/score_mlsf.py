"""
RF-Score (ODDT) on the CASF-2016 core set — Windows-compatible.

Three non-obvious workarounds are encoded here; without them the script fails:

1.  ODDT must be installed with ``pip install oddt --no-deps``.  Its doc
    dependency ``sphinxcontrib-jsmath`` ships a ``.pth`` file that executes on
    every interpreter start, which the Claude Science sandbox refuses.  ODDT's
    real runtime dependencies (numpy, scipy, scikit-learn, pandas, rdkit) are
    already present in the ``docking`` environment.
2.  ``np.in1d`` was removed in NumPy 2.0; ODDT 0.4 still calls it.  The shim
    below must run BEFORE ``import oddt``.
3.  ODDT's Open Babel backend imports the legacy top-level ``pybel`` module
    (Open Babel 2.x API).  Open Babel 3.x exposes ``openbabel.pybel`` instead,
    so the RDKit backend is used — it reads PDB and SDF but NOT PDBQT.

Performance: ODDT rebuilds the protein ``atom_dict`` by iterating residues in
Python, which costs ~60 s for a whole protein.  Trimming the receptor to
residues within 16 A of any ligand atom cuts this to ~1.5 s and leaves every
prediction bit-identical (verified: max |diff| = 1.8e-15), because the largest
descriptor cutoff is 12 A.

RF-Score v3 is deliberately excluded: its AutoDock Vina term uses
``np.hstack`` on a generator, which NumPy 2 rejects.  v2 scores highest of the
three versions anyway.
"""
import json, os, time
import numpy as np
np.in1d = np.isin                      # workaround 2 — must precede import oddt
import pandas as pd
from rdkit import Chem
import oddt
from oddt.scoring.functions import RFScore

MARGIN = 16.0                          # A around the ligand, >= 12 A cutoff + slack
LIG_OVERRIDE = {                       # see REPRODUCE.md "Two corrected complexes"
    "3syr": ("sdf", "casf/lig_fix/3syr.sdf"),
    "1a30": ("pdb", "casf/lig_fix/1a30_pep.pdb"),
}


def ligand_coords(fmt, path):
    mol = (Chem.MolFromMolFile(path, removeHs=False, sanitize=False) if fmt == "sdf"
           else Chem.MolFromPDBFile(path, removeHs=False, sanitize=False))
    return mol.GetConformer().GetPositions()


def trim_receptor(code, lig_xyz, margin=MARGIN):
    """Keep whole residues with any atom within `margin` of any ligand atom."""
    out = f"casf/prot_trim/{code}.pdb"
    if os.path.exists(out):
        return out
    lines = [l for l in open(f"casf/prot/{code}.pdb") if l.startswith("ATOM  ")]
    xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lines])
    res = np.array([l[21:27] for l in lines])          # chain + resSeq + iCode
    dmin = np.sqrt(((xyz[:, None, :] - lig_xyz[None, :, :]) ** 2).sum(-1)).min(1)
    keep = set(res[dmin <= margin])
    os.makedirs("casf/prot_trim", exist_ok=True)
    open(out, "w").write("".join(l for l, k in zip(lines, res) if k in keep) + "END\n")
    return out


def main(log="casf2016_mlsf_log.jsonl"):
    # train_set='refined' & ~core  ->  the CASF-2016 core set is NEVER in training.
    # First call downloads PDBbind v2016 descriptors and trains (~25 s per version),
    # then caches RFScore_v{1,2}_pdbbind2016.pickle in the working directory.
    models = {v: RFScore.rfscore.load(version=v) for v in (1, 2)}

    done = set()
    if os.path.exists(log):
        done = {json.loads(l)["code"] for l in open(log) if l.strip()}
    codes = pd.read_csv("casf/index.csv").code.astype(str).tolist()

    t0 = time.time()
    with open(log, "a") as fh:
        for code in codes:
            if code in done:
                continue
            rec = {"code": code}
            try:
                fmt, path = LIG_OVERRIDE.get(code, ("sdf", f"casf/lig_sdf/{code}.sdf"))
                receptor = trim_receptor(code, ligand_coords(fmt, path))
                prot = next(oddt.toolkit.readfile("pdb", receptor))
                prot.protein = True
                lig = next(oddt.toolkit.readfile(fmt, path))
                for v, model in models.items():
                    model.set_protein(prot)        # same object -> atom_dict built once
                    rec[f"rfscore_v{v}"] = float(np.ravel(model.predict([lig]))[0])
            except Exception as exc:
                rec["error"] = f"{type(exc).__name__}: {exc}"[:160]
            fh.write(json.dumps(rec) + "\n")
            fh.flush()                              # resumable: kill and rerun any time
    df = pd.DataFrame([json.loads(l) for l in open(log)])
    df.to_csv("casf2016_mlsf_raw.csv", index=False)
    print(f"{len(df)} complexes in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
