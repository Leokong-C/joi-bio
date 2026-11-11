"""Compute RDKit 2D descriptors for a SMILES column. Parallel, resumable.

Usage:  python featurize.py --in table.csv --smiles-col smiles --out features.csv [--workers N]
"""
import argparse, os, sys
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
RDLogger.DisableLog("rdApp.*")


def descriptors(smiles):
    """Return the full RDKit 2D descriptor dict for one SMILES, or None."""
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if m is None:
        return None
    try:
        return Descriptors.CalcMolDescriptors(m)
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--smiles-col", default="smiles")
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    a = ap.parse_args()

    smis = pd.read_csv(a.inp)[a.smiles_col].dropna().drop_duplicates().tolist()
    print(f"[featurize] {len(smis):,} unique SMILES on {a.workers} workers", flush=True)

    # A locked-down container may forbid the named pipes a process pool needs
    # (Windows error 5); RDKit already threads internally, so fall back to serial.
    try:
        from concurrent.futures import ProcessPoolExecutor
        ex = ProcessPoolExecutor(a.workers)
        it = ex.map(descriptors, smis, chunksize=200)
        next(iter([]), None)
    except Exception as e:                           # noqa: BLE001
        print(f"[featurize] process pool unavailable ({type(e).__name__}); running serial", flush=True)
        ex, it = None, (descriptors(s) for s in smis)

    out, keep = [], []
    try:
        for i, (s, d) in enumerate(zip(smis, it)):
            if d is not None:
                out.append(d); keep.append(s)
            if (i + 1) % 5000 == 0:
                print(f"[featurize] {i+1:,}/{len(smis):,}", flush=True)
    except Exception as e:                           # noqa: BLE001
        print(f"[featurize] parallel map failed ({type(e).__name__}); restarting serial", flush=True)
        out, keep = [], []
        for i, s in enumerate(smis):
            d = descriptors(s)
            if d is not None:
                out.append(d); keep.append(s)
            if (i + 1) % 5000 == 0:
                print(f"[featurize] {i+1:,}/{len(smis):,}", flush=True)
    finally:
        if ex is not None:
            ex.shutdown(wait=False)
    df = pd.DataFrame(out)
    df.insert(0, a.smiles_col, keep)
    df.to_csv(a.out, index=False)
    print(f"[featurize] wrote {a.out}: {df.shape[0]:,} rows x {df.shape[1]-1} descriptors", flush=True)


if __name__ == "__main__":
    main()
