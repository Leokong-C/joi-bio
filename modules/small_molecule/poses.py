#!/usr/bin/env python3
"""
modules.small_molecule.poses — 선택한 화합물의 도킹 포즈 좌표를 저장한다.

Stage 03은 스코어만 남기고 포즈를 버린다. 2,305개 포즈를 전부 보관하면 용량이
커지고 대부분은 볼 일이 없기 때문이다. 이 단계는 실제로 들여다볼 화합물만
골라 다시 도킹하고 좌표를 SDF로 남긴다.

고르는 기준
  1. 벤치마크 리간드 (belzutifan 등) — 검증의 기준점이므로 반드시 포함
  2. 도킹 상위 N개 — 눈으로 확인할 후보

출력
  poses.sdf       리간드 포즈 (SDF, 각 분자에 스코어·이름 태그)
  receptor.pdb    수용체 (Stage 01 결과를 그대로 복사)
  poses_meta.json 뷰어가 읽을 메타데이터 (포켓 박스, 화합물 목록)
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy, RDKitMolCreate
from vina import Vina

from modules.small_molecule.dock import prep_receptor

RDLogger.DisableLog("rdApp.*")


def dock_with_pose(
    rec_pdbqt: str, center: list, size: list, smiles: str, exhaustiveness: int, n_poses: int
) -> tuple[list, list]:
    """도킹 후 포즈를 RDKit mol 리스트로 돌려준다."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [], []
    mol = Chem.AddHs(mol)
    p = AllChem.ETKDGv3()
    p.randomSeed = 42
    if AllChem.EmbedMolecule(mol, p) != 0:
        p.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, p) != 0:
            return [], []
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=600)
    except Exception:
        pass

    setups = MoleculePreparation().prepare(mol)
    if not setups:
        return [], []
    lig_pdbqt, ok, _ = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        return [], []

    v = Vina(sf_name="vina", cpu=0, seed=42, verbosity=0)
    v.set_receptor(rec_pdbqt)
    v.set_ligand_from_string(lig_pdbqt)
    v.compute_vina_maps(center=list(center), box_size=list(size))
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    energies = [float(e[0]) for e in v.energies(n_poses=n_poses)]

    with tempfile.NamedTemporaryFile("w+", suffix=".pdbqt", delete=False) as fh:
        tmp = Path(fh.name)
    try:
        v.write_poses(str(tmp), n_poses=n_poses, overwrite=True)
        pmol = PDBQTMolecule.from_file(str(tmp), skip_typing=True)
        mols = RDKitMolCreate.from_pdbqt_mol(pmol)
        mols = [m for m in mols if m is not None]
    finally:
        tmp.unlink(missing_ok=True)
    return mols, energies


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--targetdir", default=None)
    ap.add_argument("--top-n", type=int, default=15, help="상위 후보 몇 개를 저장할지")
    ap.add_argument("--poses-per-ligand", type=int, default=3)
    ap.add_argument("--exhaustiveness", type=int, default=16, help="포즈용은 더 꼼꼼히")
    args = ap.parse_args()

    tdir = Path(args.targetdir) if args.targetdir else root / "targets" / args.gene
    meta = json.loads((tdir / "target_meta.json").read_text())
    pocket = meta["pocket"]
    rec_pdbqt = tdir / "receptor.pdbqt"
    if not rec_pdbqt.exists() and not prep_receptor(tdir / "receptor.pdb", rec_pdbqt):
        print("[!] 수용체 PDBQT 변환 실패", file=sys.stderr)
        return 1

    rows = [
        json.loads(l)
        for l in (tdir / "docked_all.jsonl").read_text().splitlines()
        if l.strip()
    ]
    bench_names = {b.lower() for b in meta.get("benchmark_ligands", [])}

    selected, seen = [], set()
    # 1) 벤치마크 — 검증의 기준점
    for r in rows:
        label = f"{r.get('pref_name','')} {r.get('query_name','')}".lower()
        if any(b in label for b in bench_names) and r["chembl_id"] not in seen:
            r["_group"] = "benchmark"
            selected.append(r)
            seen.add(r["chembl_id"])
    # 2) 상위 후보
    # selected 와 seen 은 항상 같이 늘어나므로 len(selected) >= len(seen) + top_n 은
    # 영원히 참이 되지 않는다. 그래서 2026-10-16 실행에서 15개가 아니라 142개를
    # 돌다 타임아웃이 났다. 벤치마크 수를 따로 세어 기준으로 삼는다.
    n_bench = len(selected)
    for r in rows:
        if len(selected) >= n_bench + args.top_n:
            break
        if r["chembl_id"] in seen:
            continue
        r["_group"] = "top_hit"
        selected.append(r)
        seen.add(r["chembl_id"])

    print(f"[*] {args.gene}: {len(selected)}개 화합물의 포즈 생성")
    writer = Chem.SDWriter(str(tdir / "poses.sdf"))
    manifest = []
    for i, r in enumerate(selected, 1):
        name = r.get("pref_name") or r["chembl_id"]
        mols, energies = dock_with_pose(
            str(rec_pdbqt),
            pocket["center"],
            pocket["size"],
            r["smiles"],
            args.exhaustiveness,
            args.poses_per_ligand,
        )
        if not mols:
            print(f"    [{i}/{len(selected)}] {name}: 포즈 생성 실패")
            continue
        for j, m in enumerate(mols[: args.poses_per_ligand]):
            e = energies[j] if j < len(energies) else None
            m.SetProp("_Name", f"{name}_pose{j + 1}")
            m.SetProp("chembl_id", r["chembl_id"])
            m.SetProp("group", r["_group"])
            m.SetProp("pose_rank", str(j + 1))
            if e is not None:
                m.SetProp("vina_score", f"{e:.2f}")
            writer.write(m)
        manifest.append(
            {
                "name": name,
                "chembl_id": r["chembl_id"],
                "group": r["_group"],
                "screen_score": r.get("score"),
                "pose_scores": [round(e, 2) for e in energies[: args.poses_per_ligand]],
                "n_poses": len(mols[: args.poses_per_ligand]),
                "smiles": r["smiles"],
            }
        )
        best = energies[0] if energies else float("nan")
        print(f"    [{i}/{len(selected)}] {name}: {len(mols)} 포즈, 최고 {best:.2f}")
    writer.close()

    (tdir / "poses_meta.json").write_text(
        json.dumps(
            {
                "gene": args.gene,
                "pdb_id": meta.get("pdb_id"),
                "resolution": meta.get("resolution"),
                "pocket": pocket,
                "reference_ligand": pocket.get("reference_ligand"),
                "compounds": manifest,
            },
            indent=2,
            ensure_ascii=False,
        )
    )
    print(f"[+] poses.sdf ({len(manifest)} 화합물) + poses_meta.json 생성")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
