#!/usr/bin/env python3
"""
modules.small_molecule.dock — AutoDock Vina 도킹.

리간드: SMILES → ETKDGv3 3D 임베딩 → MMFF 최적화 → meeko PDBQT
수용체: obabel로 PDBQT 변환 (강체, Gasteiger 전하, pH 7.4 수소 부가)
박스:   core.structure.fetch가 정한 center/size

GitHub Actions 6시간 제한을 고려해 --n-shards로 매트릭스 병렬화한다.

알려진 한계 (보고서에 명시할 것)
  - 리지드 수용체. 측쇄 유연성과 유도 적합(induced fit)을 반영하지 않는다.
  - Vina 스코어링은 금속 배위와 PPI 계면을 과소평가한다.
  - 리간드 프로토네이션은 ChEMBL 표기 그대로 — pH 7.4 상태를 계산하지 않는다.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")


def prep_receptor(pdb: Path, pdbqt: Path) -> bool:
    if pdbqt.exists() and pdbqt.stat().st_size > 0:
        return True
    cmd = [
        "obabel", str(pdb), "-O", str(pdbqt),
        "-xr", "-p", "7.4", "--partialcharge", "gasteiger",
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=600)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
        print(f"[!] 수용체 변환 실패: {e}", file=sys.stderr)
        return False
    return pdbqt.exists() and pdbqt.stat().st_size > 0


def ligand_pdbqt(smiles: str, seed: int = 42) -> str | None:
    """SMILES → 3D → PDBQT 문자열."""
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
    except ImportError:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        params.useSmallRingTorsions = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            params.useRandomCoords = True
            if AllChem.EmbedMolecule(mol, params) != 0:
                return None
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=600)
        except Exception:
            pass
        prep = MoleculePreparation()
        setups = prep.prepare(mol)
        if not setups:
            return None
        pdbqt, ok, _err = PDBQTWriterLegacy.write_string(setups[0])
        return pdbqt if ok else None
    except Exception:
        return None


def dock_one(job: tuple) -> dict:
    rec_path, center, size, rec, exhaustiveness = job
    out = {k: rec[k] for k in rec if k != "protected"}
    try:
        lig = ligand_pdbqt(rec["smiles"])
        if lig is None:
            out.update(score=None, error="ligand_prep_failed")
            return out

        from vina import Vina

        v = Vina(sf_name="vina", cpu=1, seed=42, verbosity=0)
        v.set_receptor(rec_path)
        v.set_ligand_from_string(lig)
        v.compute_vina_maps(center=list(center), box_size=list(size))
        v.dock(exhaustiveness=exhaustiveness, n_poses=5)
        energies = v.energies(n_poses=5)
        best = float(energies[0][0])
        out.update(
            score=best,
            pose_spread=round(float(energies[-1][0]) - best, 3),
            ligand_efficiency=None,
        )
        mol = Chem.MolFromSmiles(rec["smiles"])
        if mol is not None and mol.GetNumHeavyAtoms():
            out["ligand_efficiency"] = round(-best / mol.GetNumHeavyAtoms(), 4)
        return out
    except Exception as e:  # noqa: BLE001
        out.update(score=None, error=f"{type(e).__name__}: {e}"[:200])
        return out


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--targetdir", default=None)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--max-ligands", type=int, default=0, help="샤드당 상한, 0=전체")
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--workers", type=int, default=0)
    args = ap.parse_args()

    tdir = Path(args.targetdir) if args.targetdir else root / "targets" / args.gene
    meta = json.loads((tdir / "target_meta.json").read_text())
    pocket = meta["pocket"]

    if pocket["confidence"] == "low":
        print(f"[!] {args.gene}: blind docking (포켓 신뢰도 낮음). 결과 해석 주의.")

    rec_pdbqt = tdir / "receptor.pdbqt"
    if not prep_receptor(tdir / "receptor.pdb", rec_pdbqt):
        return 1

    lib = [
        json.loads(l) for l in (tdir / "library.jsonl").read_text().splitlines() if l.strip()
    ]
    # 벤치마크·활성이 모든 샤드에 고르게 퍼지도록 인덱스 기준으로 자른다
    shard_lib = [r for i, r in enumerate(lib) if i % args.n_shards == args.shard]
    if args.max_ligands:
        shard_lib = shard_lib[: args.max_ligands]

    workers = args.workers or (os.cpu_count() or 2)
    n_act = sum(1 for r in shard_lib if r["label"] == "known_active")
    print(
        f"[*] {args.gene} 도킹: {len(shard_lib)} 리간드 (활성 {n_act}) "
        f"샤드 {args.shard + 1}/{args.n_shards}, 워커 {workers}"
    )
    print(f"    박스 center={pocket['center']} size={pocket['size']} ({pocket['mode']})")

    jobs = [
        (str(rec_pdbqt), pocket["center"], pocket["size"], r, args.exhaustiveness)
        for r in shard_lib
    ]
    results, t0 = [], time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = [ex.submit(dock_one, j) for j in jobs]
        for i, f in enumerate(as_completed(futures), 1):
            results.append(f.result())
            if i % 100 == 0 or i == len(futures):
                el = time.time() - t0
                rate = i / el if el else 0
                eta = (len(futures) - i) / rate if rate else 0
                print(
                    f"    {i}/{len(futures)}  {el / 60:.1f}분 경과, "
                    f"잔여 ~{eta / 60:.1f}분",
                    flush=True,
                )

    scored = [r for r in results if r.get("score") is not None]
    scored.sort(key=lambda r: r["score"])
    out = tdir / f"docked_shard{args.shard}.jsonl"
    out.write_text("\n".join(json.dumps(r, ensure_ascii=False) for r in scored))

    errs: dict[str, int] = {}
    for r in results:
        if r.get("score") is None:
            e = r.get("error", "unknown").split(":")[0]
            errs[e] = errs.get(e, 0) + 1

    print(f"[+] 성공 {len(scored)}/{len(results)} → {out.name}")
    if errs:
        print(f"    실패 내역: {errs}")
    if scored:
        print(f"    최고 스코어: {scored[0]['score']:.2f} kcal/mol ({scored[0]['chembl_id']})")
    fail_rate = 1 - len(scored) / max(1, len(results))
    if fail_rate > 0.10:
        print(f"    [!] 실패율 {fail_rate:.0%} — 10% 초과. 리간드 준비 단계 점검 필요.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
