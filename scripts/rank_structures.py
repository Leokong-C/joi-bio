#!/usr/bin/env python3
"""
한 타깃의 홀로 구조를 **파이프라인 자신의 선택 기준으로** 전부 순위 매긴다.

## 왜 이게 필요한가

우리 간판 측정은 "구조 선택이 상위권 enrichment 의 얼마를 움직이는가"다.
그 측정이 성립하려면 **치료군은 우리 기준이 고른 구조**, **대조군은 같은
기준에서 가장 낮은 구조**여야 한다. 내가 눈으로 고르면 측정이 아니라 연출이다.

`select_structure` 는 충분히 닮은 구조를 찾으면 거기서 멈춘다(실행 시간을
아끼려고). 양팔을 정하려면 **멈추지 않고 전부** 봐야 한다. 이 스크립트가
그 역할이다.

선택 기준은 `core/structure/fetch.py` 그대로다 — 벤치마크 리간드와의 유사도로
구간을 나누고, 같은 구간 안에서 해상도로 고른다.

RCSB·ChEMBL 이 Claude 컨테이너에서 프록시 403 이라 Actions 에서만 돈다.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main() -> int:
    import yaml

    from core.structure.fetch import (evaluate_candidate, score_candidate,
                                      search_rcsb)
    from modules.small_molecule.library import fetch_benchmark_ligands

    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--max-download", type=int, default=30)
    ap.add_argument("--out", default="structure_ranking.json")
    args = ap.parse_args()

    cfg = yaml.safe_load((ROOT / "config" / "targets.yaml").read_text())
    tgt = next((t for t in cfg["targets"] if t["gene"] == args.gene), None)
    if tgt is None:
        print(f"[!] config/targets.yaml 에 {args.gene} 가 없다")
        return 1

    names = tgt.get("benchmark_ligands", [])
    print(f"[1] 벤치마크 리간드 {names}")
    bench = fetch_benchmark_ligands(names)
    bench_smiles = [b["smiles"] for b in bench]
    for b in bench:
        print(f"    {b['query_name']:16s} → {b['chembl_id']}  {b['smiles'][:60]}")
    if not bench_smiles:
        print("    [!] 하나도 못 찾았다 — 유사도가 전부 0 이 되어 순위가 해상도만 "
              "반영한다. 이름을 고치거나 SMILES 를 직접 넣을 것")

    print(f"\n[2] {tgt['uniprot']} 구조 전수 평가 (최대 {args.max_download}개)")
    ids = search_rcsb(tgt["uniprot"])
    print(f"    RCSB 후보 {len(ids)}개")

    holo, apo, failed = [], [], []
    for pid in ids[: args.max_download]:
        c = evaluate_candidate(pid)
        if c is None:
            failed.append(pid)
            continue
        if not c["has_ligand"]:
            apo.append({"pdb_id": pid, "resolution": c["resolution"]})
            continue
        m = score_candidate(c, bench_smiles)
        lig = c["ligands"][0]
        holo.append({"pdb_id": pid, "resolution": c["resolution"],
                     "ligand": lig["resname"], "n_heavy": lig["n_heavy"],
                     "similarity": round(m["similarity"], 4),
                     "ligand_smiles": m["ligand_smiles"]})
        print(f"      {pid}  {lig['resname']:5s} heavy {lig['n_heavy']:3d}  "
              f"res {c['resolution']}  유사도 {m['similarity']:.3f}")

    if not holo:
        print("    [!] 홀로 구조가 없다")
        Path(args.out).write_text(json.dumps(
            {"gene": args.gene, "holo": [], "apo": apo}, ensure_ascii=False, indent=2))
        return 1

    def bucket(c: dict) -> int:          # fetch.select_structure 와 동일
        s = c["similarity"]
        return 0 if s >= 0.45 else 1 if s >= 0.25 else 2

    holo.sort(key=lambda c: (bucket(c), round(c["resolution"] or 99.0, 1),
                             -c["similarity"]))

    print(f"\n[3] 순위 — 파이프라인 선택 기준 (유사도 구간 → 해상도 → 유사도)")
    print(f"    {'순위':>4s} {'PDB':6s} {'리간드':6s} {'해상도':>7s} {'유사도':>7s} {'구간':>5s}")
    for i, c in enumerate(holo, 1):
        mark = "  ← 치료군" if i == 1 else ("  ← 대조군 후보" if i == len(holo) else "")
        print(f"    {i:4d} {c['pdb_id']:6s} {c['ligand']:6s} "
              f"{str(c['resolution']):>7s} {c['similarity']:7.3f} "
              f"{bucket(c):5d}{mark}")

    treat, ctrl = holo[0], holo[-1]
    out = {"gene": args.gene, "uniprot": tgt["uniprot"],
           "benchmark_ligands": names,
           "benchmark_resolved": [b["chembl_id"] for b in bench],
           "n_rcsb": len(ids), "n_holo": len(holo), "n_apo": len(apo),
           "failed": failed, "holo_ranked": holo, "apo": apo,
           "treatment": treat, "control": ctrl,
           "similarity_gap": round(treat["similarity"] - ctrl["similarity"], 4)}
    Path(args.out).write_text(json.dumps(out, ensure_ascii=False, indent=2))

    print(f"\n[4] 양팔 제안")
    print(f"    치료군  {treat['pdb_id']}  {treat['ligand']}  "
          f"res {treat['resolution']}  유사도 {treat['similarity']:.3f}")
    print(f"    대조군  {ctrl['pdb_id']}  {ctrl['ligand']}  "
          f"res {ctrl['resolution']}  유사도 {ctrl['similarity']:.3f}")
    print(f"    유사도 차이 {out['similarity_gap']:.3f}")
    if out["similarity_gap"] < 0.1:
        print(f"    [!] 차이가 작다. 구조가 균질하면 구조 선택 편향도 작게 나온다 —")
        print(f"        그 자체가 결과다. EPAS1(0.400 대 0.071)과 다를 수 있다")
    print(f"    → {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
