#!/usr/bin/env python3
"""
modules.small_molecule.admet — 물성·ADMET 프로파일링과 후보 선별.

검증 게이트를 통과한 타깃에 한해 신규 화합물(= known_active가 아닌 것)의
물성을 계산한다. 신장 타깃이므로 신배설(renal clearance) 경향에 주목한다.

여기 계산은 전부 경험적 rule 기반이며 학습된 예측 모델이 아니다. 어떤 화합물도
실험 검증 전에는 '히트'가 아니며 candidate_for_testing으로만 표기한다.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from collections import Counter, defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, QED, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")


def scaffold(smiles: str) -> str | None:
    """Murcko 골격. 곁사슬을 떼고 고리-링커 뼈대만 남긴다."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
    except Exception:
        return None


def scaffold_base_rates(all_rows: list[dict], top_rows: list[dict]) -> list[dict]:
    """골격이 상위권을 차지한 것이 농축인가, 라이브러리 기저율인가.

    왜 이 검사가 필요한가
      2026-10-17 진단. EPAS1 후보 20개 중 10개가 8-아릴-1,3-다이메틸잔틴 하나였다.
      잔틴은 평평하고 단단하며 수소결합기가 많아 작고 닫힌 포켓에 잘 맞고,
      Vina 가 어디서든 높은 점수를 준다 — 문헌에 알려진 promiscuous 골격이다.

      그런데 '상위권의 절반이 한 골격' 자체는 두 가지로 설명된다.
        (1) 그 골격이 이 타깃에 실제로 유리하다 (농축)
        (2) 라이브러리에 원래 그만큼 많다 (기저율)
      ChEMBL 무작위 추출은 의약화학 프로그램의 동족 시리즈로 빽빽하므로 (2) 가
      흔하다. 둘을 구분하지 않으면 '라이브러리가 편중됐다' 를 '히트를 찾았다' 로
      착각한다. 농축비 = (상위권 비율) / (라이브러리 비율) 로 가른다.
    """
    lib = Counter()
    for r in all_rows:
        s = scaffold(r.get("smiles", ""))
        if s:
            lib[s] += 1
    top = Counter()
    for r in top_rows:
        s = r.get("scaffold")
        if s:
            top[s] += 1

    n_lib, n_top = sum(lib.values()) or 1, sum(top.values()) or 1
    out = []
    for s, c in top.most_common():
        lib_c = lib.get(s, 0)
        lib_frac = lib_c / n_lib
        top_frac = c / n_top
        out.append({
            "scaffold": s,
            "n_top": c,
            "n_library": lib_c,
            "top_pct": round(100 * top_frac, 1),
            "library_pct": round(100 * lib_frac, 2),
            "enrichment": round(top_frac / lib_frac, 1) if lib_frac else None,
        })
    return out


def diversify(rows: list[dict], per_scaffold: int, limit: int) -> list[dict]:
    """골격당 상한을 두고 고른다.

    같은 골격에서 상위 몇 개만 남기고 나머지는 다른 골격으로 채운다. 동족
    시리즈 하나가 후보 목록을 독점하는 것을 막는다. 버려진 화합물이 나쁘다는
    뜻은 아니며, 한 시리즈를 깊게 보는 것은 히트가 확정된 뒤의 일이다.
    """
    by: dict[str, list[dict]] = defaultdict(list)
    for r in rows:  # rows 는 이미 composite 내림차순
        by[r.get("scaffold") or "_none"].append(r)

    picked, rank = [], 0
    while len(picked) < limit and rank < per_scaffold:
        added = False
        # 그 순위에서 가장 좋은 골격부터
        cands = [(g[rank], g[rank]["composite_score"]) for g in by.values() if len(g) > rank]
        for r, _ in sorted(cands, key=lambda x: -x[1]):
            if len(picked) >= limit:
                break
            r["scaffold_rank"] = rank + 1
            picked.append(r)
            added = True
        if not added:
            break
        rank += 1
    return picked


def profile(smiles: str) -> dict | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)

    lipinski = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    return {
        "mw": round(mw, 1),
        "clogp": round(logp, 2),
        "tpsa": round(tpsa, 1),
        "hbd": hbd,
        "hba": hba,
        "rotb": rotb,
        "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "fsp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "qed": round(QED.qed(mol), 3),
        "lipinski_violations": lipinski,
        "veber_pass": rotb <= 10 and tpsa <= 140,
        # 신배설 경향: 저분자 + 저친유성 + 고극성
        "renal_clearance_favored": mw < 400 and logp < 2.0 and tpsa > 75,
        # 신장 타깃은 CNS 노출이 불필요하므로 뇌 투과가 낮은 편이 유리
        "bbb_penetration_unlikely": tpsa > 90 or hbd > 3,
    }


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--targetdir", default=None)
    ap.add_argument("--top-n", type=int, default=300)
    ap.add_argument("--per-scaffold", type=int, default=2,
                    help="한 골격에서 올릴 최대 후보 수")
    ap.add_argument("--report-n", type=int, default=20,
                    help="골격 다양화를 적용해 추릴 최종 후보 수")
    ap.add_argument("--min-heavy-frac", type=float, default=0.75,
                    help="알려진 활성 중원자 중앙값 대비 하한. 작고 평평한 분자가 "
                         "닫힌 포켓에서 과대평가되는 것을 막는다. 0 이면 끈다")
    args = ap.parse_args()

    tdir = Path(args.targetdir) if args.targetdir else root / "targets" / args.gene
    val = json.loads((tdir / "validation.json").read_text())
    rows = [
        json.loads(l)
        for l in (tdir / "docked_all.jsonl").read_text().splitlines()
        if l.strip()
    ]

    # 알려진 활성의 크기를 기준선으로 — 이 포켓은 작고 닫혀 있어 큰 분자가
    # 아예 들어가지 못하고(EPAS1 라이브러리의 22.6%), 그 결과 상위권이 활성보다
    # 작은 분자로 채워진다(중원자 20 대 27). 활성보다 지나치게 작은 것은 거른다.
    act_heavy = sorted(
        Chem.MolFromSmiles(r["smiles"]).GetNumHeavyAtoms()
        for r in rows
        if r.get("label") == "known_active" and Chem.MolFromSmiles(r.get("smiles", ""))
    )
    heavy_floor = 0
    if act_heavy and args.min_heavy_frac > 0:
        med = act_heavy[len(act_heavy) // 2]
        heavy_floor = int(med * args.min_heavy_frac)

    novel = [r for r in rows if r.get("label") != "known_active"][: args.top_n]
    out, n_small = [], 0
    for r in novel:
        p = profile(r["smiles"])
        if p is None:
            continue
        mol = Chem.MolFromSmiles(r["smiles"])
        nheavy = mol.GetNumHeavyAtoms()
        if heavy_floor and nheavy < heavy_floor:
            n_small += 1
            continue
        # 종합 점수: 도킹 스코어 0.7 + 약물유사성 0.3. 가중치는 임의이며
        # 실험 데이터가 쌓이면 재보정해야 한다.
        composite = (-r["score"]) * 0.7 + p["qed"] * 10 * 0.3
        out.append(
            {
                **r,
                "admet": p,
                "n_heavy": nheavy,
                "ligand_efficiency": round(-r["score"] / nheavy, 3) if nheavy else None,
                "scaffold": scaffold(r["smiles"]),
                "composite_score": round(composite, 3),
                "status": "candidate_for_testing",
                "validation_verdict": val["verdict"],
            }
        )

    out.sort(key=lambda r: -r["composite_score"])
    (tdir / "candidates.jsonl").write_text(
        "\n".join(json.dumps(r, ensure_ascii=False) for r in out)
    )

    # 골격 다양화 — 동족 시리즈 하나가 후보 목록을 독점하는 것을 막는다
    picked = diversify(out, args.per_scaffold, args.report_n)
    naive = out[: args.report_n]
    rates = scaffold_base_rates(rows, naive)
    report = {
        "gene": args.gene,
        "validation_verdict": val["verdict"],
        "heavy_atom_floor": heavy_floor,
        "n_rejected_too_small": n_small,
        "n_profiled": len(out),
        "per_scaffold_cap": args.per_scaffold,
        "naive_top": {
            "n": len(naive),
            "n_scaffolds": len({r["scaffold"] for r in naive if r["scaffold"]}),
            "scaffolds": rates[:8],
        },
        "diversified": {
            "n": len(picked),
            "n_scaffolds": len({r["scaffold"] for r in picked if r["scaffold"]}),
        },
        "compounds": [
            {k: r[k] for k in ("chembl_id", "score", "composite_score", "n_heavy",
                               "ligand_efficiency", "scaffold", "scaffold_rank", "smiles")
             if k in r}
            for r in picked
        ],
    }
    (tdir / "scaffold_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False)
    )

    print(f"[*] {args.gene}: 후보 {len(out)}개 프로파일링 (검증={val['verdict']})")
    if val["verdict"] != "PASS":
        print("    [!] 검증 미통과 — 이 후보군은 보고서에 히트로 실리지 않는다.")
    if heavy_floor:
        print(f"    중원자 하한 {heavy_floor} (활성 중앙값의 {args.min_heavy_frac:.0%}) "
              f"— {n_small}개 제외")
    print()
    print(f"[*] 골격 기저율 검사 — 단순 상위 {len(naive)}개의 골격 구성")
    print(f"    서로 다른 골격 {report['naive_top']['n_scaffolds']}개")
    print(f"    {'상위':>5s} {'라이브러리':>9s} {'농축비':>6s}  골격")
    for s in rates[:6]:
        e = f"{s['enrichment']}x" if s["enrichment"] else "-"
        flag = ""
        if s["enrichment"] and s["enrichment"] < 2.0 and s["n_top"] >= 3:
            flag = "  ← 기저율. 라이브러리에 원래 많다"
        print(f"    {s['n_top']:3d}개 {s['top_pct']:4.0f}%  {s['library_pct']:7.2f}% "
              f"{e:>6s}  {s['scaffold'][:40]}{flag}")
    print()
    print(f"[*] 골격당 최대 {args.per_scaffold}개로 다양화 → "
          f"{report['diversified']['n']}개 / 골격 {report['diversified']['n_scaffolds']}종")
    for r in picked[:8]:
        a = r["admet"]
        print(f"    {r['chembl_id']:16s} vina={r['score']:6.2f} "
              f"LE={r['ligand_efficiency']} QED={a['qed']:.2f} "
              f"중원자={r['n_heavy']:2d} 골격순위={r.get('scaffold_rank')}")
    print()
    print(f"[+] {tdir / 'scaffold_report.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
