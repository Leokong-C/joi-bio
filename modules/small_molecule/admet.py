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

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, QED, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")


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
    args = ap.parse_args()

    tdir = Path(args.targetdir) if args.targetdir else root / "targets" / args.gene
    val = json.loads((tdir / "validation.json").read_text())
    rows = [
        json.loads(l)
        for l in (tdir / "docked_all.jsonl").read_text().splitlines()
        if l.strip()
    ]

    novel = [r for r in rows if r.get("label") != "known_active"][: args.top_n]
    out = []
    for r in novel:
        p = profile(r["smiles"])
        if p is None:
            continue
        # 종합 점수: 도킹 스코어 0.7 + 약물유사성 0.3. 가중치는 임의이며
        # 실험 데이터가 쌓이면 재보정해야 한다.
        composite = (-r["score"]) * 0.7 + p["qed"] * 10 * 0.3
        out.append(
            {
                **r,
                "admet": p,
                "composite_score": round(composite, 3),
                "status": "candidate_for_testing",
                "validation_verdict": val["verdict"],
            }
        )

    out.sort(key=lambda r: -r["composite_score"])
    (tdir / "candidates.jsonl").write_text(
        "\n".join(json.dumps(r, ensure_ascii=False) for r in out)
    )

    print(f"[*] {args.gene}: 후보 {len(out)}개 프로파일링 (검증={val['verdict']})")
    if val["verdict"] != "PASS":
        print("    [!] 검증 미통과 — 이 후보군은 보고서에 히트로 실리지 않는다.")
    for r in out[:5]:
        a = r["admet"]
        print(
            f"    {r['chembl_id']:16s} vina={r['score']:6.2f} LE={r.get('ligand_efficiency')} "
            f"QED={a['qed']:.2f} MW={a['mw']:.0f} cLogP={a['clogp']:.1f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
