#!/usr/bin/env python3
"""
사소한 베이스라인 — 우리 벤치마크가 얼마나 쉬운지 잰다.

단백질도 구조도 도킹도 쓰지 않고, 분자에서 원자 하나를 세는 것만으로
순위를 매겼을 때 지표가 얼마나 나오는지 본다.

왜 이게 필요한가
  2026-10-24, 구조를 안 쓰는 리간드 전용 모델이 EF1% 천장(35.95)에 닿고
  AUC 0.99 를 냈다. 도킹(32.682 / 0.7717)을 모든 지표에서 이겼다.
  좋은 소식처럼 보였지만 숫자가 너무 좋아 의심했고, 파보니 모델이 배운 것은
  결합이 아니라 **불소를 세는 법**이었다.

  EPAS1 활성은 PT 계열 — 불소가 많이 붙은 설폰·다이아릴 에테르다.
  decoy 는 "무작위 약물유사 ChEMBL 화합물" 이라 불소가 거의 없다.
  활성 중앙값 3개, decoy 중앙값 0개. 불소 3개 이상이 활성의 80%, decoy 의 4%.

  즉 우리 지표의 상당 부분은 "활성이 도킹하기 좋아서"가 아니라
  "활성이 다른 화합물이라서" 나온 값이다.

계산 비용 0 — 도킹하지 않는다.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")

PROBES = {
    "fluorine_count": ("[F]", "불소 원자 수"),
    "sulfone": ("[SX4](=O)(=O)", "설폰기 수"),
    "nitrile": ("[C]#[N]", "나이트릴기 수"),
    "halogen_any": ("[F,Cl,Br,I]", "할로겐 수"),
}


def rank_metrics(scores: np.ndarray, labels: np.ndarray, seed: int = 0) -> dict:
    """높을수록 좋은 점수로 내림차순 정렬해 지표를 낸다. 동점은 무작위."""
    from core.validation.gate import bedroc, enrichment_factor, roc_auc

    jitter = np.random.RandomState(seed).rand(len(scores)) * 1e-6
    order = np.argsort(-(scores + jitter))
    lab = [int(labels[i]) for i in order]
    return {
        "ef1": enrichment_factor(lab, 0.01),
        "ef5": enrichment_factor(lab, 0.05),
        "auc": roc_auc(lab),
        "bedroc": bedroc(lab),
    }


def main() -> int:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

    ap = argparse.ArgumentParser()
    ap.add_argument("--docked", required=True, help="치료군 docked_all.jsonl")
    ap.add_argument("--control", default=None, help="대조군 docked_all.jsonl")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    rows = [
        json.loads(l)
        for l in Path(args.docked).read_text().splitlines()
        if l.strip()
    ]
    rows = [r for r in rows if r.get("smiles") and r.get("score") is not None]

    mols, keep = [], []
    for r in rows:
        m = Chem.MolFromSmiles(r["smiles"])
        if m is not None:
            mols.append(m)
            keep.append(r)
    y = np.array([1 if r["label"] == "known_active" else 0 for r in keep])
    print(f"[*] {len(y)} 화합물, 활성 {int(y.sum())}개")

    out: dict = {"n": len(y), "n_active": int(y.sum()), "baselines": {}}

    # 도킹 (낮을수록 좋으므로 부호를 뒤집는다)
    dock = rank_metrics(-np.array([r["score"] for r in keep]), y)
    out["docking_treatment"] = dock
    print(f"\n    {'베이스라인':22s} {'EF1%':>8s} {'EF5%':>8s} {'AUC':>8s} {'BEDROC':>8s}")
    print(f"    {'도킹 (치료군)':22s} {dock['ef1']:>8} {dock['ef5']:>8} "
          f"{dock['auc']:>8} {dock['bedroc']:>8}")

    if args.control:
        crows = {
            r["chembl_id"]: r
            for r in (
                json.loads(l)
                for l in Path(args.control).read_text().splitlines()
                if l.strip()
            )
            if r.get("score") is not None
        }
        cs = np.array([crows[r["chembl_id"]]["score"] if r["chembl_id"] in crows
                       else np.nan for r in keep])
        ok = ~np.isnan(cs)
        cd = rank_metrics(-cs[ok], y[ok])
        out["docking_control"] = cd
        print(f"    {'도킹 (대조군)':22s} {cd['ef1']:>8} {cd['ef5']:>8} "
              f"{cd['auc']:>8} {cd['bedroc']:>8}")

    # 작용기 세기
    for name, (smarts, label) in PROBES.items():
        patt = Chem.MolFromSmarts(smarts)
        cnt = np.array([len(m.GetSubstructMatches(patt)) for m in mols], dtype=float)
        mm = rank_metrics(cnt, y)
        out["baselines"][name] = {
            "description": label,
            "median_active": float(np.median(cnt[y == 1])),
            "median_decoy": float(np.median(cnt[y == 0])),
            **mm,
        }
        print(f"    {label:22s} {mm['ef1']:>8} {mm['ef5']:>8} "
              f"{mm['auc']:>8} {mm['bedroc']:>8}")

    # 물성
    for name, fn, label in [
        ("heavy_atoms", lambda m: m.GetNumHeavyAtoms(), "중원자 수"),
        ("tpsa", Descriptors.TPSA, "TPSA"),
        ("mw", Descriptors.MolWt, "분자량"),
    ]:
        v = np.array([float(fn(m)) for m in mols])
        mm = rank_metrics(v, y)
        out["baselines"][name] = {
            "description": label,
            "median_active": float(np.median(v[y == 1])),
            "median_decoy": float(np.median(v[y == 0])),
            **mm,
        }
        print(f"    {label:22s} {mm['ef1']:>8} {mm['ef5']:>8} "
              f"{mm['auc']:>8} {mm['bedroc']:>8}")

    best = max(out["baselines"].items(), key=lambda kv: kv[1]["auc"] or 0)
    out["worst_case"] = {
        "probe": best[0],
        "auc": best[1]["auc"],
        "note": (
            f"'{best[1]['description']}' 하나만 세도 AUC {best[1]['auc']} 가 나온다. "
            f"도킹 치료군은 {dock['auc']}. 단백질도 구조도 쓰지 않은 순위가 "
            "도킹과 비슷하거나 낫다면, 지표의 상당 부분은 decoy 가 활성과 "
            "물성·화학형에서 너무 다르다는 사실에서 나온 것이다."
        ),
    }
    print(f"\n[!] {out['worst_case']['note']}")

    if args.out:
        Path(args.out).write_text(json.dumps(out, ensure_ascii=False, indent=2))
        print(f"\n    → {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
