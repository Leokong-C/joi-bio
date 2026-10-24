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
    ap.add_argument("--docked", default=None, help="치료군 docked_all.jsonl")
    ap.add_argument(
        "--library",
        default=None,
        help="도킹 전 library.jsonl. 스코어 없이 라벨만으로 벤치마크 난이도를 잰다. "
        "돈을 쓰기 전에 라이브러리가 쓸만한지 확인하는 용도다",
    )
    ap.add_argument("--control", default=None, help="대조군 docked_all.jsonl")
    ap.add_argument(
        "--max-auc",
        type=float,
        default=0.60,
        help="수용 기준. 단일 기술자 최대 AUC 가 이보다 크면 비정상 종료한다",
    )
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if not args.docked and not args.library:
        print("[!] --docked 또는 --library 중 하나는 있어야 한다")
        return 2

    src = args.docked or args.library
    rows = [json.loads(l) for l in Path(src).read_text().splitlines() if l.strip()]
    rows = [r for r in rows if r.get("smiles")]
    if args.docked:
        rows = [r for r in rows if r.get("score") is not None]
    # 스크리닝 전용 화합물은 채점 대상이 아니다
    rows = [r for r in rows if r.get("label") != "screen"]

    mols, keep = [], []
    for r in rows:
        m = Chem.MolFromSmiles(r["smiles"])
        if m is not None:
            mols.append(m)
            keep.append(r)
    y = np.array([1 if r["label"] == "known_active" else 0 for r in keep])
    print(f"[*] {len(y)} 화합물, 활성 {int(y.sum())}개")

    out: dict = {"n": len(y), "n_active": int(y.sum()), "baselines": {}}

    print(f"\n    {'베이스라인':22s} {'EF1%':>8s} {'EF5%':>8s} {'AUC':>8s} {'BEDROC':>8s}")
    dock = None
    if args.docked:
        # 도킹 (낮을수록 좋으므로 부호를 뒤집는다)
        dock = rank_metrics(-np.array([r["score"] for r in keep]), y)
        out["docking_treatment"] = dock
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

    # 진짜 수용 기준 — 사소한 기술자가 우리 게이트를 통과하는가.
    #
    # 임의의 AUC 문턱보다 이쪽이 낫다. 게이트 문턱(EF1 5.0 / AUC 0.70 /
    # BEDROC 0.30)은 이미 이 프로젝트가 쓰는 값이고, 단백질도 구조도 안 쓴
    # 순위가 그걸 통과한다면 그 벤치마크로는 아무것도 증명할 수 없다.
    #
    # 기존 ChEMBL 무작위 라이브러리는 불소 세기만으로 3/3 을 통과했다.
    GATE = {"ef1": 5.0, "auc": 0.70, "bedroc": 0.30}
    mx = {k: max((v[k] or 0) for v in out["baselines"].values()) for k in GATE}
    gate_hits = [k for k, t in GATE.items() if mx[k] >= t]
    out["gate_check"] = {
        "thresholds": GATE,
        "max_over_probes": mx,
        "passed": gate_hits,
        "n_passed": len(gate_hits),
        "trivial_would_pass_gate": len(gate_hits) >= 2,
    }
    print(f"\n[*] 사소한 기술자 최고치 — EF1% {mx['ef1']:.3f} / AUC {mx['auc']:.4f} "
          f"/ BEDROC {mx['bedroc']:.4f}")
    if gate_hits:
        print(f"    게이트 문턱을 넘는 것: {', '.join(gate_hits)} ({len(gate_hits)}/3)")
    else:
        print("    게이트 문턱을 넘는 것이 하나도 없다 — 벤치마크가 제 구실을 한다")

    best = max(out["baselines"].items(), key=lambda kv: kv[1]["auc"] or 0)
    worst_auc = best[1]["auc"] or 0
    note = (
        f"'{best[1]['description']}' 하나만 세도 AUC {worst_auc} 가 나온다. "
        "단백질도 구조도 쓰지 않은 순위가 이 정도라면, 지표의 상당 부분은 "
        "decoy 가 활성과 물성·화학형에서 너무 다르다는 사실에서 나온 것이다."
    )
    if dock:
        note += f" 참고로 도킹 치료군은 {dock['auc']}."
    # 게이트를 2개 이상 넘으면 무조건 탈락. 아니면 보조 기준(AUC)으로 본다.
    passed = (not out["gate_check"]["trivial_would_pass_gate"]) and worst_auc <= args.max_auc
    out["worst_case"] = {
        "probe": best[0],
        "auc": worst_auc,
        "threshold": args.max_auc,
        "accepted": passed,
        "gate_check": out["gate_check"],
        "note": note,
    }
    print(f"\n[{'+' if passed else '!'}] 단일 기술자 최대 AUC {worst_auc} "
          f"(기준 {args.max_auc}) — {'수용' if passed else '탈락'}")
    print(f"    {note}")

    if args.out:
        Path(args.out).write_text(json.dumps(out, ensure_ascii=False, indent=2))
        print(f"\n    → {args.out}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
