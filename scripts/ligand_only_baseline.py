#!/usr/bin/env python3
"""
구조를 전혀 쓰지 않는 베이스라인 — 도킹이 구조로 무엇을 벌었는지 재는 대조.

왜 이걸 먼저 돌리는가
  EPAS1 편향 측정에서 상위권 순위의 75%가 구조 선택에서 나온다는 결과가 나왔다
  (치료군 EF1% 32.682, 대조군 8.170). 그렇다면 구조를 아예 안 쓰면 어떻게 되는가?

  리간드 지문만으로 학습한 모델은 수용체 구조를 한 번도 보지 않으므로
  **구조 선택 편향이 있을 수 없다.** 이 모델이 도킹과 비슷하거나 더 낫다면,
  그 구조가 사실상 아무것도 사주지 않았다는 뜻이다.

  ESM 을 쓰지 않는 이유: 타깃이 하나면 단백질 임베딩이 모든 행에서 같은
  상수 벡터라 모델이 쓸 수가 없다. ESM 은 여러 타깃을 가로질러 일반화할 때
  비로소 값을 한다 (scripts/README 의 2단계 참조).

평가 설계
  Murcko 골격 기준 교차검증(cold-scaffold). 같은 골격이 학습과 평가에 동시에
  들어가면 암기로 점수가 오른다 — Wallach & Heifets 가 리간드 기반 벤치마크의
  고질병으로 지목한 것이고, 우리가 이미 잔틴 사건에서 본 현상이다.
  out-of-fold 예측만 모아 도킹과 같은 지표로 잰다.

계산 비용 0 — 도킹하지 않는다.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.linear_model import LogisticRegression

RDLogger.DisableLog("rdApp.*")


def scaffold(smi: str) -> str:
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return ""
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=m)
    except Exception:
        return ""


def main() -> int:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from core.validation.gate import bedroc, enrichment_factor, roc_auc

    ap = argparse.ArgumentParser()
    ap.add_argument("--docked", required=True)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--radius", type=int, default=2)
    ap.add_argument("--bits", type=int, default=2048)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    rows = [json.loads(l) for l in Path(args.docked).read_text().splitlines() if l.strip()]
    rows = [r for r in rows if r.get("smiles") and r.get("score") is not None]

    gen = rdFingerprintGenerator.GetMorganGenerator(
        radius=args.radius, fpSize=args.bits
    )
    X, y, scaf, keep = [], [], [], []
    for r in rows:
        m = Chem.MolFromSmiles(r["smiles"])
        if m is None:
            continue
        X.append(np.array(gen.GetFingerprint(m), dtype=np.int8))
        y.append(1 if r.get("label") == "known_active" else 0)
        scaf.append(scaffold(r["smiles"]))
        keep.append(r)
    X = np.array(X)
    y = np.array(y)
    print(f"[*] {len(y)} 화합물, 활성 {int(y.sum())}개, 골격 {len(set(scaf))}종")

    # 골격을 통째로 폴드에 배정 — 같은 골격이 학습·평가에 걸치지 않게
    by_scaf: dict[str, list[int]] = defaultdict(list)
    for i, s in enumerate(scaf):
        by_scaf[s].append(i)
    groups = sorted(by_scaf.values(), key=len, reverse=True)
    folds: list[list[int]] = [[] for _ in range(args.folds)]
    for g in groups:  # 큰 골격부터 가장 작은 폴드에
        folds.sort(key=len)
        folds[0].extend(g)

    oof = np.zeros(len(y), dtype=float)
    for k, test_idx in enumerate(folds):
        train_idx = [i for i in range(len(y)) if i not in set(test_idx)]
        if len(set(y[train_idx])) < 2:
            print(f"    fold {k}: 학습셋에 한 클래스뿐 — 건너뜀")
            continue
        clf = LogisticRegression(max_iter=2000, C=1.0, class_weight="balanced")
        clf.fit(X[train_idx], y[train_idx])
        oof[test_idx] = clf.predict_proba(X[test_idx])[:, 1]
        n_act = int(y[test_idx].sum())
        print(f"    fold {k}: 평가 {len(test_idx):4d}개 (활성 {n_act:2d})")

    order = np.argsort(-oof)  # 확률 높은 순
    labels = [int(y[i]) for i in order]
    ml = {
        "ef1": enrichment_factor(labels, 0.01),
        "ef5": enrichment_factor(labels, 0.05),
        "auc": roc_auc(labels),
        "bedroc": bedroc(labels),
    }

    dock_order = sorted(range(len(keep)), key=lambda i: keep[i]["score"])
    dl = [int(y[i]) for i in dock_order]
    dock = {
        "ef1": enrichment_factor(dl, 0.01),
        "ef5": enrichment_factor(dl, 0.05),
        "auc": roc_auc(dl),
        "bedroc": bedroc(dl),
    }

    print("\n    지표        리간드만(구조 안 씀)   도킹(같은 라이브러리)")
    for k in ("ef1", "ef5", "auc", "bedroc"):
        print(f"    {k:8s}  {ml[k]:>14}   {dock[k]:>14}")

    out = {
        "n": len(y),
        "n_active": int(y.sum()),
        "n_scaffolds": len(set(scaf)),
        "folds": args.folds,
        "split": "cold_scaffold",
        "ligand_only": ml,
        "docking_same_library": dock,
        "note": (
            "리간드 지문만 쓴 모델은 수용체 구조를 보지 않으므로 구조 선택 편향이 "
            "있을 수 없다. 다만 ChEMBL 유사체 암기라는 자기 편향이 있어 "
            "골격 분할로 완화했을 뿐 제거하지는 못했다."
        ),
    }
    if args.out:
        Path(args.out).write_text(json.dumps(out, ensure_ascii=False, indent=2))
        print(f"\n    → {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
