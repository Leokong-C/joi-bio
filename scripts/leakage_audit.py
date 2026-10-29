#!/usr/bin/env python3
"""
누수 감사 — 라이브러리가 기억만으로 풀리는지 본다.

## 왜 이 검사인가

`trivial_baselines.py` 는 **단일 기술자**(불소 수, 설폰 수, 분자량…)로 순위를
매겨 봤다. 2026-10-24 에 그걸로 불소 누수를 잡았고, PubChem 매칭으로 고쳤다.

그런데 그 검사가 못 잡는 누수가 하나 더 있다. **활성끼리 서로 닮아 있으면**,
"다른 활성과 얼마나 닮았나" 만 보고도 활성을 찾아낼 수 있다. 단백질도, 구조도,
도킹도, 심지어 어떤 기술자가 중요한지도 몰라도 된다. 그냥 외우면 된다.

SieveStack 이 2025년 LIT-PCBA 감사(arXiv:2507.21404)에서 이걸로 그 벤치마크를
깼다. 학습 없는 ECFP4 최대유사도 베이스라인이 중위 EF1% 4.15 를 냈고, 같은
조건의 CHEESE 3D 인코더 중위값과 같았다. 정교한 방법이 "외우기" 를 못 이긴 것이다.

우리 EPAS1 활성은 PT 계열 한 화학형이다. **자기유사도가 높을 가능성이 크고,
그러면 PubChem 으로 decoy 를 고쳐도 벤치마크는 여전히 기억을 보상한다.**

## 무엇을 재는가

1. **최대유사도 베이스라인 (핵심)** — 각 화합물을 "자기 자신을 뺀 나머지 활성과의
   최대 Tanimoto" 로 순위 매긴다. 활성에 대해 leave-one-out 이라 자기 자신을
   맞히는 반칙이 없다. 이 순위가 게이트를 통과하면 벤치마크가 기억을 보상한다.
2. **동일 분자 중복** — 정규 SMILES 가 같은 화합물. 라벨이 엇갈리면(활성이면서
   decoy) 라벨링 오류다.
3. **유사쌍 누수** — 활성-decoy 쌍 중 Tanimoto >= 0.6. `decoys.select` 가
   0.35 로 자르므로 이건 그 상한이 실제로 지켜졌는지 확인하는 검사다.
4. **활성 자기유사도 분포** — 활성군이 한 골격에 뭉쳐 있는가.

지문은 SieveStack 과 같은 ECFP4(Morgan radius=2, 4096 bits)를 쓴다. 숫자를
그들 것과 나란히 놓으려면 같은 지문이어야 한다.

비용 0 — 도킹하지 않는다.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")

# SieveStack 과 동일 — 숫자를 비교하려면 지문이 같아야 한다
FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=4096)

ANALOG_T = 0.6  # SieveStack 이 유사쌍 누수로 센 문턱
DECOY_CAP = 0.35  # decoys.select 가 거는 상한


def load(path: Path) -> list[dict]:
    rows = [json.loads(l) for l in path.read_text().splitlines() if l.strip()]
    out = []
    for r in rows:
        if not r.get("smiles") or r.get("label") == "screen":
            continue
        m = Chem.MolFromSmiles(r["smiles"])
        if m is None:
            continue
        r["_mol"] = m
        r["_fp"] = FPGEN.GetFingerprint(m)
        r["_canon"] = Chem.MolToSmiles(m)
        out.append(r)
    return out


def max_sim_baseline(rows: list[dict]) -> dict:
    """각 화합물 점수 = 자기 자신을 뺀 활성과의 최대 Tanimoto."""
    act_idx = [i for i, r in enumerate(rows) if r["label"] == "known_active"]
    act_fps = [rows[i]["_fp"] for i in act_idx]
    act_pos = {i: j for j, i in enumerate(act_idx)}

    scores = np.zeros(len(rows))
    for i, r in enumerate(rows):
        sims = DataStructs.BulkTanimotoSimilarity(r["_fp"], act_fps)
        if i in act_pos:
            sims[act_pos[i]] = -1.0  # leave-one-out — 자기 자신 제외
        scores[i] = max(sims) if sims else 0.0
    return {"scores": scores,
            "labels": np.array([1 if r["label"] == "known_active" else 0 for r in rows])}


def scaffold_of(mol) -> str:
    from rdkit.Chem.Scaffolds import MurckoScaffold
    try:
        return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
    except Exception:
        return ""


def max_sim_baseline_scaffold_split(rows: list[dict]) -> dict:
    """같은 Murcko 골격의 활성은 참조에서 뺀다.

    왜 이 변형이 결정적인가. 위 leave-one-out 은 자기 자신만 뺀다. 그래서
    같은 시리즈의 유사체가 참조에 남아 있고, "내 쌍둥이를 찾기" 로 풀린다.
    골격째로 빼면 **다른 골격의 활성을 보고 이 활성을 알아볼 수 있는가** 를
    묻게 된다 — 그게 일반화고, 앞의 것은 기억이다.

    여기서도 통과하면 활성군 자체가 벤치마크로 쓸 수 없다는 뜻이고,
    여기서 무너지면 골격 분할 평가로 고칠 수 있다는 뜻이다.
    """
    for r in rows:
        if "_scaf" not in r:
            r["_scaf"] = scaffold_of(r["_mol"])
    acts = [r for r in rows if r["label"] == "known_active"]
    scores = np.zeros(len(rows))
    for i, r in enumerate(rows):
        pool = [a["_fp"] for a in acts if a["_scaf"] != r["_scaf"]]
        if not pool:
            scores[i] = 0.0
            continue
        scores[i] = max(DataStructs.BulkTanimotoSimilarity(r["_fp"], pool))
    return {"scores": scores,
            "labels": np.array([1 if r["label"] == "known_active" else 0 for r in rows]),
            "n_active_scaffolds": len({a["_scaf"] for a in acts})}


def duplicates(rows: list[dict]) -> dict:
    seen: dict[str, list[int]] = {}
    for i, r in enumerate(rows):
        seen.setdefault(r["_canon"], []).append(i)
    dup = {k: v for k, v in seen.items() if len(v) > 1}
    cross = {k: v for k, v in dup.items()
             if len({rows[i]["label"] for i in v}) > 1}
    return {
        "n_unique": len(seen),
        "n_duplicate_structures": len(dup),
        "n_rows_involved": sum(len(v) for v in dup.values()),
        "n_cross_label": len(cross),
        "cross_label_examples": [
            {"smiles": k, "labels": [rows[i]["label"] for i in v],
             "ids": [rows[i].get("chembl_id") for i in v]}
            for k, v in list(cross.items())[:5]
        ],
    }


def analog_leak(rows: list[dict]) -> dict:
    acts = [r for r in rows if r["label"] == "known_active"]
    decs = [r for r in rows if r["label"] != "known_active"]
    if not acts or not decs:
        return {"n_pairs_over": 0, "note": "한쪽이 비었다"}
    dfps = [r["_fp"] for r in decs]
    n_over = n_over_cap = 0
    worst = 0.0
    worst_pair = None
    for a in acts:
        sims = DataStructs.BulkTanimotoSimilarity(a["_fp"], dfps)
        arr = np.asarray(sims)
        n_over += int((arr >= ANALOG_T).sum())
        n_over_cap += int((arr > DECOY_CAP).sum())
        j = int(arr.argmax())
        if arr[j] > worst:
            worst, worst_pair = float(arr[j]), (a.get("chembl_id"), decs[j].get("chembl_id"))
    return {
        "threshold": ANALOG_T,
        "n_active_decoy_pairs_over_threshold": n_over,
        "decoy_cap": DECOY_CAP,
        "n_pairs_over_decoy_cap": n_over_cap,
        "max_active_decoy_tanimoto": round(worst, 4),
        "max_pair": worst_pair,
    }


def active_self_similarity(rows: list[dict]) -> dict:
    acts = [r for r in rows if r["label"] == "known_active"]
    if len(acts) < 2:
        return {"note": "활성이 2개 미만"}
    fps = [r["_fp"] for r in acts]
    best = []
    for i, r in enumerate(acts):
        sims = DataStructs.BulkTanimotoSimilarity(r["_fp"], fps)
        sims[i] = -1.0
        best.append(max(sims))
    b = np.array(best)
    return {
        "n_actives": len(acts),
        "median_nearest_neighbour": round(float(np.median(b)), 4),
        "mean_nearest_neighbour": round(float(b.mean()), 4),
        "frac_with_analog_over_0.6": round(float((b >= 0.6).mean()), 4),
        "frac_with_analog_over_0.8": round(float((b >= 0.8).mean()), 4),
    }


def main() -> int:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from core.validation.gate import (THRESHOLDS, bedroc, enrichment_factor,
                                      normalized_ef, roc_auc)

    ap = argparse.ArgumentParser()
    ap.add_argument("--library", required=True, help="library.jsonl 또는 docked_all.jsonl")
    ap.add_argument("--out", default=None)
    ap.add_argument("--label", default=None, help="보고용 이름")
    args = ap.parse_args()

    rows = load(Path(args.library))
    name = args.label or Path(args.library).parent.name
    n_act = sum(1 for r in rows if r["label"] == "known_active")
    print(f"[*] {name} — {len(rows):,}개 (활성 {n_act})\n")

    out: dict = {"library": args.library, "label": name,
                 "n": len(rows), "n_active": n_act}

    # 1. 최대유사도 베이스라인 — 핵심 검사
    mb = max_sim_baseline(rows)
    jitter = np.random.RandomState(0).rand(len(rows)) * 1e-9
    order = np.argsort(-(mb["scores"] + jitter))
    lab = [int(mb["labels"][i]) for i in order]
    m = {
        "ef1": enrichment_factor(lab, 0.01),
        "ef5": enrichment_factor(lab, 0.05),
        "auc": roc_auc(lab),
        "bedroc": bedroc(lab),
    }
    m["nef1"] = normalized_ef(m["ef1"], len(rows), n_act, 0.01)
    hits = [k for k, t in THRESHOLDS.items() if (m.get(k) or 0) >= t]
    m["gate_thresholds"] = THRESHOLDS
    m["gate_passed"] = hits
    m["would_pass_gate"] = len(hits) >= 2
    out["max_similarity_baseline"] = m

    print("[1] 최대유사도 베이스라인 — 단백질도 구조도 도킹도 없이, 활성끼리의 닮음만")
    print(f"      EF1%   {m['ef1']}")
    print(f"      nEF1%  {m['nef1']}   (1.0 = 천장)")
    print(f"      AUC    {m['auc']}")
    print(f"      BEDROC {m['bedroc']}")
    if hits:
        print(f"    → 게이트 문턱 통과: {', '.join(hits)} ({len(hits)}/3)")
        if len(hits) >= 2:
            print("    → 이 벤치마크는 기억을 보상한다. 도킹 결과를 신뢰할 수 없다.")
    else:
        print("    → 게이트 문턱을 하나도 넘지 못한다 — 기억만으로는 안 풀린다")

    # 1b. 골격 분할 — 기억이 아니라 일반화를 묻는다
    sb = max_sim_baseline_scaffold_split(rows)
    order2 = np.argsort(-(sb["scores"] + jitter))
    lab2 = [int(sb["labels"][i]) for i in order2]
    m2 = {
        "ef1": enrichment_factor(lab2, 0.01),
        "auc": roc_auc(lab2),
        "bedroc": bedroc(lab2),
        "n_active_scaffolds": sb["n_active_scaffolds"],
    }
    m2["nef1"] = normalized_ef(m2["ef1"], len(rows), n_act, 0.01)
    hits2 = [k for k, t in THRESHOLDS.items() if (m2.get(k) or 0) >= t]
    m2["gate_passed"] = hits2
    m2["would_pass_gate"] = len(hits2) >= 2
    out["max_similarity_scaffold_split"] = m2

    print(f"\n[1b] 같은 골격을 참조에서 뺀 경우 — 활성 골격 {m2['n_active_scaffolds']}종")
    print(f"      EF1% {m2['ef1']}   nEF1% {m2['nef1']}   AUC {m2['auc']}   BEDROC {m2['bedroc']}")
    if len(hits2) >= 2:
        print("    → 골격을 빼도 통과한다. 활성군 자체가 벤치마크로 못 쓴다.")
    elif hits:
        print("    → 골격을 빼니 무너진다. [1] 의 통과는 같은 시리즈 유사체를")
        print("       외운 것이었다 — 골격 분할 평가로 고칠 수 있다.")
    else:
        print("    → 양쪽 다 통과 못 한다")

    # 2. 중복
    d = duplicates(rows)
    out["duplicates"] = d
    print(f"\n[2] 동일 분자 중복 — 고유 구조 {d['n_unique']:,} / 행 {len(rows):,}")
    print(f"      중복 구조 {d['n_duplicate_structures']}개 (행 {d['n_rows_involved']}개)")
    print(f"      라벨이 엇갈린 중복 {d['n_cross_label']}개", end="")
    print("  ← 라벨링 오류다" if d["n_cross_label"] else "")

    # 3. 유사쌍 누수
    a = analog_leak(rows)
    out["analog_leakage"] = a
    print(f"\n[3] 활성-decoy 유사쌍")
    print(f"      Tanimoto >= {ANALOG_T} 인 쌍  {a.get('n_active_decoy_pairs_over_threshold')}")
    print(f"      decoy 상한 {DECOY_CAP} 초과 쌍  {a.get('n_pairs_over_decoy_cap')}", end="")
    print("   ← 상한이 지켜지지 않았다" if a.get("n_pairs_over_decoy_cap") else "   (상한 지켜짐)")
    print(f"      최대 유사도 {a.get('max_active_decoy_tanimoto')}  {a.get('max_pair')}")

    # 4. 활성 자기유사도
    s = active_self_similarity(rows)
    out["active_self_similarity"] = s
    print(f"\n[4] 활성 자기유사도 — 활성군이 한 골격에 뭉쳐 있는가")
    print(f"      최근접 이웃 중위값 {s.get('median_nearest_neighbour')}")
    print(f"      유사체(>=0.6) 가진 활성 비율 {s.get('frac_with_analog_over_0.6')}")
    print(f"      유사체(>=0.8) 가진 활성 비율 {s.get('frac_with_analog_over_0.8')}")

    if args.out:
        Path(args.out).write_text(json.dumps(out, ensure_ascii=False, indent=2))
        print(f"\n    → {args.out}")

    # 핵심 검사가 통과하면 비정상 종료 — CI 에서 막히도록
    return 1 if m["would_pass_gate"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
