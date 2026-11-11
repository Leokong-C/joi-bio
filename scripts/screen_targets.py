#!/usr/bin/env python3
"""
새 파이프라인에 쓸 **타깃을 고른다.** 화합물이 아니라 타깃이 산출물이다.

## 왜 이 스크립트가 필요한가

2026-11-07 에 EPAS1 이 벤치마크로 실격했다. 이유는 decoy 규칙이 아니라
**활성군 자체**였다.

    EPAS1 활성 60개 — 상호 유사도 중위 0.795, 87% 가 0.6 이상 유사체 보유
    → 활성끼리 닮았다는 것만 세는 기억 베이스라인이 AUC 0.9637 로
      도킹(0.7717)을 이긴다. 도킹이 작동하는지 물을 수 없는 벤치마크다.

그래서 타깃을 옮겨야 하고, **옮길 타깃은 같은 잣대로 미리 재야 한다.**
안 재고 옮기면 같은 함정에 다시 들어간다.

## 합격 조건 — EPAS1 의 탈락 사유를 뒤집은 것

| | 조건 | EPAS1 |
| --- | --- | --- |
| 1 | 기억 베이스라인이 게이트를 통과하지 **못할 것** | PASS 3/3 → 실격 |
| 2 | 활성군 중위 NN Tanimoto 가 낮을 것 | 0.795 → 실격 |
| 3 | 비활성이 **실측**일 것 (추정 decoy 금지) | 비활성 0개 → 실격 |
| 4 | 활성 골격이 여러 종일 것 | 한 시리즈 |

3번이 TDC Butkiewicz 셋을 쓰는 이유다. 거기 비활성은 추정이 아니라
실제로 시험해서 음성으로 나온 화합물이다.

## 두 가지를 함께 낸다

- **맞춘 부분표본** (2,157 / 활성 60) — EPAS1 과 n·활성비율이 같아야
  EF 천장이 같고, 그래야 숫자를 나란히 놓을 수 있다
- **전체 활성군의 응집도** — 부분표본 60개가 아니라 활성 전수로 잰다.
  타깃의 성질은 전수에서 나온다

음성 결과도 적는다. 떨어진 타깃과 그 이유가 다음 사람에게 필요하다.
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.leakage_audit import (  # noqa: E402
    FPGEN, active_self_similarity, max_sim_baseline,
    max_sim_baseline_scaffold_split, scaffold_of)

EPAS1_N = 2157
EPAS1_ACTIVES = 60


def build_rows(smiles_labels: list[tuple[str, int]]) -> list[dict]:
    """SMILES 를 감사 함수가 먹는 형태로. 실패한 분자는 조용히 버리지 않고 센다."""
    from rdkit import Chem

    rows, bad = [], 0
    for smi, y in smiles_labels:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            bad += 1
            continue
        rows.append({
            "smiles": smi,
            "label": "known_active" if y == 1 else "decoy",
            "_mol": m,
            "_fp": FPGEN.GetFingerprint(m),
        })
    if bad:
        print(f"    [!] 파싱 실패 {bad:,}개는 제외했다")
    return rows


def gate_metrics(rows: list[dict]) -> dict:
    import numpy as np

    from core.validation.gate import (THRESHOLDS, bedroc, enrichment_factor,
                                      normalized_ef, roc_auc)

    n_act = sum(1 for r in rows if r["label"] == "known_active")
    mb = max_sim_baseline(rows)
    jitter = np.random.RandomState(0).rand(len(rows)) * 1e-9
    order = np.argsort(-(mb["scores"] + jitter))
    lab = [int(mb["labels"][i]) for i in order]
    m = {
        "ef1": enrichment_factor(lab, 0.01),
        "auc": roc_auc(lab),
        "bedroc": bedroc(lab),
    }
    m["nef1"] = normalized_ef(m["ef1"], len(rows), n_act, 0.01)
    hits = [k for k, t in THRESHOLDS.items() if (m.get(k) or 0) >= t]
    m["gate_passed"] = hits
    m["would_pass_gate"] = len(hits) >= 2
    return m


def screen_one(name: str, data_dir: str, seed: int) -> dict:
    import warnings
    warnings.filterwarnings("ignore")
    from tdc.single_pred import HTS

    print(f"\n{'=' * 72}\n[*] {name}\n{'=' * 72}")
    df = HTS(name=name, path=data_dir).get_data()
    smi_col = "Drug" if "Drug" in df.columns else df.columns[1]

    pairs = [(str(r[smi_col]), int(r["Y"])) for _, r in df.iterrows()]
    n_total = len(pairs)
    n_act_total = sum(1 for _, y in pairs if y == 1)
    print(f"    전수 {n_total:,}개, 활성 {n_act_total:,} "
          f"({100 * n_act_total / n_total:.3f}%), 비활성은 실측")

    res: dict = {
        "target": name,
        "n_total": n_total,
        "n_active_total": n_act_total,
        "active_fraction_total": round(n_act_total / n_total, 6),
    }

    # --- 활성 전수의 응집도. 타깃의 성질은 여기서 나온다 ---
    act_rows = build_rows([(s, 1) for s, y in pairs if y == 1])
    sim = active_self_similarity(act_rows)
    res["active_self_similarity_all"] = sim
    res["n_active_scaffolds_all"] = len({scaffold_of(r["_mol"]) for r in act_rows})
    print(f"    활성 전수 응집도 — 중위 NN {sim.get('median_nearest_neighbour')}, "
          f"0.6 이상 {sim.get('frac_with_analog_over_0.6')}, "
          f"골격 {res['n_active_scaffolds_all']}종")

    # --- EPAS1 과 나란히 놓을 수 있는 부분표본 ---
    rng = random.Random(seed)
    acts = [p for p in pairs if p[1] == 1]
    decs = [p for p in pairs if p[1] == 0]
    na = min(EPAS1_ACTIVES, len(acts))
    nd = min(EPAS1_N - na, len(decs))
    sub = rng.sample(acts, na) + rng.sample(decs, nd)
    res["matched"] = {"n": na + nd, "n_active": na,
                      "short_of_epas1": na < EPAS1_ACTIVES or nd < EPAS1_N - na}

    rows = build_rows(sub)
    res["matched"]["memorisation"] = gate_metrics(rows)
    sb = max_sim_baseline_scaffold_split(rows)
    res["matched"]["n_active_scaffolds"] = sb["n_active_scaffolds"]

    m = res["matched"]["memorisation"]
    verdict = "실격 — 기억이 게이트 통과" if m["would_pass_gate"] else "합격 후보"
    res["verdict"] = verdict
    print(f"    기억 베이스라인 — nEF1% {m['nef1']}  AUC {m['auc']}  "
          f"BEDROC {m['bedroc']}  → {verdict}")
    return res


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", default="ALL_BUTKIEWICZ",
                    help="쉼표 구분 TDC HTS 이름, 또는 ALL_BUTKIEWICZ")
    ap.add_argument("--data-dir", default="tdc_raw")
    ap.add_argument("--out", default="screen_targets.json")
    ap.add_argument("--seed", type=int, default=1023)
    args = ap.parse_args()

    if args.names.strip() == "ALL_BUTKIEWICZ":
        from tdc.utils import retrieve_dataset_names
        names = sorted(n for n in retrieve_dataset_names("HTS")
                       if "butkiewicz" in n.lower())
        print(f"[*] Butkiewicz HTS {len(names)}종을 찾았다")
    else:
        names = [s.strip() for s in args.names.split(",") if s.strip()]

    results, failed = [], []
    for n in names:
        try:
            results.append(screen_one(n, args.data_dir, args.seed))
        except Exception as e:          # 하나가 죽어도 나머지는 잰다
            print(f"    [!] {n} 실패: {type(e).__name__}: {e}")
            failed.append({"target": n, "error": f"{type(e).__name__}: {e}"})

    # EPAS1 은 실측 기준선으로 항상 같이 싣는다 — 비교 대상이 없으면 숫자가 안 읽힌다
    ref = {"target": "EPAS1 (현행, Tanimoto 상한 0.35 decoy)", "reference": True,
           "n_total": 2160, "n_active_total": 60,
           "active_self_similarity_all": {"median_nearest_neighbour": 0.795,
                                          "frac_with_analog_over_0.6": 0.87},
           "matched": {"memorisation": {"nef1": 1.000, "auc": 0.9637,
                                        "would_pass_gate": True}},
           "verdict": "실격 — 기억이 게이트 통과"}

    out = {"reference": ref, "screened": results, "failed": failed}
    Path(args.out).write_text(json.dumps(out, ensure_ascii=False, indent=2))

    print(f"\n\n{'=' * 100}\n순위 — 기억 베이스라인 AUC 가 낮을수록 좋은 타깃이다\n{'=' * 100}")
    print(f"{'타깃':42s} {'활성':>7s} {'골격':>6s} {'중위NN':>7s} "
          f"{'≥0.6':>6s} {'기억AUC':>8s} {'판정'}")
    print("-" * 100)
    rank = sorted(results, key=lambda r: r["matched"]["memorisation"]["auc"] or 0)
    for r in [ref] + rank:
        s = r["active_self_similarity_all"]
        mm = r["matched"]["memorisation"]
        print(f"{r['target'][:42]:42s} {r['n_active_total']:7,d} "
              f"{r.get('n_active_scaffolds_all', '-'):>6} "
              f"{s.get('median_nearest_neighbour', 0):7.3f} "
              f"{s.get('frac_with_analog_over_0.6', 0):6.2f} "
              f"{(mm['auc'] or 0):8.4f}  {r['verdict']}")
    for f in failed:
        print(f"{f['target'][:42]:42s} {'—':>7} {'—':>6} {'—':>7} {'—':>6} "
              f"{'—':>8}  측정 실패: {f['error'][:40]}")
    print(f"\n합격 후보 {sum(1 for r in results if not r['matched']['memorisation']['would_pass_gate'])}종 "
          f"/ 측정 {len(results)}종 / 실패 {len(failed)}종")
    print(f"→ {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
