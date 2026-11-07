#!/usr/bin/env python3
"""
TDC 의 HTS 데이터셋을 우리 라이브러리 형식으로 바꾼다.

## 왜 필요한가

2026-10-29 에 다섯 번째 교란을 찾았다. 활성과의 최대 ECFP4 Tanimoto 만으로
순위를 매기는 **기억 베이스라인이 AUC 0.9637 로 도킹(0.7717)을 이긴다.**

원인 진단은 이랬다: `decoys.select` 가 decoy-활성 Tanimoto 를 0.35 로 자르고,
활성의 87% 는 0.6 이상의 유사체를 갖는다. **상한 아래에 decoy 를 묶어두면
그 위는 전부 활성**이므로 문턱 하나로 갈린다(0.478 에서 99.77% 정확도).

이건 **진단**이지 증명이 아니다. 증명하려면 **상한이 없는 벤치마크**에서
기억 베이스라인이 무너지는지 봐야 한다.

TDC 의 Butkiewicz HTS 셋이 그것이다. decoy 가 추정이 아니라 **실제로 시험해서
비활성으로 나온 화합물**이라, 활성과 닮았다는 이유로 걸러낼 필요가 없었다.
여기서도 기억이 이기면 내 진단이 틀린 것이고, 무너지면 맞은 것이다.

## 비교를 성립시키는 것

TDC HTS 셋은 수십만 개에 활성 비율이 0.1% 수준이다. 우리 EPAS1 은 2,157개에
2.78% 다. 규모와 활성 비율이 다르면 EF 천장이 달라 비교가 안 된다
(그래서 nEF1% 를 쓰지만, 기억 베이스라인의 AUC 자체도 표본 크기에 흔들린다).

그래서 두 가지를 다 낸다:
  1. 전수 — 그 벤치마크의 실제 성질
  2. **맞춘 부분표본** — n 과 활성 비율을 우리 EPAS1 과 같게 (2,157 / 60)

2번이 우리 숫자와 나란히 놓을 수 있는 쪽이다.
"""
from __future__ import annotations

import argparse
import json
import random
from pathlib import Path


def main() -> int:
    import warnings
    warnings.filterwarnings("ignore")
    from tdc.single_pred import HTS

    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True, help="TDC HTS 데이터셋 이름")
    ap.add_argument("--data-dir", default="tdc_raw")
    ap.add_argument("--out-full", required=True)
    ap.add_argument("--out-matched", default=None)
    ap.add_argument("--match-n", type=int, default=2157, help="EPAS1 평가셋 크기")
    ap.add_argument("--match-actives", type=int, default=60, help="EPAS1 활성 수")
    ap.add_argument("--seed", type=int, default=1023)
    args = ap.parse_args()

    d = HTS(name=args.name, path=args.data_dir)
    df = d.get_data()
    # TDC 의 열 이름은 Drug_ID / Drug / Y
    smi_col = "Drug" if "Drug" in df.columns else df.columns[1]
    id_col = "Drug_ID" if "Drug_ID" in df.columns else df.columns[0]

    rows = []
    for _, r in df.iterrows():
        y = int(r["Y"])
        rows.append({
            "smiles": str(r[smi_col]),
            "chembl_id": f"TDC{r[id_col]}",   # 파이프라인이 이 키를 쓴다
            "pref_name": "",
            # 비활성은 추정이 아니라 실측이다 — 그래도 라벨 이름은 파이프라인과 맞춘다
            "label": "known_active" if y == 1 else "decoy",
            "source": f"tdc_hts:{args.name}",
            "measured_inactive": y == 0,
            "protected": False,
        })

    n_act = sum(1 for r in rows if r["label"] == "known_active")
    print(f"[*] {args.name} — {len(rows):,}개, 활성 {n_act:,} "
          f"({100 * n_act / len(rows):.3f}%)")
    print(f"    decoy 는 추정이 아니라 실측 비활성이다")

    Path(args.out_full).write_text(
        "\n".join(json.dumps(r, ensure_ascii=False) for r in rows) + "\n")
    print(f"    전수 → {args.out_full}")

    if args.out_matched:
        rng = random.Random(args.seed)
        acts = [r for r in rows if r["label"] == "known_active"]
        decs = [r for r in rows if r["label"] != "known_active"]
        na = min(args.match_actives, len(acts))
        nd = args.match_n - na
        if nd > len(decs):
            print(f"    [!] 비활성이 {len(decs):,}개뿐이라 {nd:,}개를 못 뽑는다")
            nd = len(decs)
        sub = rng.sample(acts, na) + rng.sample(decs, nd)
        rng.shuffle(sub)
        Path(args.out_matched).write_text(
            "\n".join(json.dumps(r, ensure_ascii=False) for r in sub) + "\n")
        print(f"    맞춘 부분표본 {len(sub):,}개 (활성 {na}, "
              f"{100 * na / len(sub):.2f}%) → {args.out_matched}")
        print(f"    ↑ 우리 EPAS1 (2,157 / 60 / 2.78%) 와 나란히 놓을 수 있는 쪽")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
