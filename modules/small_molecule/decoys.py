#!/usr/bin/env python3
"""
물성 매칭 decoy 선택 — 벤치마크를 어렵게 만든다.

## 왜 필요한가

2026-10-24, 단백질도 구조도 도킹도 쓰지 않고 **불소 원자만 센** 순위가
EPAS1 라이브러리에서 AUC 0.9114, BEDROC 0.687 을 냈다. 도킹 치료군이
0.7717 / 0.6697 이다. 분자에서 원자 하나 세는 것이 도킹을 이긴다.

원인은 decoy 다. EPAS1 활성은 PT 계열 — 불소가 많이 붙은 설폰·다이아릴
에테르인데, 기존 decoy 는 "분자량 220~550 인 무작위 ChEMBL 화합물" 이라
불소가 거의 없다. 활성 중앙값 3개, decoy 0개. 불소 3개 이상이 활성의 80%,
decoy 의 4%. 물성 9개만으로도 AUC 0.905 가 나온다.

즉 EF1%·AUC·BEDROC 의 상당 부분은 "도킹이 활성을 잘 찾아서"가 아니라
"활성이 decoy 와 애초에 다른 분자라서" 나온 값이다.

## 어떻게 고치는가

DUD-E 방식을 따른다 — 활성 하나하나에 대해 **물성은 닮고 토폴로지는 다른**
화합물을 골라 decoy 로 쓴다.

  물성 매칭   분자량, cLogP, 회전결합, HBD, HBA, 형식전하, 할로겐 수
  토폴로지 분리  모든 활성에 대한 ECFP4 Tanimoto 최대값이 문턱 미만

할로겐 수를 매칭 항목에 넣은 것은 DUD-E 표준에는 없는 추가다. 우리 실패가
구체적으로 불소였기 때문이고, 분자량과 logP 만 맞춰서는 비불소 화합물로도
맞출 수 있어 지름길이 남는다.

## 성공을 어떻게 아는가

추측하지 않는다. `scripts/trivial_baselines.py` 를 새 라이브러리에 돌려
**단일 기술자 하나가 내는 최대 AUC** 를 본다. 이것이 수용 기준이다.

  AUC <= 0.60   쓸만하다
  0.60 ~ 0.70   경계 — 어느 기술자가 새는지 보고 판단
  > 0.70        여전히 지름길이 있다. 더 조여야 한다

기존 EPAS1 라이브러리는 이 기준으로 **0.9114 (불소) — 탈락**이다.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")

# 매칭할 물성. 이름 → (계산 함수, 허용 오차의 기준이 되는 눈금)
#
# 눈금은 "이만큼 차이나면 1 단위 떨어진 것으로 센다" 는 뜻이다. DUD-E 가
# 쓰는 구간 폭을 참고해 잡았고, 할로겐은 개수 자체가 작아 1개를 1 단위로 둔다.
PROPS = {
    "mw": (Descriptors.MolWt, 25.0),
    "clogp": (Descriptors.MolLogP, 0.5),
    "rotb": (lambda m: float(Descriptors.NumRotatableBonds(m)), 1.0),
    "hbd": (lambda m: float(Descriptors.NumHDonors(m)), 1.0),
    "hba": (lambda m: float(Descriptors.NumHAcceptors(m)), 1.0),
    "charge": (lambda m: float(Chem.GetFormalCharge(m)), 1.0),
    "halogen": (
        lambda m: float(len(m.GetSubstructMatches(Chem.MolFromSmarts("[F,Cl,Br,I]")))),
        1.0,
    ),
}

# 작용기 매칭 — 기본 물성 7개로 부족할 때 켠다.
#
# EPAS1 에서 확인한 것 (2026-10-24): 할로겐 총수를 맞춰도 불소 AUC 가 0.816 으로
# 남는다. 염소 3개로 불소 3개를 "할로겐 3개" 로 맞춰버리기 때문이다. 설폰과
# 나이트릴은 애초에 매칭 대상이 아니어서 그대로 샌다 (AUC 0.798, 0.722).
# PT 계열처럼 작용기 서명이 뚜렷한 활성군에는 이쪽까지 맞춰야 한다.
EXTRA_PROPS = {
    "fluorine": ("[F]", 1.0),
    "sulfone": ("[SX4](=O)(=O)", 1.0),
    "nitrile": ("[C]#[N]", 1.0),
    "aryl_ether": ("c-O-c", 1.0),
}

_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def _active_props(extra: bool) -> dict:
    """매칭에 쓸 항목. extra=True 면 작용기 네 가지를 더한다."""
    p = dict(PROPS)
    if extra:
        for name, (smarts, scale) in EXTRA_PROPS.items():
            patt = Chem.MolFromSmarts(smarts)
            p[name] = ((lambda pt: lambda m: float(len(m.GetSubstructMatches(pt))))(patt), scale)
    return p


def props(mol: Chem.Mol, table: dict | None = None) -> dict[str, float]:
    return {k: float(fn(mol)) for k, (fn, _) in (table or PROPS).items()}


def distance(a: dict[str, float], b: dict[str, float], table: dict | None = None) -> float:
    """눈금으로 나눈 차이의 제곱합. 작을수록 닮았다."""
    t = table or PROPS
    return math.sqrt(sum(((a[k] - b[k]) / scale) ** 2 for k, (_, scale) in t.items()))


def select(
    actives: list[dict],
    pool: list[dict],
    per_active: int,
    max_tanimoto: float = 0.35,
    verbose: bool = True,
    match_groups: bool = True,
) -> tuple[list[dict], dict]:
    """활성마다 물성이 가장 닮고 토폴로지는 충분히 다른 화합물을 고른다."""
    table = _active_props(match_groups)
    amols, aprops, afps = [], [], []
    for a in actives:
        m = Chem.MolFromSmiles(a["smiles"])
        if m is None:
            continue
        amols.append(a)
        aprops.append(props(m, table))
        afps.append(_GEN.GetFingerprint(m))

    pmols = []
    for p in pool:
        m = Chem.MolFromSmiles(p["smiles"])
        if m is not None:
            pmols.append((p, props(m, table), _GEN.GetFingerprint(m)))

    if verbose:
        print(f"[*] 활성 {len(amols)}개, 후보 풀 {len(pmols)}개")
        print(f"    활성당 {per_active}개, Tanimoto 상한 {max_tanimoto}")

    # 모든 활성에 대한 최대 유사도를 미리 계산 — 토폴로지가 겹치는 후보를 먼저 뺀다
    eligible = []
    n_too_similar = 0
    for rec, pr, fp in pmols:
        sims = DataStructs.BulkTanimotoSimilarity(fp, afps)
        mx = max(sims) if sims else 0.0
        if mx >= max_tanimoto:
            n_too_similar += 1
            continue
        eligible.append((rec, pr, mx))
    if verbose:
        print(f"    토폴로지 제외 {n_too_similar}개 → 남은 후보 {len(eligible)}개")

    used: set[str] = set()
    picked: list[dict] = []
    dists: list[float] = []
    for a, ap in zip(amols, aprops):
        cand = sorted(
            (
                (distance(ap, pr, table), rec, mx)
                for rec, pr, mx in eligible
                if rec["chembl_id"] not in used
            ),
            key=lambda t: t[0],
        )[:per_active]
        for d, rec, mx in cand:
            used.add(rec["chembl_id"])
            r = dict(rec)
            r["label"] = "decoy"
            r["source"] = "property_matched"
            r["matched_to"] = a["chembl_id"]
            r["prop_distance"] = round(d, 3)
            r["max_tanimoto_to_actives"] = round(mx, 3)
            picked.append(r)
            dists.append(d)

    stats = {
        "n_actives": len(amols),
        "n_pool": len(pmols),
        "n_excluded_similar": n_too_similar,
        "n_eligible": len(eligible),
        "per_active_requested": per_active,
        "n_selected": len(picked),
        "max_tanimoto_threshold": max_tanimoto,
        "prop_distance": {
            "median": round(sorted(dists)[len(dists) // 2], 3) if dists else None,
            "max": round(max(dists), 3) if dists else None,
        },
        "properties_matched": list(table),
        "match_groups": match_groups,
    }
    if verbose:
        print(f"    선택 {len(picked)}개, 물성 거리 중앙값 {stats['prop_distance']['median']}")
    return picked, stats


def compare(actives: list[dict], decoys: list[dict]) -> dict:
    """활성과 decoy 의 물성 중앙값을 나란히 놓는다. 매칭이 됐는지 눈으로 본다."""
    def med(rows, key):
        vals = []
        tbl = _active_props(True)
        for r in rows:
            m = Chem.MolFromSmiles(r["smiles"])
            if m is not None:
                vals.append(tbl[key][0](m))
        vals.sort()
        return round(vals[len(vals) // 2], 2) if vals else None

    out = {}
    for k in _active_props(True):
        a, d = med(actives, k), med(decoys, k)
        out[k] = {"active": a, "decoy": d,
                  "diff": round(a - d, 2) if a is not None and d is not None else None}
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--library", required=True, help="기존 library.jsonl")
    ap.add_argument("--pool", default=None,
                    help="후보 풀 jsonl. 없으면 library 의 decoy 를 풀로 쓴다")
    ap.add_argument("--per-active", type=int, default=35)
    ap.add_argument("--max-tanimoto", type=float, default=0.35)
    ap.add_argument("--no-group-matching", action="store_true",
                    help="작용기(불소·설폰·나이트릴·아릴에테르) 매칭을 끈다")
    ap.add_argument("--out", required=True)
    ap.add_argument("--stats-out", default=None)
    args = ap.parse_args()

    lib = [json.loads(l) for l in Path(args.library).read_text().splitlines() if l.strip()]
    actives = [r for r in lib if r.get("label") == "known_active"]
    pool = (
        [json.loads(l) for l in Path(args.pool).read_text().splitlines() if l.strip()]
        if args.pool
        else [r for r in lib if r.get("label") != "known_active"]
    )

    picked, stats = select(actives, pool, args.per_active, args.max_tanimoto,
                           match_groups=not args.no_group_matching)
    stats["property_medians"] = compare(actives, picked)

    print("\n    물성 비교 (중앙값)")
    print(f"    {'항목':10s} {'활성':>9s} {'decoy':>9s} {'차이':>9s}")
    for k, v in stats["property_medians"].items():
        print(f"    {k:10s} {str(v['active']):>9s} {str(v['decoy']):>9s} {str(v['diff']):>9s}")

    with open(args.out, "w") as f:
        for r in actives + picked:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(f"\n    → {args.out}  (활성 {len(actives)} + decoy {len(picked)})")

    if args.stats_out:
        Path(args.stats_out).write_text(json.dumps(stats, ensure_ascii=False, indent=2))
        print(f"    → {args.stats_out}")

    print("\n[!] 다음: scripts/trivial_baselines.py 로 수용 기준을 확인할 것.")
    print("    단일 기술자 최대 AUC 가 0.60 이하여야 쓸만하다.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
