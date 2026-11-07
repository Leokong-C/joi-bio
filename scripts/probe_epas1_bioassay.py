#!/usr/bin/env python3
"""
EPAS1 에 쓸 수 있는 실측 비활성이 PubChem BioAssay 에 있는가.

## 왜 이걸 묻는가

다섯 번째 교란(기억이 도킹을 이긴다)의 원인은 우리 decoy 가 **추정** 비활성
이라는 것이다. 실제로 비활성인지 모르니, 활성과 너무 닮은 것은 사실 활성일까
봐 Tanimoto 상한 0.35 로 잘라낸다. 그 상한이 간격을 만들고, 그 간격 때문에
"활성과 닮았나" 가 완벽한 분류기가 된다.

**실제로 시험해서 비활성으로 나온 화합물**을 쓰면 상한을 걸 이유가 없어진다.
그러면 이 교란이 사라진다.

문제는 EPAS1 에 그런 데이터가 있느냐다. TDC 의 HTS 목록에는 EPAS1 이 없다.
PubChem BioAssay 에 직접 물어본다.

없으면 그것도 결과다 — 이 경로는 EPAS1 에서 막히고, 타깃을 옮기거나
다른 방법을 찾아야 한다는 뜻이다.

Claude 컨테이너는 PubChem 이 프록시에서 막혀 있어(403) Actions 에서만 돈다.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import requests

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
UA = {"User-Agent": "joi-bio/1.0 (benchmark construction)"}
PAUSE = 0.25

# EPAS1 = HIF-2alpha. 별칭이 여럿이라 전부 두드린다.
SYMBOLS = ["EPAS1", "HIF2A", "MOP2", "PASD2", "HLF", "ECYT4"]

# 쓸만하려면 이 정도는 돼야 한다 — 우리 평가셋이 2,157개다
MIN_TESTED = 2000


def get(url: str, params: dict | None = None, retries: int = 3):
    for i in range(retries):
        try:
            r = requests.get(url, params=params, timeout=90, headers=UA)
            if r.status_code == 200:
                return r
            if r.status_code in (503, 504):
                time.sleep(2 ** i)
                continue
            return None
        except Exception as e:
            print(f"    [!] {type(e).__name__} {str(e)[:90]}")
            time.sleep(2 ** i)
    return None


def aids_for(symbol: str) -> list[int]:
    r = get(f"{PUG}/assay/target/genesymbol/{symbol}/aids/JSON")
    time.sleep(PAUSE)
    if r is None:
        return []
    try:
        return r.json().get("IdentifierList", {}).get("AID", [])
    except Exception:
        return []


def summarize(aid: int) -> dict | None:
    r = get(f"{PUG}/assay/aid/{aid}/summary/JSON")
    time.sleep(PAUSE)
    if r is None:
        return None
    try:
        s = r.json()["AssaySummaries"]["AssaySummary"][0]
    except Exception:
        return None
    return {
        "aid": aid,
        "name": str(s.get("Name", ""))[:120],
        "n_tested": s.get("CIDCountAll", 0),
        "n_active": s.get("CIDCountActive", 0),
        "n_inactive": s.get("CIDCountInactive", 0),
        "activity_outcome_method": s.get("Method", ""),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="epas1_bioassays.json")
    ap.add_argument("--max-aids", type=int, default=40)
    args = ap.parse_args()

    print("[1] 유전자명·별칭으로 어세이 찾기")
    found: dict[int, str] = {}
    per_symbol = {}
    for sym in SYMBOLS:
        a = aids_for(sym)
        per_symbol[sym] = len(a)
        print(f"    {sym:8s} AID {len(a)}개")
        for x in a:
            found.setdefault(x, sym)

    if not found:
        out = {"symbols": per_symbol, "n_aids": 0, "assays": [],
               "verdict": "EPAS1 어세이를 찾지 못했다 — 이 경로는 EPAS1 에서 막힌다"}
        Path(args.out).write_text(json.dumps(out, ensure_ascii=False, indent=2))
        print(f"\n[=] {out['verdict']}")
        return 0

    print(f"\n[2] 고유 AID {len(found)}개 — 규모 확인 (최대 {args.max_aids}개)")
    rows = []
    for aid in list(found)[: args.max_aids]:
        s = summarize(aid)
        if s is None:
            continue
        s["symbol"] = found[aid]
        rows.append(s)
        print(f"    AID {aid:>9}  시험 {s['n_tested']:>8,}  활성 {s['n_active']:>6,}  "
              f"비활성 {s['n_inactive']:>8,}  {s['name'][:44]}")

    usable = [r for r in rows if r["n_tested"] >= MIN_TESTED and r["n_inactive"] > 0]
    rows.sort(key=lambda r: -r["n_inactive"])

    print(f"\n[3] 판정 — 시험 {MIN_TESTED:,}개 이상이고 비활성이 있는 어세이")
    if usable:
        best = max(usable, key=lambda r: r["n_inactive"])
        verdict = (f"쓸 수 있다 — AID {best['aid']} 에 실측 비활성 "
                   f"{best['n_inactive']:,}개")
        print(f"    {len(usable)}개 있다. 가장 큰 것: AID {best['aid']}, "
              f"비활성 {best['n_inactive']:,}개")
        print(f"    → 이걸로 decoy 를 바꾸면 Tanimoto 상한을 걸 이유가 없어진다")
    else:
        biggest = rows[0] if rows else None
        verdict = ("쓸 만한 규모의 어세이가 없다 — EPAS1 에서는 이 경로가 막힌다. "
                   "타깃을 옮기거나 다른 방법을 찾아야 한다")
        print(f"    없다. 가장 큰 것도 비활성 "
              f"{biggest['n_inactive'] if biggest else 0:,}개")
        print(f"    → {verdict}")

    out = {"symbols": per_symbol, "n_aids": len(found),
           "min_tested_required": MIN_TESTED,
           "assays": rows, "n_usable": len(usable), "verdict": verdict}
    Path(args.out).write_text(json.dumps(out, ensure_ascii=False, indent=2))
    print(f"\n    → {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
