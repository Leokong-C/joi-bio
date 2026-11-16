#!/usr/bin/env python3
"""
임의의 타깃에 **실측 비활성**이 PubChem BioAssay 에 있는가. 선정 기준 C3.

`probe_epas1_bioassay.py` 를 타깃 하나에 묶이지 않게 일반화한 것이다.
EPAS1 전용본은 tdc-audit.yml 이 참조하므로 그대로 둔다.

## 왜 이걸 먼저 묻는가

C3(비활성이 실측일 것)가 타깃 풀을 정한다. 추정 decoy 를 쓰면 Tanimoto 상한이
다시 들어와 지금 재려는 것을 가린다. 그래서 **어떤 타깃을 검토하든 이 질문이
먼저**다 — 여기서 막히면 도킹이고 구조고 의미가 없다.

EPAS1 이 그렇게 막혔다. 어세이 155개가 전부 소규모 문헌 IC50 이었고
**실측 비활성이 0개**였다.

없으면 그것도 결과다. 숨기지 말고 기록한다.

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
MIN_TESTED = 2000  # 우리 평가셋이 2,157개다


def get(url: str, retries: int = 3):
    for i in range(retries):
        try:
            r = requests.get(url, timeout=90, headers=UA)
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
        "name": str(s.get("Name", ""))[:140],
        "n_tested": s.get("CIDCountAll", 0),
        "n_active": s.get("CIDCountActive", 0),
        "n_inactive": s.get("CIDCountInactive", 0),
    }


def probe(label: str, symbols: list[str], max_aids: int) -> dict:
    print(f"\n{'=' * 72}\n[*] {label} — 별칭 {', '.join(symbols)}\n{'=' * 72}")
    found: dict[int, str] = {}
    per_symbol = {}
    for sym in symbols:
        a = aids_for(sym)
        per_symbol[sym] = len(a)
        print(f"    {sym:12s} AID {len(a)}개")
        for x in a:
            found.setdefault(x, sym)

    out: dict = {"label": label, "symbols": per_symbol, "n_aids": len(found),
                 "min_tested_required": MIN_TESTED}
    if not found:
        out["assays"], out["n_usable"] = [], 0
        out["verdict"] = "어세이를 찾지 못했다 — C3 실패"
        print(f"\n    [=] {out['verdict']}")
        return out

    print(f"\n    고유 AID {len(found)}개 — 규모 확인 (최대 {max_aids}개)")
    rows = []
    for aid in list(found)[:max_aids]:
        s = summarize(aid)
        if s is None:
            continue
        s["symbol"] = found[aid]
        rows.append(s)
    rows.sort(key=lambda r: -r["n_inactive"])

    for s in rows[:12]:
        print(f"      AID {s['aid']:>9}  시험 {s['n_tested']:>8,}  "
              f"활성 {s['n_active']:>6,}  비활성 {s['n_inactive']:>8,}  "
              f"{s['name'][:46]}")
    if len(rows) > 12:
        print(f"      … 외 {len(rows) - 12}개")

    usable = [r for r in rows if r["n_tested"] >= MIN_TESTED and r["n_inactive"] > 0]
    out["assays"], out["n_usable"] = rows, len(usable)
    out["total_inactive_seen"] = sum(r["n_inactive"] for r in rows)
    out["max_inactive_single_assay"] = rows[0]["n_inactive"] if rows else 0

    if usable:
        best = max(usable, key=lambda r: r["n_inactive"])
        out["verdict"] = (f"C3 통과 — AID {best['aid']} 에 실측 비활성 "
                          f"{best['n_inactive']:,}개")
        out["c3_pass"] = True
    else:
        out["verdict"] = (f"C3 실패 — 시험 {MIN_TESTED:,}개 이상이면서 비활성이 "
                          f"있는 어세이가 없다. 최대 비활성 "
                          f"{out['max_inactive_single_assay']:,}개")
        out["c3_pass"] = False
    print(f"\n    [=] {out['verdict']}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", required=True,
                    help="'라벨=별칭1|별칭2' 를 쉼표로. 예: PH20=SPAM1|PH20,HYAL1=HYAL1")
    ap.add_argument("--max-aids", type=int, default=120)
    ap.add_argument("--out", default="bioassay.json")
    args = ap.parse_args()

    res = []
    for spec in args.targets.split(","):
        spec = spec.strip()
        if not spec:
            continue
        label, _, rest = spec.partition("=")
        syms = [s.strip() for s in (rest or label).split("|") if s.strip()]
        try:
            res.append(probe(label.strip(), syms, args.max_aids))
        except Exception as e:
            print(f"    [!] {label} 실패: {type(e).__name__}: {e}")
            res.append({"label": label.strip(), "error": f"{type(e).__name__}: {e}"})

    Path(args.out).write_text(json.dumps(res, ensure_ascii=False, indent=2))

    print(f"\n\n{'=' * 86}\nC3 — 실측 비활성이 있는가\n{'=' * 86}")
    print(f"{'타깃':14s} {'AID':>7s} {'최대 비활성/어세이':>18s} "
          f"{'쓸만한 어세이':>13s}  판정")
    print("-" * 86)
    for r in res:
        if "error" in r:
            print(f"{r['label'][:14]:14s} {'—':>7s} {'—':>18s} {'—':>13s}  "
                  f"측정 실패")
            continue
        print(f"{r['label'][:14]:14s} {r['n_aids']:7,d} "
              f"{r.get('max_inactive_single_assay', 0):18,d} "
              f"{r['n_usable']:13,d}  {r['verdict'][:40]}")
    print(f"→ {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
