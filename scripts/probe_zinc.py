#!/usr/bin/env python3
"""
ZINC 접근 가능성 탐침 — 무엇이 실제로 되는지 알아본다.

왜 탐침인가
  ChEMBL 무작위 풀 2,097개에 PT 서명(불소>=3 + 설폰 + 나이트릴 + 다이아릴에테르)을
  모두 가진 화합물이 0개였다. 활성은 35%가 가지고 있다. 맞출 상대가 없으면
  물성 매칭 decoy 는 작동할 수 없다 (2026-10-24 측정).

  DUD-E 가 ZINC 를 쓰는 이유가 이것이다. ZINC22 는 구매 가능 화합물을 수십억
  규모로 담고 있고, tranche 가 HBD·logP·MW 로 미리 구간화돼 있어 우리가 필요한
  물성 창을 그대로 집을 수 있다 — 고 3자 문서가 말한다. 확인한 적은 없다.

  Claude 가 도는 컨테이너는 아웃바운드가 막혀 있어 여기선 못 본다.
  Actions 러너에서 돌려 실제로 무엇이 되는지 기록한다.

무엇을 재는가
  1. 각 후보 엔드포인트가 응답하는가 (상태코드, 응답 크기, 지연)
  2. 무작위 표본에 우리가 필요한 화학형이 얼마나 있는가
  3. tranche 디렉터리 구조가 문서대로인가

아무것도 안 되면 그것도 결과다. 그때는 PubChem 이나 Enamine REAL 로 간다.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import requests

TIMEOUT = 45
UA = {"User-Agent": "joi-bio/1.0 (virtual screening benchmark construction)"}


def probe(name: str, url: str, params: dict | None = None, head: bool = False) -> dict:
    t = time.time()
    rec = {"name": name, "url": url, "params": params}
    try:
        fn = requests.head if head else requests.get
        r = fn(url, params=params, timeout=TIMEOUT, headers=UA,
               allow_redirects=True, stream=not head)
        body = b"" if head else r.content[:200000]
        rec.update(
            status=r.status_code,
            seconds=round(time.time() - t, 2),
            bytes=len(body),
            content_type=r.headers.get("content-type", "")[:60],
            final_url=r.url[:160],
            sample=(body[:400].decode("utf-8", "replace") if body else ""),
        )
    except Exception as e:
        rec.update(status=None, seconds=round(time.time() - t, 2),
                   error=f"{type(e).__name__}: {str(e)[:200]}")
    return rec


def chemotype_census(smiles_list: list[str]) -> dict:
    """표본에 PT 서명이 얼마나 있는가 — 매칭 상대가 존재하는지 본다."""
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    pats = {
        "fluorine_ge3": (Chem.MolFromSmarts("[F]"), 3),
        "halogen_ge3": (Chem.MolFromSmarts("[F,Cl,Br,I]"), 3),
        "sulfone": (Chem.MolFromSmarts("[SX4](=O)(=O)"), 1),
        "nitrile": (Chem.MolFromSmarts("[C]#[N]"), 1),
        "aryl_ether": (Chem.MolFromSmarts("c-O-c"), 1),
    }
    hits = {k: 0 for k in pats}
    hits["F3_and_sulfone"] = 0
    hits["all_four"] = 0
    n = 0
    for smi in smiles_list:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        n += 1
        got = {k: len(m.GetSubstructMatches(p)) >= c for k, (p, c) in pats.items()}
        for k, v in got.items():
            hits[k] += bool(v)
        if got["fluorine_ge3"] and got["sulfone"]:
            hits["F3_and_sulfone"] += 1
        if got["fluorine_ge3"] and got["sulfone"] and got["nitrile"] and got["aryl_ether"]:
            hits["all_four"] += 1
    return {"n_parsed": n,
            "counts": hits,
            "pct": {k: round(100 * v / n, 2) if n else None for k, v in hits.items()}}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="zinc_probe.json")
    ap.add_argument("--sample-count", type=int, default=2000)
    args = ap.parse_args()

    results = []

    # 1. ZINC22 CartBlanche — 무작위 표본
    results.append(probe(
        "zinc22_random_druglike",
        "https://cartblanche22.docking.org/substance/random.txt",
        {"count": args.sample_count, "subset": "drug-like",
         "output_fields": "zinc_id,smiles"},
    ))

    # 2. ZINC20 — 구버전이 살아 있는지. .txt 가 HTML 을 돌려주므로 .csv 도 본다
    results.append(probe(
        "zinc20_random_csv",
        "https://zinc20.docking.org/substances/random.csv",
        {"count": 100, "output_fields": "zinc_id,smiles"},
    ))
    results.append(probe(
        "zinc20_subsets",
        "https://zinc20.docking.org/substances/subsets/drug-like.csv",
        {"count": 100},
    ))

    # 3. tranche 파일 서버가 서 있는지 (디렉터리 인덱스)
    results.append(probe("files_zinc22_index", "https://files.docking.org/zinc22/"))
    results.append(probe("files_root", "https://files.docking.org/"))

    # 4. 대안 — PubChem 은 대용량이고 공개 API 가 안정적이다
    results.append(probe(
        "pubchem_property_search",
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsubstructure/smarts/"
        "FC(F)(F)c1ccccc1/cids/TXT",
        {"MaxRecords": 50},
    ))

    for r in results:
        s = r.get("status")
        mark = "OK " if s == 200 else "   "
        print(f"{mark}{r['name']:28s} status={s} {r.get('bytes', 0):>8,d}B "
              f"{r['seconds']:>5.1f}s {r.get('error', '')[:60]}")
        if r.get("sample"):
            first = r["sample"].strip().splitlines()[:2]
            for line in first:
                print(f"      {line[:120]}")

    # 성공한 무작위 표본이 있으면 화학형을 센다
    census = None
    best = next((r for r in results
                 if r["name"].startswith("zinc") and r.get("status") == 200
                 and r.get("bytes", 0) > 1000), None)
    if best:
        smis = []
        for line in best["sample"].splitlines():
            parts = line.split()
            if len(parts) >= 2 and parts[1].count("c") + parts[1].count("C") > 0:
                smis.append(parts[1])
        # sample 은 앞 400바이트뿐이라 다시 받아 전체를 센다
        try:
            r = requests.get(best["url"], params=best["params"],
                             timeout=TIMEOUT, headers=UA)
            smis = [l.split()[1] for l in r.text.splitlines()[1:]
                    if len(l.split()) >= 2]
        except Exception as e:
            print(f"[!] 전체 표본 재수신 실패: {e}")
        if smis:
            census = chemotype_census(smis)
            print(f"\n[*] {best['name']} 표본 {census['n_parsed']}개의 화학형 분포")
            print("    (비교: EPAS1 활성 — 불소>=3 80%, 설폰 72%, 넷 모두 35%)")
            print("    (비교: ChEMBL 무작위 — 불소>=3 4.3%, 설폰 9.5%, 넷 모두 0%)")
            for k, v in census["pct"].items():
                shown = f"{v:6.2f}%" if v is not None else "   n/a"
                print(f"      {k:18s} {shown}  ({census['counts'][k]}개)")

    out = {"probes": results, "chemotype_census": census,
           "verdict": ("ZINC 접근 가능" if best else "ZINC 접근 불가 — 대안 필요")}
    Path(args.out).write_text(json.dumps(out, ensure_ascii=False, indent=2))
    print(f"\n[=] {out['verdict']}")
    print(f"    → {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
