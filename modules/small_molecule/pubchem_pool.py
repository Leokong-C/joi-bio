#!/usr/bin/env python3
"""
PubChem 에서 물성 매칭용 decoy 후보 풀을 모은다.

## 왜 PubChem 인가

1순위는 ZINC 였다. DUD-E 가 쓰고, tranche 가 HBD·logP·MW 로 미리 구간화돼
있어 우리가 필요한 물성 창을 그대로 집을 수 있다. 그런데 2026-10-25 탐침에서
**cartblanche22.docking.org 와 files.docking.org 가 GitHub Actions 에서
연결 타임아웃**이었다(45초). zinc20 은 200 을 주지만 TSV 가 아니라 HTML 을
돌려준다.

PubChem 은 같은 탐침에서 1.2초에 200 을 줬다. 화합물 약 1.2억 개이고
**SMARTS 부분구조 검색**이 공개 API 로 열려 있다 — 이게 ZINC tranche 보다
오히려 우리 문제에 맞는다. 무작위로 퍼내서 거르는 대신, 필요한 화학형을
직접 찾아 올 수 있기 때문이다.

## 무엇을 찾는가

ChEMBL 무작위 풀의 문제는 PT 서명을 가진 화합물이 0개였다는 것이다.
그래서 활성의 화학형을 SMARTS 로 적어 그것을 가진 화합물을 직접 찾는다.

  EPAS1 활성 (PT 계열)   불소>=3 80%, 설폰 72%, 나이트릴 45%, 다이아릴에테르 60%
  ChEMBL 무작위          4.3%, 9.5%, 1.4%, 1.7% — 넷 모두 0개

여기서 모은 풀을 decoys.select 에 넘기면 물성·작용기 매칭이 비로소 맞출
상대를 갖게 된다.

## 중요 — 이건 decoy 를 쉽게 만드는 게 아니다

활성과 **닮은** 화합물을 모으는 것이 목적이지만, decoys.select 가 ECFP4
Tanimoto 상한(기본 0.35)으로 토폴로지가 겹치는 것을 걸러낸다. 즉 "작용기는
같은데 골격은 다른" 화합물만 남는다. 그래야 불소 세기 같은 지름길이 막히고
도킹이 실제로 일을 해야 한다.

decoy 가 실제로 비활성인지는 알 수 없다 — 추정 비활성이다. 이 한계는
docs/context-of-use.md 에 적혀 있다.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import requests

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
UA = {"User-Agent": "joi-bio/1.0 (virtual screening benchmark construction)"}
TIMEOUT = 120
PAUSE = 0.25  # PubChem 권고 — 초당 5건 이하


def _get(url: str, params: dict | None = None, retries: int = 3) -> requests.Response | None:
    for i in range(retries):
        try:
            r = requests.get(url, params=params, timeout=TIMEOUT, headers=UA)
            if r.status_code == 200:
                return r
            if r.status_code in (503, 504):  # 혼잡 — 물러났다 다시
                time.sleep(2 ** i)
                continue
            print(f"    [!] {r.status_code} {url[:80]}")
            return None
        except Exception as e:
            print(f"    [!] {type(e).__name__} {str(e)[:120]}")
            time.sleep(2 ** i)
    return None


def smarts_cids(smarts: str, max_records: int, label: str) -> list[str]:
    """SMARTS 부분구조로 CID 를 찾는다.

    POST 로 보낸다. SMARTS 에 흔한 '#'(삼중결합)이 URL 에서는 조각 구분자라
    경로에 넣으면 잘려 400 이 난다 — 2026-10-25 aryl_nitrile_F 가 그랬다.
    """
    r = None
    for attempt in range(3):
        try:
            r = requests.post(
                f"{PUG}/compound/fastsubstructure/smarts/cids/TXT",
                data={"smarts": smarts, "MaxRecords": max_records},
                timeout=TIMEOUT, headers=UA,
            )
            if r.status_code == 200:
                break
            if r.status_code in (503, 504):
                time.sleep(2 ** attempt); continue
            print(f"    [!] {label}: {r.status_code} {r.text[:120]}")
            r = None
            break
        except Exception as e:
            print(f"    [!] {label}: {type(e).__name__} {str(e)[:100]}")
            time.sleep(2 ** attempt)
            r = None
    time.sleep(PAUSE)
    if r is None:
        return []
    cids = [l.strip() for l in r.text.splitlines() if l.strip().isdigit()]
    print(f"    {label:28s} CID {len(cids):,}개")
    return cids


SMILES_FIELDS = ["SMILES", "CanonicalSMILES", "ConnectivitySMILES", "IsomericSMILES"]


def _working_smiles_field() -> str:
    """PubChem 이 2025년에 CanonicalSMILES 를 바꿨다. 무엇이 사는지 직접 본다."""
    for f in SMILES_FIELDS:
        try:
            r = requests.get(f"{PUG}/compound/cid/2244/property/{f}/JSON",
                             timeout=30, headers=UA)
            if r.status_code == 200 and f in r.text:
                print(f"    SMILES 필드: {f}")
                return f
        except Exception:
            pass
        time.sleep(PAUSE)
    print("    [!] 쓸 수 있는 SMILES 필드를 못 찾았다 — SMILES 로 시도")
    return "SMILES"


def fetch_properties(cids: list[str], chunk: int = 150) -> list[dict]:
    """CID 묶음에서 SMILES 와 물성을 받는다."""
    sf = _working_smiles_field()
    want = f"{sf},MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount,RotatableBondCount"
    out: list[dict] = []
    for i in range(0, len(cids), chunk):
        part = cids[i : i + chunk]
        try:
            r = requests.post(f"{PUG}/compound/cid/property/{want}/JSON",
                              data={"cid": ",".join(part)},
                              timeout=TIMEOUT, headers=UA)
            if r.status_code != 200:
                print(f"    [!] 물성 {r.status_code}: {r.text[:140]}")
                r = None
        except Exception as e:
            print(f"    [!] 물성 {type(e).__name__}: {str(e)[:100]}")
            r = None
        time.sleep(PAUSE)
        if r is None:
            continue
        try:
            for p in r.json()["PropertyTable"]["Properties"]:
                smi = p.get(sf) or p.get("SMILES") or p.get("CanonicalSMILES")
                if not smi:
                    continue
                out.append(
                    {
                        "smiles": smi,
                        "chembl_id": f"PUBCHEM{p['CID']}",  # 파이프라인이 이 키를 쓴다
                        "pubchem_cid": p["CID"],
                        "pref_name": "",
                        "label": "decoy",
                        "source": "pubchem_smarts",
                        "protected": False,
                    }
                )
        except Exception as e:
            print(f"    [!] 응답 파싱 실패: {str(e)[:100]}")
        if (i // chunk) % 10 == 0:
            print(f"      {len(out):,}개 수집", flush=True)
    return out


def build_pool(
    smarts_spec: dict[str, str], per_pattern: int, verbose: bool = True
) -> tuple[list[dict], dict]:
    """각 SMARTS 마다 CID 를 모아 합치고 중복을 없앤다."""
    all_cids: dict[str, str] = {}
    per_label: dict[str, int] = {}
    for label, smarts in smarts_spec.items():
        cids = smarts_cids(smarts, per_pattern, label)
        per_label[label] = len(cids)
        for c in cids:
            all_cids.setdefault(c, label)

    if verbose:
        print(f"\n[*] 고유 CID {len(all_cids):,}개 — 물성 수신 시작")
    recs = fetch_properties(list(all_cids))
    for r in recs:
        r["matched_pattern"] = all_cids.get(str(r["pubchem_cid"]), "")

    stats = {
        "patterns": smarts_spec,
        "per_pattern_cids": per_label,
        "unique_cids": len(all_cids),
        "fetched": len(recs),
    }
    return recs, stats


# EPAS1 활성(PT 계열)의 화학형. 다른 타깃은 활성에서 뽑아 바꿔 쓴다.
EPAS1_PATTERNS = {
    "aryl_CF2_or_CF3": "c1ccccc1C(F)(F)",
    "aryl_sulfone": "c1ccccc1S(=O)(=O)",
    "aryl_nitrile_F": "N#Cc1cc(F)ccc1",
    "diaryl_ether_halo": "c1ccccc1Oc1cc(F)ccc1",
    "fluorinated_sulfone": "FC(F)S(=O)(=O)",
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-pattern", type=int, default=3000)
    ap.add_argument("--out", required=True)
    ap.add_argument("--stats-out", default=None)
    ap.add_argument("--patterns", default=None,
                    help="SMARTS JSON 파일. 없으면 EPAS1 기본값")
    args = ap.parse_args()

    spec = (
        json.loads(Path(args.patterns).read_text())
        if args.patterns
        else EPAS1_PATTERNS
    )
    print(f"[*] SMARTS {len(spec)}개, 패턴당 최대 {args.per_pattern:,}개")

    recs, stats = build_pool(spec, args.per_pattern)

    with open(args.out, "w") as f:
        for r in recs:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(f"\n[+] 풀 {len(recs):,}개 → {args.out}")

    if args.stats_out:
        Path(args.stats_out).write_text(json.dumps(stats, ensure_ascii=False, indent=2))
        print(f"    → {args.stats_out}")

    print("\n    다음: decoys.py 로 물성 매칭, 그다음 trivial_baselines.py --library 로 수용 확인")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
