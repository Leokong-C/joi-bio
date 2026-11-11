#!/usr/bin/env python3
"""
타깃의 **도킹 가능한 구조**가 있는지 두드려 본다. 선정 기준 C4.

## 왜 별도 스크립트인가

선별(C1~C3)을 통과해도 구조가 없으면 못 쓴다. 그리고 우리 설계는 구조
하나로 끝나지 않는다 — **치료군·대조군 양팔**이 필요하므로 리간드가 붙은
홀로 구조가 최소 2개 있어야 한다.

아포(리간드 없음) 구조만 있는 타깃은 포켓을 정의할 기준이 없다. 10-16 에
구조 선택이 상위권 enrichment 의 75% 를 움직인다는 것을 쟀으므로, 어떤
구조가 있는지는 부수 정보가 아니라 **결과를 좌우하는 변수**다.

## 무엇을 세는가

- 홀로 구조 수 — 분자량 150 이상의 비완충 리간드가 붙은 것
- 해상도 분포 — X선 기준. 저해상도만 있으면 측쇄 배치를 못 믿는다
- 실험 방법 — cryo-EM 만 있는 막단백질은 포켓 정의가 어렵다

Claude 컨테이너는 UniProt·RCSB 가 프록시 403 이라 **Actions 에서만** 된다.
"""
from __future__ import annotations

import argparse
import json
import sys

import requests

UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/"

# 결정화·동결 보호제. 분자량만으로는 안 걸러진다
BUFFERS = {
    "HOH", "SO4", "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "2PE", "P6G",
    "ACT", "FMT", "CIT", "FLC", "TRS", "MES", "EPE", "DMS", "IMD", "MPD",
    "PO4", "NO3", "ACY", "BME", "DTT", "TLA", "MLI", "BTB", "NHE", "SIN",
    "OLC", "OLA", "PLM", "STE", "CLR", "Y01", "HEX", "LDA", "D10", "12P",
    "C8E", "LMT", "DD9", "UNL", "UNX",
}
MIN_LIGAND_MW = 150.0


def uniprot_accession(gene: str) -> dict:
    r = requests.get(UNIPROT, timeout=60, params={
        "query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
        "fields": "accession,protein_name", "format": "json", "size": 5})
    r.raise_for_status()
    hits = r.json().get("results", [])
    if not hits:
        return {}
    h = hits[0]
    name = (h.get("proteinDescription", {}).get("recommendedName", {})
             .get("fullName", {}).get("value", ""))
    return {"accession": h["primaryAccession"], "protein": name,
            "n_reviewed_hits": len(hits)}


def pdb_ids(acc: str) -> list[str]:
    q = {
        "query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers."
                         "reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": acc}},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 500}},
    }
    r = requests.post(RCSB_SEARCH, json=q, timeout=90)
    if r.status_code == 204:          # 검색 결과 없음은 오류가 아니다
        return []
    r.raise_for_status()
    return [h["identifier"] for h in r.json().get("result_set", [])]


def entry_detail(pdb: str) -> dict:
    r = requests.get(RCSB_ENTRY + pdb, timeout=60)
    r.raise_for_status()
    j = r.json()
    res = (j.get("rcsb_entry_info", {}).get("resolution_combined") or [None])[0]
    method = ",".join(j.get("rcsb_entry_info", {}).get("experimental_method", "")
                      if isinstance(j.get("rcsb_entry_info", {})
                                    .get("experimental_method"), str)
                      else j.get("rcsb_entry_info", {}).get("experimental_method", []) or [])
    ligs = [c for c in (j.get("rcsb_entry_info", {})
                        .get("nonpolymer_bound_components") or [])
            if c not in BUFFERS]
    return {"pdb_id": pdb, "resolution": res,
            "method": method or j.get("exptl", [{}])[0].get("method", ""),
            "ligands": ligs}


def probe(gene: str) -> dict:
    out: dict = {"gene": gene}
    u = uniprot_accession(gene)
    if not u:
        out["error"] = "UniProt 에서 사람 reviewed 항목을 못 찾았다"
        return out
    out.update(u)
    ids = pdb_ids(u["accession"])
    out["n_pdb_total"] = len(ids)
    print(f"[*] {gene} ({u['accession']}) — PDB {len(ids)}개: {u['protein'][:60]}")

    entries, errs = [], 0
    for p in ids[:120]:               # 상한을 둔다. 전수가 필요한 단계가 아니다
        try:
            entries.append(entry_detail(p))
        except Exception:
            errs += 1
    holo = [e for e in entries if e["ligands"]]
    xray = [e for e in holo if e["resolution"] and "X-RAY" in (e["method"] or "").upper()]
    good = [e for e in xray if e["resolution"] <= 3.0]

    out.update({
        "n_examined": len(entries), "n_fetch_errors": errs,
        "n_holo": len(holo), "n_holo_xray": len(xray),
        "n_holo_xray_le_3A": len(good),
        "best": sorted(
            ({"pdb_id": e["pdb_id"], "resolution": e["resolution"],
              "method": e["method"], "ligands": e["ligands"][:4]} for e in holo),
            key=lambda e: (e["resolution"] is None, e["resolution"] or 99))[:8],
    })
    out["c4_pass"] = out["n_holo"] >= 2
    print(f"    홀로 {out['n_holo']}개 (X선 {out['n_holo_xray']}, 3.0 A 이하 "
          f"{out['n_holo_xray_le_3A']}) → C4 {'통과' if out['c4_pass'] else '실패'}")
    for e in out["best"][:4]:
        print(f"      {e['pdb_id']}  {e['resolution']} A  {e['method'][:18]:18s} "
              f"{','.join(e['ligands'])}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", required=True, help="사람 유전자 기호, 쉼표 구분")
    ap.add_argument("--out", default="structures.json")
    args = ap.parse_args()

    res = []
    for g in [s.strip() for s in args.genes.split(",") if s.strip()]:
        try:
            res.append(probe(g))
        except Exception as e:
            print(f"    [!] {g} 실패: {type(e).__name__}: {e}")
            res.append({"gene": g, "error": f"{type(e).__name__}: {e}"})

    with open(args.out, "w") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)

    print(f"\n{'유전자':10s} {'UniProt':10s} {'PDB':>5s} {'홀로':>5s} "
          f"{'X선':>5s} {'<=3A':>5s}  C4")
    print("-" * 62)
    for r in res:
        if "error" in r:
            print(f"{r['gene']:10s} {'—':10s} {'—':>5s} {'—':>5s} {'—':>5s} "
                  f"{'—':>5s}  측정 실패")
            continue
        print(f"{r['gene']:10s} {r['accession']:10s} {r['n_pdb_total']:5d} "
              f"{r['n_holo']:5d} {r['n_holo_xray']:5d} {r['n_holo_xray_le_3A']:5d}  "
              f"{'통과' if r['c4_pass'] else '실패'}")
    print(f"→ {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
