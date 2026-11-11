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
RCSB_CHEMCOMP = "https://data.rcsb.org/rest/v1/core/chemcomp/"
RCSB_NONPOLY = "https://data.rcsb.org/rest/v1/core/nonpolymer_entity/"

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


_MW: dict[str, float | None] = {}


def ligand_mw(comp: str) -> float | None:
    """화학성분의 분자량. 이온·완충제를 이름으로만 거르면 샌다 —
    11-11 1차 탐침에서 나트륨(NA)이 리간드로 세어져 오렉신1 의 홀로 수가
    부풀었다. 선언만 해 두고 쓰지 않은 하한을 실제로 적용한다."""
    if comp in _MW:
        return _MW[comp]
    try:
        r = requests.get(RCSB_CHEMCOMP + comp, timeout=30)
        r.raise_for_status()
        _MW[comp] = r.json().get("chem_comp", {}).get("formula_weight")
    except Exception:
        _MW[comp] = None
    return _MW[comp]


def entry_detail(pdb: str) -> dict:
    r = requests.get(RCSB_ENTRY + pdb, timeout=60)
    r.raise_for_status()
    j = r.json()
    info = j.get("rcsb_entry_info", {})
    res = (info.get("resolution_combined") or [None])[0]

    # experimental_method 는 문자열이고 experimental_method_count 와 짝이다.
    # 리스트로 보고 join 하면 글자가 쉼표로 쪼개진다 (1차 탐침의 'X,-,r,a,y')
    method = info.get("experimental_method")
    if not isinstance(method, str):
        method = ", ".join(e.get("method", "") for e in j.get("exptl", []) or [])

    # 리간드 성분 ID 를 모은다. nonpolymer_bound_components 는 비어 있는
    # 항목이 많아 (11-11 2차 탐침에서 KCNQ2 58개가 전부 0으로 나왔다)
    # 비중합체 엔티티를 직접 훑는 쪽을 1순위로 둔다.
    comps: list[str] = []
    ids = (j.get("rcsb_entry_container_identifiers", {})
            .get("non_polymer_entity_ids") or [])
    for eid in ids:
        try:
            e = requests.get(f"{RCSB_NONPOLY}{pdb}/{eid}", timeout=30).json()
            c = (e.get("pdbx_entity_nonpoly", {}) or {}).get("comp_id")
            if c:
                comps.append(c)
        except Exception:
            pass
    if not comps:
        comps = list(info.get("nonpolymer_bound_components") or [])

    ligs = []
    for c in dict.fromkeys(comps):
        if c in BUFFERS:
            continue
        mw = ligand_mw(c)
        if mw is None or mw < MIN_LIGAND_MW:
            continue
        ligs.append(c)
    return {"pdb_id": pdb, "resolution": res, "method": method or "",
            "n_nonpoly_entities": len(ids), "ligands": ligs}


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
    xray = [e for e in holo if "X-RAY" in (e["method"] or "").upper()]
    em = [e for e in holo if "ELECTRON" in (e["method"] or "").upper()]
    good = [e for e in holo if e["resolution"] and e["resolution"] <= 3.0]

    out.update({
        "n_examined": len(entries), "n_fetch_errors": errs,
        "n_holo": len(holo), "n_holo_xray": len(xray), "n_holo_em": len(em),
        "n_holo_le_3A": len(good),
        "best": sorted(
            ({"pdb_id": e["pdb_id"], "resolution": e["resolution"],
              "method": e["method"], "ligands": e["ligands"][:4]} for e in holo),
            key=lambda e: (e["resolution"] is None, e["resolution"] or 99))[:8],
    })
    out["c4_pass"] = out["n_holo"] >= 2
    print(f"    홀로 {out['n_holo']}개 (X선 {out['n_holo_xray']}, cryo-EM "
          f"{out['n_holo_em']}, 3.0 A 이하 {out['n_holo_le_3A']}) "
          f"→ C4 {'통과' if out['c4_pass'] else '실패'}")
    for e in out["best"][:4]:
        print(f"      {e['pdb_id']}  {e['resolution']} A  {e['method'][:18]:18s} "
              f"{','.join(e['ligands'])}")
    return out


# 리간드가 붙어 있다고 알려진 항목. 여기서 0 이 나오면 탐침이 고장난 것이지
# 구조가 없는 것이 아니다. 11-11 에 두 번 연속 조용히 틀린 수를 냈다.
SELFTEST = {"6WJC": "CHRM1 아세틸콜린", "1RFF": "TDP1 바나듐산",
            "7CR0": "KCNQ2 리간드 결합"}


def selftest() -> bool:
    print("[*] 자체 점검 — 리간드가 있다고 알려진 항목에서 리간드가 나오는가")
    ok = True
    for pdb, what in SELFTEST.items():
        try:
            e = entry_detail(pdb)
            found = bool(e["ligands"])
            print(f"    {pdb} ({what}): 비중합체 엔티티 "
                  f"{e.get('n_nonpoly_entities', '?')}개 → "
                  f"리간드 {e['ligands'] or '없음'}  {'OK' if found else '[!] 0'}")
            ok &= found
        except Exception as ex:
            print(f"    {pdb}: 조회 실패 {type(ex).__name__}: {ex}")
            ok = False
    if not ok:
        print("    [!] 점검 실패 — 아래 숫자를 믿지 말 것")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", required=True, help="사람 유전자 기호, 쉼표 구분")
    ap.add_argument("--out", default="structures.json")
    args = ap.parse_args()

    healthy = selftest()

    res = []
    for g in [s.strip() for s in args.genes.split(",") if s.strip()]:
        try:
            res.append(probe(g))
        except Exception as e:
            print(f"    [!] {g} 실패: {type(e).__name__}: {e}")
            res.append({"gene": g, "error": f"{type(e).__name__}: {e}"})

    with open(args.out, "w") as f:
        json.dump({"selftest_passed": healthy, "targets": res},
                  f, ensure_ascii=False, indent=2)

    print(f"\n{'유전자':10s} {'UniProt':10s} {'PDB':>5s} {'홀로':>5s} "
          f"{'X선':>5s} {'cryoEM':>7s} {'<=3A':>5s}  C4")
    print("-" * 70)
    for r in res:
        if "error" in r:
            print(f"{r['gene']:10s} {'—':10s} {'—':>5s} {'—':>5s} {'—':>5s} "
                  f"{'—':>7s} {'—':>5s}  측정 실패")
            continue
        print(f"{r['gene']:10s} {r['accession']:10s} {r['n_pdb_total']:5d} "
              f"{r['n_holo']:5d} {r['n_holo_xray']:5d} {r['n_holo_em']:7d} "
              f"{r['n_holo_le_3A']:5d}  {'통과' if r['c4_pass'] else '실패'}")
    print(f"→ {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
