#!/usr/bin/env python3
"""
core.structure.fetch — 타깃 구조 확보와 포켓 정의.

이 모듈의 판단이 파이프라인 전체의 상한을 결정한다. 포켓이 틀리면 뒤의 도킹·
검증·보고가 전부 무의미해지므로, 공결정 리간드가 있는 구조를 최우선으로 고른다.

전략
  1. RCSB에서 UniProt accession에 매핑된 X-ray 구조를 해상도순으로 조회
  2. 상위 후보를 내려받아 '약물유사 공결정 리간드'가 있는지 검사
  3. 리간드가 있는 구조 중 해상도가 가장 좋은 것을 선택
  4. 그 리간드를 품은 폴리머 체인만 남기고, 리간드 무게중심을 도킹 박스 중심으로
  5. 리간드가 있는 구조가 없으면 apo 구조 → blind docking (confidence=low)
  6. 실험 구조가 전혀 없으면 AlphaFold 모델 → blind docking (confidence=low)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import requests

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_FILE = "https://files.rcsb.org/download/{pdb_id}.pdb"
AFDB_TEMPLATES = [
    "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v4.pdb",
    "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v3.pdb",
    "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v2.pdb",
]
TIMEOUT = 90

# 용매·완충제·결정화 첨가물 — 리간드로 오인하면 포켓이 엉뚱한 곳에 잡힌다.
SOLVENT = {
    "HOH", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "2PE",
    "MPD", "ACT", "ACY", "FMT", "CIT", "TRS", "EPE", "MES", "BME", "DTT", "IMD",
    "DMS", "IPA", "ETA", "BU3", "P6G", "12P", "15P", "SIN", "TLA", "MLI", "NO3",
    "AZI", "CO3", "BCT", "FLC", "OXL", "MYR", "PLM", "OLA", "LDA", "C8E", "BOG",
    "SDS", "UNX", "UNL", "NH4", "CAC", "BTB", "NHE", "POL", "PE4", "PEU",
}
# 금속·보조인자 — 활성부위 구성 요소이므로 수용체에 남긴다(리간드로는 세지 않음).
COFACTOR_KEEP = {
    "ZN", "MG", "MN", "FE", "FE2", "CA", "CU", "CU1", "NI", "CO", "K", "NA",
    "CL", "BR", "IOD", "F", "CD", "HG", "SR", "CS", "RB", "YB", "GD",
}


def _post(url: str, payload: dict, retries: int = 3) -> dict | None:
    last = ""
    for attempt in range(retries):
        try:
            r = requests.post(url, json=payload, timeout=TIMEOUT)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 204:  # 조건에 맞는 구조 없음
                return {}
            last = f"HTTP {r.status_code}: {r.text[:300]}"
            if 400 <= r.status_code < 500:
                break  # 쿼리 자체가 잘못됨 — 재시도해도 같다
        except requests.RequestException as e:
            last = f"{type(e).__name__}: {e}"
        time.sleep(2 * (attempt + 1))
    if last:
        print(f"    [rcsb] {last}", file=sys.stderr)
    return None


def _get(url: str, retries: int = 3) -> requests.Response | None:
    last = ""
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=TIMEOUT)
            if r.status_code == 200:
                return r
            if r.status_code == 404:
                return None
            last = f"HTTP {r.status_code} {url}"
        except requests.RequestException as e:
            last = f"{type(e).__name__}: {e} ({url})"
        time.sleep(2 * (attempt + 1))
    if last:
        print(f"    [http] {last}", file=sys.stderr)
    return None


ACCESSION_ATTR = (
    "rcsb_polymer_entity_container_identifiers"
    ".reference_sequence_identifiers.database_accession"
)


def _accession_node(uniprot: str) -> dict:
    return {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": ACCESSION_ATTR,
            "operator": "exact_match",
            "value": uniprot,
        },
    }


def _resolution_node() -> dict:
    return {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": "rcsb_entry_info.resolution_combined",
            "operator": "less_or_equal",
            "value": 3.5,
        },
    }


def search_rcsb(uniprot: str, max_hits: int = 25) -> list[str]:
    """UniProt accession에 매핑된 실험 구조를 반환.

    RCSB Search API는 버전에 따라 sort 키와 필터 허용 범위가 달라진다.
    쿼리를 단순한 쪽으로 단계적으로 낮추며 재시도한다.
    """
    sort_asc = [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}]
    attempts = [
        {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [_accession_node(uniprot), _resolution_node()],
            },
            "return_type": "entry",
            "request_options": {
                "sort": sort_asc,
                "paginate": {"start": 0, "rows": max_hits},
            },
        },
        {
            "query": _accession_node(uniprot),
            "return_type": "entry",
            "request_options": {
                "sort": sort_asc,
                "paginate": {"start": 0, "rows": max_hits},
            },
        },
        {
            "query": _accession_node(uniprot),
            "return_type": "entry",
            "request_options": {"paginate": {"start": 0, "rows": max_hits}},
        },
    ]

    for i, payload in enumerate(attempts, 1):
        data = _post(RCSB_SEARCH, payload)
        if data:
            ids = [hit["identifier"] for hit in data.get("result_set", [])]
            if ids:
                print(f"    RCSB 검색 {len(ids)}건 (쿼리 형태 {i})")
                return ids
        if data == {}:
            print(f"    RCSB 검색 0건 (쿼리 형태 {i})")
    print("    [!] RCSB 검색 전부 실패")
    return []


def parse_pdb(text: str) -> dict:
    """PDB 텍스트에서 첫 모델의 원자와 해상도를 뽑는다."""
    resolution = None
    atoms, hetatms = [], []
    in_model = True
    for line in text.splitlines():
        if line.startswith("REMARK   2 RESOLUTION."):
            for tok in line.split():
                try:
                    resolution = float(tok)
                    break
                except ValueError:
                    continue
        rec = line[:6].strip()
        if rec == "MODEL":
            parts = line.split()
            in_model = len(parts) > 1 and parts[1] == "1"
            continue
        if rec == "ENDMDL":
            in_model = False
            continue
        if not in_model or rec not in ("ATOM", "HETATM"):
            continue
        try:
            rec_d = {
                "line": line,
                "name": line[12:16].strip(),
                "altloc": line[16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21].strip(),
                "resseq": line[22:27].strip(),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
                "element": line[76:78].strip().upper() or line[12:14].strip()[0],
            }
        except (ValueError, IndexError):
            continue
        (atoms if rec == "ATOM" else hetatms).append(rec_d)
    return {"resolution": resolution, "atoms": atoms, "hetatms": hetatms}


def find_druglike_ligands(hetatms: list[dict], min_heavy: int = 12) -> list[dict]:
    """공결정 리간드 후보를 무거운 것부터 반환. 용매·금속은 제외."""
    groups: dict[tuple, list[dict]] = {}
    for a in hetatms:
        rn = a["resname"]
        if rn in SOLVENT or rn in COFACTOR_KEEP:
            continue
        groups.setdefault((rn, a["chain"], a["resseq"]), []).append(a)

    out = []
    for (resname, chain, resseq), ats in groups.items():
        heavy = [a for a in ats if a["element"] != "H"]
        if len(heavy) < min_heavy:
            continue
        n = len(heavy)
        center = [sum(a[k] for a in heavy) / n for k in ("x", "y", "z")]
        out.append(
            {
                "resname": resname,
                "chain": chain,
                "resseq": resseq,
                "n_heavy": n,
                "center": [round(c, 3) for c in center],
                "atoms": heavy,
            }
        )
    out.sort(key=lambda g: -g["n_heavy"])
    return out


def evaluate_candidate(pdb_id: str) -> dict | None:
    """구조 하나를 내려받아 리간드 유무와 해상도로 점수를 매긴다."""
    r = _get(RCSB_FILE.format(pdb_id=pdb_id))
    if r is None:
        return None
    parsed = parse_pdb(r.text)
    if not parsed["atoms"]:
        return None
    ligands = find_druglike_ligands(parsed["hetatms"])
    return {
        "pdb_id": pdb_id,
        "text": r.text,
        "resolution": parsed["resolution"],
        "parsed": parsed,
        "ligands": ligands,
        "has_ligand": bool(ligands),
    }


def select_structure(uniprot: str, max_download: int = 12) -> dict | None:
    """리간드를 품은 최고 해상도 구조를 고른다. 없으면 최고 해상도 apo."""
    pdb_ids = search_rcsb(uniprot)
    if not pdb_ids:
        return None

    holo, apo = [], []
    for pdb_id in pdb_ids[:max_download]:
        cand = evaluate_candidate(pdb_id)
        if cand is None:
            continue
        (holo if cand["has_ligand"] else apo).append(cand)
        # 리간드 구조를 3개 확보하면 충분 — 해상도순 조회라 이미 상위권이다.
        if len(holo) >= 3:
            break

    pool = holo or apo
    if not pool:
        return None
    pool.sort(key=lambda c: (c["resolution"] if c["resolution"] else 99.0))
    return pool[0]


def build_receptor(cand: dict, out_pdb: Path) -> dict:
    """리간드를 품은 체인만 남겨 수용체 PDB를 쓴다. 금속 보조인자는 유지."""
    parsed = cand["parsed"]
    ligands = cand["ligands"]

    if ligands:
        lig = ligands[0]
        keep_chain = lig["chain"]
    else:
        lig = None
        chains = {}
        for a in parsed["atoms"]:
            chains[a["chain"]] = chains.get(a["chain"], 0) + 1
        keep_chain = max(chains, key=chains.get) if chains else ""

    lines = []
    for a in parsed["atoms"]:
        if a["chain"] != keep_chain:
            continue
        if a["altloc"] not in ("", "A"):  # 대체 위치는 A만
            continue
        lines.append(a["line"])

    # 선택한 체인 근처의 금속/보조인자 유지 (8 Å 이내)
    prot_xyz = [
        (a["x"], a["y"], a["z"])
        for a in parsed["atoms"]
        if a["chain"] == keep_chain
    ]
    kept_cofactors = []
    for a in parsed["hetatms"]:
        if a["resname"] not in COFACTOR_KEEP:
            continue
        if any(
            (a["x"] - p[0]) ** 2 + (a["y"] - p[1]) ** 2 + (a["z"] - p[2]) ** 2 < 64
            for p in prot_xyz
        ):
            lines.append(a["line"])
            kept_cofactors.append(a["resname"])

    out_pdb.write_text("\n".join(lines) + "\nEND\n")

    return {
        "chain": keep_chain,
        "n_atoms": len(lines),
        "cofactors": sorted(set(kept_cofactors)),
        "cocrystal_ligand": (
            {
                "resname": lig["resname"],
                "n_heavy": lig["n_heavy"],
                "center": lig["center"],
            }
            if lig
            else None
        ),
    }


def define_pocket(cand: dict, receptor_meta: dict, out_pdb: Path) -> dict:
    """도킹 박스를 정의한다. 공결정 리간드 기준이 가장 신뢰도 높다."""
    ligands = cand["ligands"]
    if ligands:
        lig = ligands[0]
        center = lig["center"]
        spans = [
            max(abs(a[k] - center[i]) for a in lig["atoms"])
            for i, k in enumerate(("x", "y", "z"))
        ]
        size = [round(max(20.0, min(30.0, s * 2 + 9.0)), 1) for s in spans]
        return {
            "mode": "cocrystal_ligand",
            "reference_ligand": lig["resname"],
            "center": center,
            "size": size,
            "confidence": "high",
            "note": f"{lig['resname']} ({lig['n_heavy']} heavy atoms) 기준",
        }

    # apo / 예측 모델 — 단백질 무게중심에 큰 박스 (blind docking)
    coords = []
    for line in out_pdb.read_text().splitlines():
        if line[:6].strip() in ("ATOM", "HETATM"):
            try:
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                continue
    n = len(coords) or 1
    center = [round(sum(c[i] for c in coords) / n, 3) for i in range(3)]
    return {
        "mode": "blind_centroid",
        "reference_ligand": None,
        "center": center,
        "size": [30.0, 30.0, 30.0],
        "confidence": "low",
        "note": "공결정 리간드 없음 — blind docking. 결과 해석 주의.",
    }


def prepare_target(target: dict, outdir: Path) -> dict:
    """타깃 하나의 구조·포켓을 준비하고 메타데이터를 반환."""
    outdir.mkdir(parents=True, exist_ok=True)
    gene, uniprot = target["gene"], target["uniprot"]
    meta = {"gene": gene, "uniprot": uniprot, "declared_pocket_mode": target.get("pocket_mode")}

    cand = select_structure(uniprot)

    if cand is None:
        print("    실험 구조 확보 실패 → AlphaFold 모델 시도")
        r = None
        for tpl in AFDB_TEMPLATES:
            r = _get(tpl.format(acc=uniprot))
            if r is not None:
                break
        if r is None:
            meta["error"] = "실험 구조·AlphaFold 모델 모두 확보 실패"
            return meta
        (outdir / "raw.pdb").write_text(r.text)
        parsed = parse_pdb(r.text)
        cand = {
            "pdb_id": None,
            "text": r.text,
            "resolution": None,
            "parsed": parsed,
            "ligands": [],
            "has_ligand": False,
        }
        meta["source"] = "AlphaFoldDB"
        meta["warning"] = "예측 모델 — 측쇄 배치가 실제와 다를 수 있음"
    else:
        (outdir / "raw.pdb").write_text(cand["text"])
        meta["source"] = "RCSB"
        meta["pdb_id"] = cand["pdb_id"]
        meta["resolution"] = cand["resolution"]

    receptor_pdb = outdir / "receptor.pdb"
    meta["receptor"] = build_receptor(cand, receptor_pdb)
    meta["pocket"] = define_pocket(cand, meta["receptor"], receptor_pdb)
    meta["benchmark_ligands"] = target.get("benchmark_ligands", [])
    if target.get("risk"):
        meta["declared_risk"] = target["risk"]
    return meta


def main() -> int:
    import argparse

    import yaml

    root = Path(__file__).resolve().parents[2]
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", help="단일 타깃 (미지정 시 전체)")
    ap.add_argument("--config", default=str(root / "config" / "targets.yaml"))
    ap.add_argument("--outroot", default=str(root / "targets"))
    args = ap.parse_args()

    targets = yaml.safe_load(Path(args.config).read_text())["targets"]
    if args.gene:
        targets = [t for t in targets if t["gene"].upper() == args.gene.upper()]
        if not targets:
            print(f"[!] 알 수 없는 유전자: {args.gene}", file=sys.stderr)
            return 1

    failures = 0
    for t in targets:
        outdir = Path(args.outroot) / t["gene"]
        print(f"[*] {t['gene']} ({t['uniprot']}) 구조 준비")
        meta = prepare_target(t, outdir)
        (outdir / "target_meta.json").write_text(
            json.dumps(meta, indent=2, ensure_ascii=False)
        )

        if meta.get("error"):
            print(f"    [!] {meta['error']}")
            failures += 1
            continue
        p, rc = meta["pocket"], meta["receptor"]
        print(
            f"    source={meta['source']} pdb={meta.get('pdb_id')} "
            f"res={meta.get('resolution')} chain={rc['chain']} atoms={rc['n_atoms']}"
        )
        print(
            f"    pocket={p['mode']} ({p['confidence']}) "
            f"center={p['center']} size={p['size']}"
        )
        print(f"    {p['note']}")
        if rc["cofactors"]:
            print(f"    cofactors: {rc['cofactors']}")

    return 1 if failures == len(targets) else 0


if __name__ == "__main__":
    raise SystemExit(main())
