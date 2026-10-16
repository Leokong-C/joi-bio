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
import re
import sys
import time
from pathlib import Path

import requests

from core import chem

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
            # "REMARK   2 RESOLUTION.    1.60 ANGSTROMS."
            # 단순 split 으로 첫 숫자를 집으면 REMARK 번호 2 를 해상도로 읽는다.
            # EPAS1 1차 실행에서 모든 후보가 2.0 Å 로 나온 원인.
            m = re.search(r"RESOLUTION\.\s+([0-9]*\.?[0-9]+)", line)
            if m:
                try:
                    resolution = float(m.group(1))
                except ValueError:
                    pass
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


def score_candidate(cand: dict, bench_smiles: list[str]) -> dict:
    """후보 구조의 공결정 리간드가 벤치마크 화합물과 얼마나 닮았는지."""
    info = {"similarity": 0.0, "ligand_smiles": None, "matched_benchmark": None}
    if not cand["ligands"]:
        return info
    lig = cand["ligands"][0]
    smi = chem.rcsb_ligand_smiles(lig["resname"])
    info["ligand_smiles"] = smi
    if smi and bench_smiles:
        sim, idx = chem.max_similarity(smi, bench_smiles)
        info["similarity"] = sim
        info["matched_benchmark"] = idx
    return info


def select_structure(
    uniprot: str,
    bench_smiles: list[str] | None = None,
    max_download: int = 16,
    force_pdb_id: str | None = None,
) -> dict | None:
    """벤치마크 화합물과 가장 닮은 리간드가 결합한 구조를 고른다.

    해상도만 보고 고르면 포켓이 다른 모양으로 닫혀 있는 구조를 집을 수 있다.
    리지드 도킹은 측쇄 재배열을 모델링하지 못하므로, 벤치마크와 닮은 리간드가
    이미 들어앉아 있는 구조라야 그 화합물을 제대로 수용한다.
    """
    bench_smiles = bench_smiles or []
    if force_pdb_id:
        # 음성 대조군: 벤치마크와 닮지 않은 구조를 일부러 지정한다.
        # 유사도 선택이 만드는 편향의 크기를 재기 위한 것이다.
        print(f"    구조 강제 지정: {force_pdb_id}")
        cand = evaluate_candidate(force_pdb_id.upper())
        if cand is None:
            return None
        cand["match"] = score_candidate(cand, bench_smiles)
        print(f"      유사도={cand['match']['similarity']:.3f} res={cand['resolution']}")
        return cand

    pdb_ids = search_rcsb(uniprot)
    if not pdb_ids:
        return None

    holo, apo = [], []
    for pdb_id in pdb_ids[:max_download]:
        cand = evaluate_candidate(pdb_id)
        if cand is None:
            continue
        if not cand["has_ligand"]:
            apo.append(cand)
            continue
        cand["match"] = score_candidate(cand, bench_smiles)
        holo.append(cand)
        sim = cand["match"]["similarity"]
        lig = cand["ligands"][0]
        print(
            f"      {pdb_id}: {lig['resname']} ({lig['n_heavy']} heavy) "
            f"res={cand['resolution']} 유사도={sim:.3f}"
        )
        if sim >= 0.45:  # 충분히 닮았다 — 더 볼 필요 없다
            break
        if len(holo) >= 8:
            break

    if not holo:
        if not apo:
            return None
        apo.sort(key=lambda c: c["resolution"] or 99.0)
        apo[0]["match"] = {"similarity": 0.0, "ligand_smiles": None, "matched_benchmark": None}
        return apo[0]

    # 유사도 구간으로 먼저 나누고, 같은 구간 안에서 해상도로 고른다.
    # 유사도가 조금 높다고 해상도 3.5 Å 구조를 집는 일을 막는다.
    def bucket(c: dict) -> int:
        sim = c["match"]["similarity"]
        return 0 if sim >= 0.45 else 1 if sim >= 0.25 else 2

    # 같은 구간·같은 해상도면 더 닮은 리간드를 고른다. 해상도를 0.1 Å 단위로
    # 뭉개야 2.0 과 2.0 이 동점으로 잡히고 유사도가 결정권을 갖는다.
    holo.sort(
        key=lambda c: (
            bucket(c),
            round(c["resolution"] or 99.0, 1),
            -c["match"]["similarity"],
        )
    )
    return holo[0]


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


def prepare_target(target: dict, outdir: Path, force_pdb_id: str | None = None) -> dict:
    """타깃 하나의 구조·포켓을 준비하고 메타데이터를 반환."""
    outdir.mkdir(parents=True, exist_ok=True)
    gene, uniprot = target["gene"], target["uniprot"]
    meta = {"gene": gene, "uniprot": uniprot, "declared_pocket_mode": target.get("pocket_mode")}

    # 벤치마크 화합물의 SMILES를 먼저 확보한다. 구조 선택의 기준이 된다.
    bench_names = target.get("benchmark_ligands", [])
    bench_smiles, bench_found = [], []
    for name in bench_names:
        rec = chem.chembl_by_name(name)
        if rec:
            bench_smiles.append(rec["smiles"])
            bench_found.append(name)
    meta["benchmark_smiles_found"] = bench_found
    if bench_names and not bench_smiles:
        print("    [!] 벤치마크 SMILES 확보 실패 — 해상도만으로 구조를 고른다")
    else:
        print(f"    벤치마크 SMILES {len(bench_smiles)}개 확보 → 구조 선택 기준으로 사용")

    cand = select_structure(uniprot, bench_smiles, force_pdb_id=force_pdb_id)

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
        meta["ligand_match"] = cand.get("match", {})
        meta["structure_selection"] = "forced" if force_pdb_id else "similarity"

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
    ap.add_argument(
        "--pdb-id",
        default=None,
        help="구조를 강제 지정 (음성 대조군용). 유사도 선택을 건너뛴다",
    )
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
        meta = prepare_target(t, outdir, force_pdb_id=args.pdb_id)
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
