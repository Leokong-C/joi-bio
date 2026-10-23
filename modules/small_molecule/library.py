#!/usr/bin/env python3
"""
modules.small_molecule.library — 스크리닝 라이브러리 구성.

세 종류를 섞는다.
  1. 벤치마크 리간드 (승인약·임상 화합물) — 이름으로 명시 수집. 반드시 포함된다.
  2. ChEMBL 활성 화합물 — 해당 타깃에 10 uM 이하 활성이 보고된 것
  3. Decoy — 무관한 약물유사 화합물

1·2가 양성 대조군이고, 이들이 도킹 상위에 올라오는지로 파이프라인을 검증한다.
벤치마크 리간드는 RDKit 필터에서 면제한다. 검증의 기준점을 필터가 지워버리면
시험 자체가 성립하지 않기 때문이다.
"""
from __future__ import annotations

import argparse
import json
import random
import sys
import time
from pathlib import Path

import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, FilterCatalog

RDLogger.DisableLog("rdApp.*")

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
TIMEOUT = 120
ACTIVITY_CUTOFF_NM = 10000.0  # 10 uM 이하만 활성으로 인정
SCREEN_SEED = 1023  # 스크리닝 전용 집합을 평가셋과 다른 씨드로 뽑는다

_params = FilterCatalog.FilterCatalogParams()
_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
PAINS = FilterCatalog.FilterCatalog(_params)


def _get(url: str, params: dict, retries: int = 3) -> dict | None:
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=TIMEOUT)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 404:
                return None
        except requests.RequestException:
            pass
        time.sleep(2 * (attempt + 1))
    return None


def fetch_benchmark_ligands(names: list[str]) -> list[dict]:
    """승인약·임상 화합물을 이름으로 찾아 명시적으로 수집한다 (spike-in)."""
    out = []
    for name in names:
        rec = None
        # pref_name 정확 일치 → 동의어 검색 순
        for params in (
            {"pref_name__iexact": name, "limit": 5},
            {"molecule_synonyms__molecule_synonym__icontains": name, "limit": 5},
        ):
            data = _get(f"{CHEMBL}/molecule.json", params)
            mols = (data or {}).get("molecules", [])
            for m in mols:
                smi = (m.get("molecule_structures") or {}).get("canonical_smiles")
                if smi:
                    rec = {
                        "smiles": smi,
                        "chembl_id": m["molecule_chembl_id"],
                        "pref_name": m.get("pref_name") or name,
                        "query_name": name,
                        "label": "known_active",
                        "source": "benchmark_spikein",
                        "protected": True,  # RDKit 필터 면제
                    }
                    break
            if rec:
                break
        if rec:
            out.append(rec)
            print(f"    [+] 벤치마크 확보: {name} → {rec['chembl_id']}")
        else:
            print(f"    [!] 벤치마크 미확보: {name} (ChEMBL에서 찾지 못함)")
    return out


def fetch_chembl_actives(uniprot: str, limit: int = 1500) -> list[dict]:
    """UniProt → ChEMBL target → 활성 데이터."""
    data = _get(
        f"{CHEMBL}/target.json",
        {"target_components__accession": uniprot, "limit": 20},
    )
    tids = [t["target_chembl_id"] for t in (data or {}).get("targets", [])]
    if not tids:
        return []
    print(f"    ChEMBL 타깃: {tids[:5]}")

    rows, offset = [], 0
    while len(rows) < limit and offset < 10000:
        page = _get(
            f"{CHEMBL}/activity.json",
            {
                "target_chembl_id__in": ",".join(tids[:5]),
                "standard_type__in": "IC50,Ki,Kd,EC50",
                "standard_relation": "=",
                "limit": 1000,
                "offset": offset,
            },
        )
        acts = (page or {}).get("activities", [])
        if not acts:
            break
        for a in acts:
            smi = a.get("canonical_smiles")
            val = a.get("standard_value")
            if not smi or val is None:
                continue
            try:
                val = float(val)
            except (TypeError, ValueError):
                continue
            if (a.get("standard_units") or "nM") != "nM" or val > ACTIVITY_CUTOFF_NM:
                continue
            rows.append(
                {
                    "smiles": smi,
                    "chembl_id": a["molecule_chembl_id"],
                    "pref_name": "",
                    "activity_type": a.get("standard_type"),
                    "activity_nM": val,
                    "label": "known_active",
                    "source": "chembl_activity",
                    "protected": False,
                }
            )
        offset += 1000
    return rows[:limit]


def annotate_names(rows: list[dict], batch: int = 50) -> None:
    """활성 화합물에 pref_name을 붙인다 (벤치마크 회수 판정에 필요)."""
    ids = sorted({r["chembl_id"] for r in rows if not r.get("pref_name")})
    names: dict[str, str] = {}
    for i in range(0, len(ids), batch):
        chunk = ids[i : i + batch]
        data = _get(
            f"{CHEMBL}/molecule.json",
            {"molecule_chembl_id__in": ",".join(chunk), "limit": batch},
        )
        for m in (data or {}).get("molecules", []):
            if m.get("pref_name"):
                names[m["molecule_chembl_id"]] = m["pref_name"]
    for r in rows:
        if not r.get("pref_name"):
            r["pref_name"] = names.get(r["chembl_id"], "")


def fetch_decoys(n: int, seed: int = 42) -> list[dict]:
    """무관한 약물유사 화합물 = decoy."""
    random.seed(seed)
    out, offset = [], random.randint(0, 20000)
    while len(out) < n:
        data = _get(
            f"{CHEMBL}/molecule.json",
            {
                "molecule_properties__mw_freebase__lte": 550,
                "molecule_properties__mw_freebase__gte": 220,
                "max_phase__isnull": "true",
                "limit": 1000,
                "offset": offset,
            },
        )
        mols = (data or {}).get("molecules", [])
        if not mols:
            break
        for m in mols:
            smi = (m.get("molecule_structures") or {}).get("canonical_smiles")
            if smi:
                out.append(
                    {
                        "smiles": smi,
                        "chembl_id": m["molecule_chembl_id"],
                        "pref_name": m.get("pref_name") or "",
                        "label": "decoy",
                        "source": "chembl_random",
                        "protected": False,
                    }
                )
        offset += 1000
    return out[:n]


def passes_filters(smiles: str) -> tuple[bool, str]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "unparseable"
    if Chem.GetFormalCharge(mol) not in range(-3, 4):
        return False, "charge"
    mw = Descriptors.MolWt(mol)
    if not 150 <= mw <= 650:
        return False, "mw"
    logp = Descriptors.MolLogP(mol)
    if not -2.0 <= logp <= 6.0:
        return False, "logp"
    if Descriptors.NumHDonors(mol) > 6 or Descriptors.NumHAcceptors(mol) > 12:
        return False, "hbond"
    if Descriptors.NumRotatableBonds(mol) > 12:
        return False, "rotbonds"
    if Descriptors.TPSA(mol) > 180:
        return False, "tpsa"
    if mol.GetNumHeavyAtoms() < 12:
        return False, "too_small"
    if PAINS.HasMatch(mol):
        return False, "PAINS"
    return True, "ok"


def canonical(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol else None


def cap_actives(lib: list[dict], cap: int) -> tuple[list[dict], int]:
    """활성을 가장 강력한 것 위주로 cap 개까지만 남긴다.

    활성이 라이브러리의 24%를 차지하면 EF1% 의 이론적 최대가 4.1 로 떨어져
    임계값 5.0 을 어떤 성능으로도 넘을 수 없다. AUC 와 BEDROC 도 같은 이유로
    부풀려진다. 가상 스크리닝 벤치마크의 통상 비율은 활성 2% 안팎이다.
    """
    bench = [r for r in lib if r["source"] == "benchmark_spikein"]
    others = [r for r in lib if r["source"] != "benchmark_spikein"]
    actives = [r for r in others if r["label"] == "known_active"]
    decoys = [r for r in others if r["label"] != "known_active"]

    room = max(0, cap - len(bench))
    # 활성도가 센 순(값이 작은 순). 값 없는 것은 뒤로.
    actives.sort(key=lambda r: r.get("activity_nM") or float("inf"))
    dropped = max(0, len(actives) - room)
    return bench + actives[:room] + decoys, dropped


def build(
    uniprot: str,
    benchmarks: list[str],
    n_decoys: int,
    max_actives: int,
    n_screen: int = 0,
) -> tuple[list[dict], dict]:
    print("  [1/4] 벤치마크 리간드 수집")
    bench = fetch_benchmark_ligands(benchmarks)

    print("  [2/4] ChEMBL 활성 데이터 수집")
    actives = fetch_chembl_actives(uniprot)
    print(f"    활성 {len(actives)}건 (10 uM 이하)")

    print("  [3/4] 화합물명 보강")
    annotate_names(actives)

    print("  [4/4] Decoy 수집")
    decoys = fetch_decoys(n_decoys)
    print(f"    decoy {len(decoys)}개")

    # 스크리닝 전용 집합 — 평가셋과 겹치지 않는 별도 화합물.
    #
    # 2026-10-23 이전에는 이 구분이 없어 decoy 전체가 채점의 음성이면서 동시에
    # 후보 추출원이었다. EPAS1·MET 양 타깃 모두 상위 20 후보가 전원 decoy 로
    # 나왔는데, 파이프라인이 제 성능을 채점할 때 음성으로 센 바로 그 화합물을
    # 후보로 내놓았다는 뜻이다. EF1% 가 천장에 닿으면 후보가 0개가 되므로
    # 후보는 지표의 오차 안에만 존재했다.
    #
    # decoy 를 쪼개 쓰는 길도 있었으나 평가셋이 작아지면서 활성 비율이 오르고
    # EF 천장이 내려간다(50% 분할 시 2.78%→5.41%, 천장 35.95→18.48). 기존
    # 결과와 비교가 끊기므로 평가셋은 그대로 두고 screen 을 위에 더한다.
    #
    # 기본값 0 — 측정에는 screen 이 필요 없다. 후보를 뽑을 때만 명시적으로
    # 켜고, 그만큼 도킹 비용이 더 든다.
    screen: list[dict] = []
    if n_screen:
        screen = fetch_decoys(n_screen, seed=SCREEN_SEED)
        for r in screen:
            r["label"] = "screen"
        print(f"    스크리닝 전용 {len(screen)}개 (도킹 비용 추가)")

    # 중복 제거 — 벤치마크 > 활성 > decoy 순으로 우선
    seen: dict[str, dict] = {}
    rejected: dict[str, int] = {}
    for rec in bench + actives + decoys + screen:
        cs = canonical(rec["smiles"])
        if cs is None:
            rejected["unparseable"] = rejected.get("unparseable", 0) + 1
            continue
        if cs in seen:
            continue
        if not rec.get("protected"):
            ok, reason = passes_filters(cs)
            if not ok:
                rejected[reason] = rejected.get(reason, 0) + 1
                continue
        rec["smiles"] = cs
        rec["rotb"] = Descriptors.NumRotatableBonds(Chem.MolFromSmiles(cs))
        seen[cs] = rec

    lib, n_dropped = cap_actives(list(seen.values()), max_actives)

    n_act = sum(1 for r in lib if r["label"] == "known_active")
    n_eval = n_act + sum(1 for r in lib if r["label"] == "decoy")
    ratio = n_act / n_eval if n_eval else 0
    stats = {
        "active_fraction": round(ratio, 4),
        "max_reachable_ef1": round(1 / ratio, 2) if ratio else None,
        "n_actives_dropped": n_dropped,
        "max_actives_cap": max_actives,
        "n_total": len(lib),
        "n_eval_set": n_eval,
        "n_benchmark": sum(1 for r in lib if r["source"] == "benchmark_spikein"),
        "n_known_active": sum(1 for r in lib if r["label"] == "known_active"),
        "n_decoy": sum(1 for r in lib if r["label"] == "decoy"),
        "n_screen": sum(1 for r in lib if r["label"] == "screen"),
        "benchmarks_requested": benchmarks,
        "benchmarks_found": [
            r["query_name"] for r in lib if r["source"] == "benchmark_spikein"
        ],
        "rejected_by": rejected,
        "activity_cutoff_nM": ACTIVITY_CUTOFF_NM,
    }
    return lib, stats


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--uniprot")
    ap.add_argument("--n-decoys", type=int, default=2500)
    ap.add_argument(
        "--n-screen",
        type=int,
        default=0,
        help="스크리닝 전용 화합물 수. 0이면 후보를 뽑지 않는다(측정 전용). "
        "켜면 그만큼 도킹 비용이 늘어난다.",
    )
    ap.add_argument(
        "--max-actives",
        type=int,
        default=60,
        help="활성 상한. 활성 비율이 높으면 EF 임계값이 도달 불가능해진다",
    )
    ap.add_argument("--targetdir", default=None)
    args = ap.parse_args()

    tdir = Path(args.targetdir) if args.targetdir else root / "targets" / args.gene
    tdir.mkdir(parents=True, exist_ok=True)

    uniprot, benchmarks = args.uniprot, []
    meta_path = tdir / "target_meta.json"
    if meta_path.exists():
        meta = json.loads(meta_path.read_text())
        uniprot = uniprot or meta.get("uniprot")
        benchmarks = meta.get("benchmark_ligands", [])
    if not benchmarks or not uniprot:
        import yaml

        cfg = yaml.safe_load((root / "config" / "targets.yaml").read_text())
        for t in cfg["targets"]:
            if t["gene"].upper() == args.gene.upper():
                uniprot = uniprot or t["uniprot"]
                benchmarks = benchmarks or t.get("benchmark_ligands", [])
    if not uniprot:
        print(f"[!] {args.gene}: UniProt accession을 찾을 수 없음", file=sys.stderr)
        return 1

    print(f"[*] {args.gene} ({uniprot}) 라이브러리 구성")
    lib, stats = build(
        uniprot, benchmarks, args.n_decoys, args.max_actives, n_screen=args.n_screen
    )

    (tdir / "library.jsonl").write_text(
        "\n".join(json.dumps(r, ensure_ascii=False) for r in lib)
    )
    (tdir / "library_stats.json").write_text(
        json.dumps(stats, indent=2, ensure_ascii=False)
    )

    print(
        f"[+] 총 {stats['n_total']} "
        f"(활성 {stats['n_known_active']}, decoy {stats['n_decoy']}, "
        f"벤치마크 {stats['n_benchmark']})"
    )
    print(f"    제외: {stats['rejected_by']}")
    print(
        f"    활성 비율 {100 * stats['active_fraction']:.1f}% "
        f"→ EF1% 이론적 최대 {stats['max_reachable_ef1']}"
    )
    if stats["n_actives_dropped"]:
        print(f"    활성 {stats['n_actives_dropped']}개는 상한 초과로 제외 (약한 것부터)")
    missing = set(benchmarks) - set(stats["benchmarks_found"])
    if missing:
        print(f"    [!] 벤치마크 미포함: {sorted(missing)} — 검증 기준점이 약해진다")
    if stats["n_known_active"] < 10:
        print("    [!] 활성 10개 미만 — 검증 게이트가 INSUFFICIENT로 판정될 것")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
