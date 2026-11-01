#!/usr/bin/env python3
"""
scripts/freeze_library.py — 미리 만들어 둔 라이브러리를 도킹 입력으로 고정한다.

왜 필요한가 (2026-11-01, A안)
  modules.small_molecule.library 는 decoy 를 ChEMBL 에서 생으로 질의해
  만든다. PubChem SMARTS 풀로 가는 경로가 없다. 그래서 워크플로를 그냥
  돌리면 "불소만 세도 AUC 0.911" 이 나오던 그 라이브러리가 다시 만들어진다.

  그리고 더 중요한 이유: 두 팔이 각자 라이브러리를 다시 만들면 같은
  라이브러리라는 보장이 없다. 편향 측정의 전제는 "구조만 다르고 나머지는
  동일" 이다. 전제가 깨지면 측정값이 의미를 잃는다. 그래서 한 파일을
  고정해 양쪽에 넣고, sha256 을 찍어 같았음을 증명한다.

하는 일 세 가지
  1. 묵은 도킹 결과 필드를 지운다. 고정 파일의 활성 60개에는 10월 실행의
     score/pose_spread/ligand_efficiency 가 남아 있다. 다른 수용체에서 나온
     값이므로 입력에 남겨 두면 안 된다.
  2. rotb 를 전부 다시 센다. dock.py 는 회전결합 수로 정렬한 뒤 라운드로빈으로
     샤드를 나눠 부하를 고르게 만든다. 고정 파일의 decoy 2,100개에는 rotb 가
     없어서(활성 60개만 있다) 전부 0 으로 읽히고, 무거운 리간드가 한 샤드에
     몰린다. EPAS1 1차 실행에서 샤드 3개가 350분 타임아웃으로 취소된 바로
     그 경로다.
  3. 게이트와 같은 공식으로 통계를 다시 적는다.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from core.validation.gate import perfect_ef  # noqa: E402

RDLogger.DisableLog("rdApp.*")

# 다른 수용체에서 나온 값이다. 입력에 남기지 않는다.
STALE = ("score", "pose_spread", "ligand_efficiency", "error")


def normalize(rows: list[dict]) -> tuple[list[dict], dict]:
    out, notes = [], {"rotb_recomputed": 0, "stale_stripped": 0, "unparseable": 0}
    for r in rows:
        rec = {k: v for k, v in r.items() if k not in STALE}
        if len(rec) != len(r):
            notes["stale_stripped"] += 1
        mol = Chem.MolFromSmiles(rec.get("smiles", ""))
        if mol is None:
            notes["unparseable"] += 1
            continue
        rotb = Descriptors.NumRotatableBonds(mol)
        if rec.get("rotb") != rotb:
            notes["rotb_recomputed"] += 1
        rec["rotb"] = rotb
        out.append(rec)
    return out, notes


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--library", required=True, help="고정할 라이브러리 jsonl")
    ap.add_argument("--gene", required=True)
    ap.add_argument("--targetdir", default=None)
    ap.add_argument("--label", default="", help="기록용 이름")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[1]
    src = Path(args.library)
    if not src.exists():
        print(f"[!] 라이브러리를 찾을 수 없다: {src}", file=sys.stderr)
        return 1

    tdir = Path(args.targetdir) if args.targetdir else root / "targets" / args.gene
    tdir.mkdir(parents=True, exist_ok=True)

    rows = [json.loads(l) for l in src.read_text().splitlines() if l.strip()]
    print(f"[*] {src} 에서 {len(rows):,}행")

    lib, notes = normalize(rows)
    print(f"    묵은 도킹 필드 제거 {notes['stale_stripped']}행")
    print(f"    rotb 재계산 {notes['rotb_recomputed']}행")
    if notes["unparseable"]:
        print(f"    [!] SMILES 파싱 실패로 제외 {notes['unparseable']}행")

    n_act = sum(1 for r in lib if r["label"] == "known_active")
    n_dec = sum(1 for r in lib if r["label"] == "decoy")
    n_screen = sum(1 for r in lib if r["label"] == "screen")
    n_eval = n_act + n_dec
    if n_act < 10:
        print(f"[!] 활성 {n_act}개 — 게이트 최소 10개 미달", file=sys.stderr)
        return 1

    # 결정적 직렬화. 두 팔이 같은 파일을 썼음을 sha256 으로 증명한다.
    body = "\n".join(json.dumps(r, ensure_ascii=False, sort_keys=True) for r in lib)
    (tdir / "library.jsonl").write_text(body)
    digest = hashlib.sha256(body.encode()).hexdigest()

    rotbs = sorted((r["rotb"] for r in lib), reverse=True)
    stats = {
        "frozen_from": str(src),
        "frozen_label": args.label or src.name,
        "sha256": digest,
        "n_total": len(lib),
        "n_eval_set": n_eval,
        "n_known_active": n_act,
        "n_decoy": n_dec,
        "n_screen": n_screen,
        "active_fraction": round(n_act / n_eval, 4) if n_eval else 0,
        "max_reachable_ef1": perfect_ef(n_eval, n_act, 0.01),
        "decoy_mode": "pubchem_smarts_property_matched",
        "rotb_max": rotbs[0] if rotbs else None,
        "rotb_median": rotbs[len(rotbs) // 2] if rotbs else None,
        "normalize_notes": notes,
    }
    (tdir / "library_stats.json").write_text(json.dumps(stats, indent=2, ensure_ascii=False))

    print(f"[+] {tdir / 'library.jsonl'}")
    print(f"    총 {len(lib):,} (활성 {n_act}, decoy {n_dec}, screen {n_screen})")
    print(f"    활성비율 {stats['active_fraction']:.4f}  EF1 천장 {stats['max_reachable_ef1']}")
    print(f"    rotb 최대 {stats['rotb_max']} 중위 {stats['rotb_median']}")
    print(f"    sha256 {digest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
