#!/usr/bin/env python3
"""
core.validation.perturb — 도구 출력 교란 감사.

파이프라인에 틀린 입력을 일부러 먹이고 **눈치채는지** 본다.

왜 필요한가
  FDA·EMA Good AI Practice 8번 원칙(2026-01-14)은 "인간 상호작용을 포함한
  시스템 전체"의 성능 평가를 요구한다. 지금까지 우리가 잰 것은 도킹 스코어의
  분별력뿐이고, 입력이 오염됐을 때 파이프라인이 어떻게 되는지는 한 번도
  재지 않았다.

  설계는 Pseudo-Lab/abstain-dti 의 perturb/ 하네스에서 빌렸다
  (단위 뒤바꿈, 철회된 레코드, SMILES 교체, 타겟 ID 교체).

무엇을 재는가
  각 교란을 주입한 뒤 게이트 판정과 지표가 어떻게 변하는지 본다.
  판정이 바뀌면 탐지(detected), 그대로면 미탐지(missed).

  미탐지가 많이 나올 것으로 예상한다. 그게 결과다 — 지금 파이프라인에는
  입력 무결성을 보는 장치가 하나도 없다.

계산 비용
  0. 기존 도킹 결과를 재채점할 뿐 새로 도킹하지 않는다.
"""
from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

SEED = 1023


def _rows(path: Path) -> list[dict]:
    rows: list[dict] = []
    for shard in sorted(path.glob("docked_shard*.jsonl")):
        rows += [json.loads(l) for l in shard.read_text().splitlines() if l.strip()]
    if not rows:
        allf = path / "docked_all.jsonl"
        if allf.exists():
            rows = [json.loads(l) for l in allf.read_text().splitlines() if l.strip()]
    return [r for r in rows if r.get("score") is not None]


# ---------- 교란 ----------


def unit_flip(rows: list[dict], frac: float = 0.2) -> tuple[list[dict], str]:
    """활성 판정의 단위가 nM↔µM 로 뒤바뀐 상황.

    1000배 차이이므로, 10 µM 문턱 아래였던 화합물 일부가 사실은 위였다는 뜻.
    라벨이 known_active 에서 decoy 로 뒤집힌다.
    """
    out = [dict(r) for r in rows]
    rnd = random.Random(SEED)
    acts = [r for r in out if r.get("label") == "known_active"]
    n = int(len(acts) * frac)
    for r in rnd.sample(acts, n):
        r["label"] = "decoy"
        r["perturbed"] = "unit_flip"
    return out, f"활성 {n}/{len(acts)}개의 단위가 뒤바뀌어 decoy 로 재분류"


def swap_smiles(rows: list[dict], frac: float = 0.2) -> tuple[list[dict], str]:
    """레코드의 SMILES 가 다른 화합물 것으로 바뀐 상황.

    스코어는 원래 화합물의 것이고 구조만 남의 것이다. 파이프라인이 스코어와
    구조의 대응을 검사하지 않으면 그대로 지나간다.
    """
    out = [dict(r) for r in rows]
    rnd = random.Random(SEED + 1)
    idx = rnd.sample(range(len(out)), int(len(out) * frac))
    smis = [out[i].get("smiles") for i in idx]
    rnd.shuffle(smis)
    for i, s in zip(idx, smis):
        out[i]["smiles"] = s
        out[i]["perturbed"] = "swap_smiles"
    return out, f"{len(idx)}개 레코드의 SMILES 를 서로 섞음 (스코어는 그대로)"


def swap_target(rows_a: list[dict], rows_b: list[dict]) -> tuple[list[dict], str]:
    """타깃을 착각한 상황 — A 의 라벨에 B 의 스코어를 붙인다.

    EPAS1 라이브러리를 MET 구조에 도킹한 결과를 EPAS1 결과로 보고하는 것과 같다.
    """
    kb = {r.get("chembl_id"): r for r in rows_b}
    out = []
    n = 0
    for r in rows_a:
        c = dict(r)
        other = kb.get(r.get("chembl_id"))
        if other is not None:
            c["score"] = other["score"]
            c["perturbed"] = "swap_target"
            n += 1
        out.append(c)
    return out, f"{n}개 화합물의 스코어를 다른 타깃 실행 것으로 교체"


def stale_records(rows: list[dict], frac: float = 0.1) -> tuple[list[dict], str]:
    """철회된 ChEMBL 레코드가 섞인 상황 — 활성 라벨이 근거를 잃었다."""
    out = [dict(r) for r in rows]
    rnd = random.Random(SEED + 2)
    acts = [r for r in out if r.get("label") == "known_active"]
    n = int(len(acts) * frac)
    for r in rnd.sample(acts, n):
        r["withdrawn"] = True
        r["perturbed"] = "stale_records"
    return out, f"활성 {n}/{len(acts)}개를 철회된 레코드로 표시 (라벨은 그대로)"


def main() -> int:
    from core.validation.gate import evaluate

    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--targetdir", required=True)
    ap.add_argument("--controldir", default=None)
    ap.add_argument("--otherdir", default=None, help="swap_target 용 다른 타깃 결과")
    ap.add_argument("--n-boot", type=int, default=500)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    rows = _rows(Path(args.targetdir))
    if not rows:
        print("[!] 도킹 결과 없음")
        return 1
    ctrl = _rows(Path(args.controldir)) if args.controldir else None
    other = _rows(Path(args.otherdir)) if args.otherdir else None

    meta_p = Path(args.targetdir) / "target_meta.json"
    bench = (
        json.loads(meta_p.read_text()).get("benchmark_ligands", [])
        if meta_p.exists()
        else []
    )

    base = evaluate(rows, bench, control_rows=ctrl, n_boot=args.n_boot)
    print(f"[*] {args.gene} 기준선: {base['verdict']} "
          f"(EF1 {base['ef1']}, AUC {base['auc']}, BEDROC {base['bedroc']})")

    cases = [
        ("unit_flip", *unit_flip(rows)),
        ("swap_smiles", *swap_smiles(rows)),
        ("stale_records", *stale_records(rows)),
    ]
    if other:
        cases.append(("swap_target", *swap_target(rows, other)))
    else:
        print("    swap_target 생략 — --otherdir 없음")

    results = []
    for name, perturbed, desc in cases:
        m = evaluate(perturbed, bench, control_rows=ctrl, n_boot=args.n_boot)
        changed = m["verdict"] != base["verdict"]
        results.append(
            {
                "case": name,
                "description": desc,
                "verdict": m["verdict"],
                "baseline_verdict": base["verdict"],
                "detected": changed,
                "ef1": m["ef1"],
                "auc": m["auc"],
                "bedroc": m["bedroc"],
                "delta_ef1": (
                    round(m["ef1"] - base["ef1"], 3)
                    if m["ef1"] is not None and base["ef1"] is not None
                    else None
                ),
            }
        )
        mark = "탐지" if changed else "미탐지"
        print(f"    {name:16s} {mark}  판정 {m['verdict']:11s} "
              f"EF1 {m['ef1']} (Δ{results[-1]['delta_ef1']})")
        print(f"                     {desc}")

    n_det = sum(1 for r in results if r["detected"])
    out = {
        "gene": args.gene,
        "baseline": {
            "verdict": base["verdict"],
            "ef1": base["ef1"],
            "auc": base["auc"],
            "bedroc": base["bedroc"],
        },
        "cases": results,
        "n_cases": len(results),
        "n_detected": n_det,
        "detection_rate": round(n_det / len(results), 3) if results else None,
        "note": (
            "미탐지는 파이프라인에 입력 무결성 검사가 없다는 뜻이다. "
            "스코어 분별력과 입력 신뢰성은 다른 문제다."
        ),
    }
    print(f"\n[=] 탐지율 {n_det}/{len(results)}")

    dest = Path(args.out) if args.out else Path(args.targetdir) / "perturb_audit.json"
    dest.write_text(json.dumps(out, ensure_ascii=False, indent=2))
    print(f"    → {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
