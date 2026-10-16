#!/usr/bin/env python3
"""
modules.small_molecule.ensemble — 앙상블 도킹 결과를 합친다.

왜 필요한가
  리지드 수용체 도킹은 리간드가 들어올 때 일어나는 측쇄 재배열을 모델링하지
  못한다. EPAS1 에서 belzutifan 이 상위 1% 를 놓친 것도 이 때문으로 보인다 —
  Ren et al. 2022 (Mol Pharmacol 102:240) 은 belzutifan 결합 시 M252 측쇄가
  포켓 밖으로 뒤집힌다고 보고했다. 한 구조에 고정하면 그 구조가 수용하는
  형태의 리간드만 잘 잡는다.

  앙상블 도킹은 여러 수용체 구조에 각각 도킹하고 리간드별 최고 스코어를 취해
  수용체 유연성을 근사한다. 표준 해법이다.

주의 — 앙상블은 편향을 없애지 않는다
  벤치마크와 닮은 리간드가 결합한 구조를 앙상블에 넣으면 그 구조가 최고
  스코어를 독점하고, 결국 단일 구조 결과와 같아진다. 그래서 이 모듈은 리간드마다
  '어느 수용체가 최고 스코어를 냈는지'를 기록한다. 한 구조가 대부분을 가져가면
  앙상블이 실질적으로 작동하지 않은 것이다.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def load_arm(tdir: Path) -> tuple[list[dict], dict]:
    rows: list[dict] = []
    allf = tdir / "docked_all.jsonl"
    if allf.exists():
        rows = [json.loads(l) for l in allf.read_text().splitlines() if l.strip()]
    else:
        for shard in sorted(tdir.glob("docked_shard*.jsonl")):
            rows += [json.loads(l) for l in shard.read_text().splitlines() if l.strip()]
    meta_f = tdir / "target_meta.json"
    meta = json.loads(meta_f.read_text()) if meta_f.exists() else {}
    return [r for r in rows if r.get("score") is not None], meta


def combine(arms: dict[str, list[dict]]) -> tuple[list[dict], dict]:
    """리간드별로 가장 좋은(낮은) 스코어를 취한다."""
    best: dict[str, dict] = {}
    for label, rows in arms.items():
        for r in rows:
            key = r.get("chembl_id") or r.get("smiles")
            if key is None:
                continue
            cur = best.get(key)
            if cur is None or r["score"] < cur["score"]:
                best[key] = {**r, "best_receptor": label}
            # 어느 수용체에서 얼마였는지 전부 남긴다
            best[key].setdefault("per_receptor", {})
            best[key]["per_receptor"][label] = r["score"]

    out = sorted(best.values(), key=lambda r: r["score"])
    winner = Counter(r["best_receptor"] for r in out)
    n = len(out) or 1
    stats = {
        "n_ligands": len(out),
        "n_arms": len(arms),
        "arm_sizes": {k: len(v) for k, v in arms.items()},
        "best_receptor_share": {k: round(100 * c / n, 1) for k, c in winner.most_common()},
    }
    top = winner.most_common(1)[0] if winner else (None, 0)
    stats["dominant_receptor"] = top[0]
    stats["dominant_share"] = round(100 * top[1] / n, 1)
    stats["effective"] = stats["dominant_share"] < 85.0
    stats["effectiveness_note"] = (
        f"{stats['dominant_receptor']} 가 {stats['dominant_share']}% 의 리간드에서 최고 스코어를 냈다. "
        + (
            "앙상블이 실질적으로 작동했다 — 리간드마다 다른 수용체가 선택됐다."
            if stats["effective"]
            else "한 구조가 거의 전부를 가져갔다. 앙상블이 단일 구조와 사실상 같다."
        )
    )
    return out, stats


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument(
        "--arm",
        action="append",
        required=True,
        metavar="LABEL=DIR",
        help="앙상블 구성원. 예: 6CZW=arms/treatment/EPAS1 (여러 번 지정)",
    )
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    arms, metas = {}, {}
    for spec in args.arm:
        if "=" not in spec:
            print(f"[!] --arm 형식 오류: {spec} (LABEL=DIR 이어야 함)")
            return 1
        label, d = spec.split("=", 1)
        rows, meta = load_arm(Path(d))
        if not rows:
            print(f"[!] {label}: 도킹 결과 없음 ({d})")
            return 1
        arms[label] = rows
        metas[label] = {
            "pdb_id": meta.get("pdb_id"),
            "resolution": meta.get("resolution"),
            "ligand": ((meta.get("receptor") or {}).get("cocrystal_ligand") or {}).get("resname"),
            "similarity": (meta.get("ligand_match") or {}).get("similarity"),
        }
        print(f"  {label}: {len(rows)} 리간드  ({metas[label]['pdb_id']}, "
              f"리간드 {metas[label]['ligand']}, 유사도 {metas[label]['similarity']})")

    combined, stats = combine(arms)
    stats["receptors"] = metas

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "docked_all.jsonl").write_text(
        "\n".join(json.dumps(r, ensure_ascii=False) for r in combined)
    )
    # 게이트가 벤치마크 목록을 읽을 수 있도록 메타를 하나 만들어 둔다
    base_meta = next(iter(arms))
    _, m0 = load_arm(Path(args.arm[0].split("=", 1)[1]))
    m0["pdb_id"] = "+".join(v["pdb_id"] or "?" for v in metas.values())
    m0["structure_selection"] = "ensemble"
    m0["ensemble"] = stats
    (outdir / "target_meta.json").write_text(json.dumps(m0, indent=2, ensure_ascii=False))
    (outdir / "ensemble_stats.json").write_text(json.dumps(stats, indent=2, ensure_ascii=False))

    print()
    print(f"[*] {args.gene} 앙상블 {stats['n_arms']}개 구조, 리간드 {stats['n_ligands']}개")
    print("    최고 스코어를 낸 수용체 분포:")
    for k, pct in stats["best_receptor_share"].items():
        print(f"      {k:10s} {pct:5.1f}%")
    print(f"    {stats['effectiveness_note']}")
    print(f"[+] {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
