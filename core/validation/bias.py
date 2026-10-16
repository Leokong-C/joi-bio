#!/usr/bin/env python3
"""
core.validation.bias — 구조 선택 편향의 크기를 측정한다.

왜 필요한가
  Chaput et al. (2016, J Cheminform 8:56) 은 도킹 프로그램 4종을 벤치마크한 뒤,
  Glide 의 좋은 성능이 "활성 화합물과 결정구조 리간드의 유사성" 으로 설명된다고
  결론냈다. 편향이 큰 타깃을 빼고 47개만 남기자 성공 타깃이 Glide 5개, Gold 4개,
  FlexX·Surflex 각 2개로 급감했다.

  우리 파이프라인은 벤치마크 화합물과 닮은 리간드가 결합한 구조를 고른다.
  성능이 좋게 나와도 그것이 실력인지 편향인지 구분할 수 없다. 순환논증이다.

어떻게 재는가
  같은 라이브러리를 두 수용체에 도킹한다.
    treatment  유사도로 고른 구조 (예: EPAS1 → 6CZW, 유사도 0.400)
    control    유사도가 가장 낮은 구조 (예: EPAS1 → 3H82, 유사도 0.071)
  두 결과의 차이가 편향의 크기다.

어떤 지표로 비교하는가
  AUC            활성이 decoy 보다 상위에 랭크될 확률. 활성:decoy 비율과 무관해
                 대조군의 decoy 수가 적어도 비교가 성립한다.
  벤치마크 백분위  승인약이 상위 몇 % 에 들어왔는가. 라이브러리 크기로 정규화된다.
  EF 는 비율 의존이라 두 팔의 라이브러리 크기가 다르면 쓰지 않는다.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from core.validation.gate import benchmark_recovery, roc_auc

ROOT = Path(__file__).resolve().parents[2]


def load_rows(tdir: Path) -> list[dict]:
    rows: list[dict] = []
    allf = tdir / "docked_all.jsonl"
    if allf.exists():
        rows = [json.loads(l) for l in allf.read_text().splitlines() if l.strip()]
    else:
        for shard in sorted(tdir.glob("docked_shard*.jsonl")):
            rows += [json.loads(l) for l in shard.read_text().splitlines() if l.strip()]
    return [r for r in rows if r.get("score") is not None]


def arm_metrics(rows: list[dict], benchmarks: list[str]) -> dict:
    rows = sorted(rows, key=lambda r: r["score"])
    labels = [1 if r.get("label") == "known_active" else 0 for r in rows]
    br = benchmark_recovery(rows, benchmarks)
    return {
        "n": len(rows),
        "n_active": sum(labels),
        "auc": roc_auc(labels),
        "best_score": rows[0]["score"] if rows else None,
        "median_active": _median([r["score"] for r, l in zip(rows, labels) if l]),
        "median_decoy": _median([r["score"] for r, l in zip(rows, labels) if not l]),
        "benchmarks": {f["name"]: f["percentile"] for f in br.get("checked", [])},
    }


def _median(xs: list[float]) -> float | None:
    if not xs:
        return None
    xs = sorted(xs)
    n = len(xs)
    return round(xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2, 3)


def compare(treat: dict, ctrl: dict) -> dict:
    """두 팔을 공통 화합물로 맞춰 비교한다."""
    out = {"treatment": treat["metrics"], "control": ctrl["metrics"]}

    t_auc, c_auc = treat["metrics"]["auc"], ctrl["metrics"]["auc"]
    if t_auc is not None and c_auc is not None:
        delta = round(t_auc - c_auc, 4)
        out["auc_delta"] = delta
        # 해석 기준은 보수적으로 잡는다. AUC 0.05 는 랭킹이 눈에 띄게 달라지는 폭이다.
        if delta >= 0.15:
            verdict = "BIAS_LIKELY"
            reason = (
                f"유사도 구조가 AUC 를 {delta:.3f} 끌어올렸다. 성능의 상당 부분이 "
                "구조 선택에서 나왔을 수 있으므로, 결과를 보고할 때 두 값을 함께 낸다."
            )
        elif delta >= 0.05:
            verdict = "BIAS_MODERATE"
            reason = f"AUC 차이 {delta:.3f}. 편향이 있으나 지배적이지는 않다."
        elif delta >= -0.05:
            verdict = "BIAS_LOW"
            reason = (
                f"AUC 차이 {delta:.3f}. 구조 선택이 성능을 거의 바꾸지 않았다. "
                "농축도가 화합물 자체의 상보성에서 나왔다는 뜻이다."
            )
        else:
            verdict = "CONTROL_BETTER"
            reason = (
                f"대조군이 오히려 {-delta:.3f} 높다. 유사도 기반 구조 선택이 "
                "역효과를 냈다. 선택 로직을 재검토할 것."
            )
        out["verdict"] = verdict
        out["verdict_reason"] = reason
    else:
        out["verdict"] = "INCOMPARABLE"
        out["verdict_reason"] = "한쪽 팔의 AUC 를 계산할 수 없다 (활성 또는 decoy 부재)."

    # 벤치마크 백분위 대조 — 승인약이 어느 쪽에서 더 잘 잡히는가
    rows = []
    for name in sorted(set(treat["metrics"]["benchmarks"]) | set(ctrl["metrics"]["benchmarks"])):
        t = treat["metrics"]["benchmarks"].get(name)
        c = ctrl["metrics"]["benchmarks"].get(name)
        rows.append(
            {
                "name": name,
                "treatment_percentile": t,
                "control_percentile": c,
                "improvement": round(c - t, 3) if (t is not None and c is not None) else None,
            }
        )
    out["benchmark_comparison"] = rows
    out["structures"] = {
        "treatment": treat["structure"],
        "control": ctrl["structure"],
    }
    return out


def arm(tdir: Path, benchmarks: list[str] | None = None) -> dict:
    meta_f = tdir / "target_meta.json"
    meta = json.loads(meta_f.read_text()) if meta_f.exists() else {}
    bench = benchmarks if benchmarks is not None else meta.get("benchmark_ligands", [])
    rows = load_rows(tdir)
    return {
        "metrics": arm_metrics(rows, bench),
        "structure": {
            "pdb_id": meta.get("pdb_id"),
            "resolution": meta.get("resolution"),
            "similarity": (meta.get("ligand_match") or {}).get("similarity"),
            "selection": meta.get("structure_selection"),
            "ligand": ((meta.get("receptor") or {}).get("cocrystal_ligand") or {}).get("resname"),
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--treatment", required=True, help="유사도 선택 구조의 결과 디렉터리")
    ap.add_argument("--control", required=True, help="강제 지정 구조의 결과 디렉터리")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    t = arm(Path(args.treatment))
    c = arm(Path(args.control), benchmarks=None)
    result = {"gene": args.gene, **compare(t, c)}

    out = Path(args.out) if args.out else ROOT / "results" / f"bias_{args.gene}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False))

    ts, cs = result["structures"]["treatment"], result["structures"]["control"]
    print(f"[*] {args.gene} 구조 선택 편향 측정")
    print(
        f"    treatment {ts['pdb_id']} ({ts['ligand']}, 유사도 {ts['similarity']}, "
        f"{ts['resolution']} Å)"
    )
    print(
        f"    control   {cs['pdb_id']} ({cs['ligand']}, 유사도 {cs['similarity']}, "
        f"{cs['resolution']} Å)"
    )
    print()
    print(f"    AUC  treatment={result['treatment']['auc']}  control={result['control']['auc']}")
    print(f"    차이 {result.get('auc_delta')}")
    print()
    for b in result["benchmark_comparison"]:
        print(
            f"    {b['name']:16s} treatment {b['treatment_percentile']}%  "
            f"control {b['control_percentile']}%  개선 {b['improvement']}"
        )
    print()
    print(f"    판정: {result['verdict']}")
    print(f"    {result['verdict_reason']}")
    print(f"[+] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
