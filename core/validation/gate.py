#!/usr/bin/env python3
"""
core.validation.gate — 파이프라인 자체를 검증하는 게이트.

Joi_bio 신뢰성의 근거. 알려진 활성 화합물(known actives)을 라이브러리에 섞어
넣고, 도킹 스코어가 이들을 상위로 올리는지 측정한다.

게이트를 통과하지 못한 타깃의 화합물은 어떤 보고서에도 '히트'로 실리지 않는다.
검증 없이 나온 후보 리스트는 파는 순간 신뢰를 잃기 때문이다.

지표
  EF1%   상위 1%에서 활성이 무작위 대비 몇 배로 농축됐는가
  AUC    활성이 decoy보다 상위에 랭크될 확률
  BEDROC 상위 가중 농축도 (alpha=20 → 상위 8% 구간에 가중)
  benchmark_recovery  승인약/임상 화합물이 실제로 상위 1%에 들어왔는가

판정
  INSUFFICIENT  활성 수가 10개 미만 → 통계적으로 무의미, 지표를 믿지 않는다
  INVALID       활성 비율이 높아 EF 임계값이 도달 불가능 → 라이브러리 재구성 필요
  PASS          유효한 임계값 중 2개 이상 충족
  FAIL          그 외
"""
from __future__ import annotations

import json
import math
from pathlib import Path

THRESHOLDS = {"ef1": 5.0, "auc": 0.70, "bedroc": 0.30}
# 임계값은 활성 비율 2~3% 인 통상적 벤치마크를 전제로 한다.
# 활성 비율이 높으면 EF 상한이 임계값 아래로 내려가고 AUC/BEDROC 도 부풀려진다.
MIN_ACTIVES = 10
BEDROC_ALPHA = 20.0


def enrichment_factor(labels: list[int], frac: float) -> float | None:
    """labels는 스코어 좋은 순으로 정렬된 0/1 배열."""
    n, n_act = len(labels), sum(labels)
    if n == 0 or n_act == 0:
        return None
    top_n = max(1, int(round(n * frac)))
    hits = sum(labels[:top_n])
    return round((hits / top_n) / (n_act / n), 3)


def roc_auc(labels: list[int]) -> float | None:
    """랭킹 기반 AUC (Mann-Whitney U). labels[0]이 최상위."""
    n_act = sum(labels)
    n_dec = len(labels) - n_act
    if n_act == 0 or n_dec == 0:
        return None
    rank_sum = sum(i + 1 for i, l in enumerate(labels) if l == 1)
    u = rank_sum - n_act * (n_act + 1) / 2  # decoy가 active보다 앞선 쌍의 수
    return round(1.0 - u / (n_act * n_dec), 4)


def bedroc(labels: list[int], alpha: float = BEDROC_ALPHA) -> float | None:
    """Truchon & Bayly (2007) BEDROC. labels[0]이 최상위."""
    n = len(labels)
    n_act = sum(labels)
    if n == 0 or n_act == 0 or n_act == n:
        return None
    ra = n_act / n

    s = sum(math.exp(-alpha * (i + 1) / n) for i, l in enumerate(labels) if l == 1)
    denom = math.exp(alpha / n) - 1.0
    if denom <= 0:
        return None
    random_sum = ra * (1.0 - math.exp(-alpha)) / denom
    if random_sum <= 0:
        return None
    rie = s / random_sum

    try:
        factor = (
            ra
            * math.sinh(alpha / 2.0)
            / (math.cosh(alpha / 2.0) - math.cosh(alpha / 2.0 - alpha * ra))
        )
        offset = 1.0 / (1.0 - math.exp(alpha * (1.0 - ra)))
    except (OverflowError, ZeroDivisionError):
        return None
    return round(rie * factor + offset, 4)


def benchmark_recovery(rows: list[dict], names: list[str], frac: float = 0.01) -> dict:
    """승인약·임상 화합물이 상위 frac 안에 들어왔는지 이름으로 확인."""
    if not names:
        return {"checked": [], "note": "벤치마크 리간드 미지정"}
    n = len(rows)
    top_n = max(1, int(round(n * frac)))
    wanted = {x.lower() for x in names}

    found = []
    for idx, r in enumerate(rows):
        label = " ".join(
            str(r.get(k, "")) for k in ("pref_name", "chembl_id", "synonyms")
        ).lower()
        for w in wanted:
            if w in label:
                found.append(
                    {
                        "name": w,
                        "rank": idx + 1,
                        "percentile": round(100.0 * (idx + 1) / n, 3),
                        "score": r.get("score"),
                        "chembl_id": r.get("chembl_id"),
                        "in_top_pct": idx < top_n,
                    }
                )
                break

    # 같은 이름이 여러 번 나오면 가장 좋은 순위만 남긴다
    best: dict[str, dict] = {}
    for f in found:
        if f["name"] not in best or f["rank"] < best[f["name"]]["rank"]:
            best[f["name"]] = f

    checked = sorted(best.values(), key=lambda f: f["rank"])
    return {
        "checked": checked,
        "n_requested": len(names),
        "n_found_in_library": len(checked),
        "n_in_top_pct": sum(1 for f in checked if f["in_top_pct"]),
        "top_pct": frac * 100,
        "missing": sorted(wanted - set(best)),
    }


def evaluate(rows: list[dict], benchmark_ligands: list[str] | None = None) -> dict:
    """도킹 결과(스코어 좋은 순 정렬)를 받아 게이트 판정을 내린다."""
    rows = [r for r in rows if r.get("score") is not None]
    rows.sort(key=lambda r: r["score"])  # Vina는 낮을수록 좋다
    labels = [1 if r.get("label") == "known_active" else 0 for r in rows]

    n_act = sum(labels)
    metrics = {
        "n_compounds": len(rows),
        "n_known_active": n_act,
        "n_decoy": len(rows) - n_act,
        "ef1": enrichment_factor(labels, 0.01),
        "ef5": enrichment_factor(labels, 0.05),
        "ef10": enrichment_factor(labels, 0.10),
        "auc": roc_auc(labels),
        "bedroc": bedroc(labels),
        "best_score": rows[0]["score"] if rows else None,
        "thresholds": THRESHOLDS,
        "min_actives_required": MIN_ACTIVES,
    }
    metrics["benchmark_recovery"] = benchmark_recovery(rows, benchmark_ligands or [])

    # EF 는 비율로 나누므로 상한이 있다. 활성 비율이 높으면 임계값이
    # 어떤 성능으로도 도달 불가능해지는데, 이를 모르고 FAIL 로 읽으면
    # 파이프라인 탓을 하게 된다. 실제로 EPAS1 첫 실행에서 그랬다.
    ratio = n_act / len(rows) if rows else 0
    max_ef1 = round(1 / ratio, 2) if ratio else None
    metrics["active_fraction"] = round(ratio, 4)
    metrics["max_reachable_ef1"] = max_ef1
    ef1_unreachable = max_ef1 is not None and THRESHOLDS["ef1"] > max_ef1
    metrics["ef1_threshold_unreachable"] = ef1_unreachable

    checks = {
        k: (metrics[k] is not None and metrics[k] >= v) for k, v in THRESHOLDS.items()
    }
    if ef1_unreachable:
        checks["ef1"] = None  # 판정에서 제외 — 통과도 미달도 아니다
    metrics["checks"] = checks

    passed_checks = sum(1 for v in checks.values() if v is True)
    usable_checks = sum(1 for v in checks.values() if v is not None)

    if n_act < MIN_ACTIVES:
        metrics["verdict"] = "INSUFFICIENT"
        metrics["verdict_reason"] = (
            f"활성 화합물 {n_act}개 (최소 {MIN_ACTIVES}개 필요). "
            "지표가 계산되더라도 통계적으로 무의미하다."
        )
    elif ef1_unreachable and usable_checks < 3:
        metrics["verdict"] = "INVALID"
        metrics["verdict_reason"] = (
            f"활성 비율 {100 * ratio:.1f}% — EF1% 이론적 최대가 {max_ef1}이라 "
            f"임계값 {THRESHOLDS['ef1']}에 도달할 수 없다. 남은 지표(AUC, BEDROC)도 "
            "활성 비율이 높으면 함께 부풀려지므로 믿을 수 없다. "
            "라이브러리를 활성 2~3% 비율로 다시 구성할 것."
        )
    elif passed_checks >= 2:
        metrics["verdict"] = "PASS"
        passed = [k for k, v in checks.items() if v is True]
        metrics["verdict_reason"] = f"{passed_checks}/{usable_checks} 통과: {', '.join(passed)}"
    else:
        metrics["verdict"] = "FAIL"
        failed = [k for k, v in checks.items() if v is False]
        metrics["verdict_reason"] = (
            f"미달: {', '.join(failed)}. "
            "포켓 정의, 구조 선택, 또는 이 타깃 클래스에 대한 "
            "리지드 도킹의 한계를 의심할 것."
        )

    # 벤치마크 회수는 별도 경고. 지표가 통과해도 승인약을 못 찾으면
    # 그 파이프라인은 신뢰할 수 없다.
    br = metrics["benchmark_recovery"]
    if br.get("n_found_in_library") and not br.get("n_in_top_pct"):
        metrics["benchmark_warning"] = (
            f"벤치마크 {br['n_found_in_library']}개 중 상위 {br['top_pct']}%에 "
            "진입한 것이 없다. 지표가 통과하더라도 이 결과는 신뢰할 수 없다."
        )
    return metrics


def main() -> int:
    import argparse

    root = Path(__file__).resolve().parents[2]
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--targetdir", default=None)
    args = ap.parse_args()

    tdir = Path(args.targetdir) if args.targetdir else root / "targets" / args.gene

    # 샤드 파일이 기본이지만, 앙상블처럼 합쳐진 결과는 docked_all.jsonl 로 온다.
    # 이 폴백이 없어 앙상블 검증이 "도킹 결과 없음" 으로 멈춘 적이 있다.
    rows: list[dict] = []
    for shard in sorted(tdir.glob("docked_shard*.jsonl")):
        rows += [json.loads(l) for l in shard.read_text().splitlines() if l.strip()]
    if not rows:
        allf = tdir / "docked_all.jsonl"
        if allf.exists():
            rows = [json.loads(l) for l in allf.read_text().splitlines() if l.strip()]
    if not rows:
        print(f"[!] {args.gene}: 도킹 결과 없음")
        return 1

    meta_path = tdir / "target_meta.json"
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    bench = meta.get("benchmark_ligands", [])

    result = evaluate(rows, bench)
    result["gene"] = args.gene
    result["pocket_confidence"] = meta.get("pocket", {}).get("confidence")

    scored = sorted(
        (r for r in rows if r.get("score") is not None), key=lambda r: r["score"]
    )
    (tdir / "docked_all.jsonl").write_text(
        "\n".join(json.dumps(r, ensure_ascii=False) for r in scored)
    )
    (tdir / "validation.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False)
    )

    print(f"[*] {args.gene} 검증")
    print(f"    화합물 {result['n_compounds']} (활성 {result['n_known_active']})")
    for k in ("ef1", "ef5", "auc", "bedroc"):
        v = result[k]
        print(f"    {k:8s} = {v}" if v is not None else f"    {k:8s} = n/a")

    br = result["benchmark_recovery"]
    if br.get("checked"):
        print(f"    벤치마크 리간드 ({br['n_in_top_pct']}/{br['n_found_in_library']} 상위 1%):")
        for f in br["checked"]:
            mark = "OK " if f["in_top_pct"] else "   "
            print(
                f"      {mark}{f['name']:16s} rank {f['rank']} "
                f"({f['percentile']}%) score {f['score']}"
            )
    if br.get("missing"):
        print(f"    라이브러리에 없음: {br['missing']}")

    print(f"    판정: {result['verdict']} — {result['verdict_reason']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
