#!/usr/bin/env python3
"""
core.validation.stats — 지표에 오차범위를 붙인다.

왜 필요한가
  우리는 EF1% 를 32.682 처럼 소수점 셋째 자리까지 보고해 왔다. 그런데 이 값은
  2,157개 중 **상위 21개** 에만 의존한다. 활성 하나가 21위에서 22위로 밀리면
  EF1 은 32.68 에서 30.95 로 떨어진다. 소수점 셋째 자리는 의미가 없고, 더
  중요하게는 treatment 32.68 대 control 8.17 의 4배 차이가 진짜인지 추출 운인지
  지금까지 확인한 적이 없다.

무엇을 하는가
  부트스트랩  화합물을 복원추출해 지표를 다시 계산한다. 2,000회 반복해
              백분위수로 95% 신뢰구간을 얻는다.
  순열검정    라벨을 무작위로 섞어 지표를 다시 계산한다. 관측값이 섞은 분포의
              어디에 있는지가 p 값이다. '우연히 이만큼 나올 확률' 을 직접 센다.
  짝지은 비교 두 팔이 같은 라이브러리를 쓰므로 비교는 짝지어져 있다. 같은
              복원추출 인덱스를 양쪽에 적용해 차이의 신뢰구간을 구한다.
              각 팔을 따로 부트스트랩해 빼면 상관을 무시해 구간이 과대해진다.

한계
  부트스트랩은 우리가 가진 라이브러리 안에서의 불확실성만 잰다. 다른 decoy 를
  뽑았거나 활성 60개를 다르게 골랐을 때의 변동은 측정하지 않는다. 그것은
  라이브러리를 다시 만들어야 알 수 있다.
"""
from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

from core.validation.gate import bedroc, enrichment_factor, roc_auc

ROOT = Path(__file__).resolve().parents[2]
METRICS = ("ef1", "ef5", "auc", "bedroc")


def compute(labels: list[int]) -> dict[str, float | None]:
    return {
        "ef1": enrichment_factor(labels, 0.01),
        "ef5": enrichment_factor(labels, 0.05),
        "auc": roc_auc(labels),
        "bedroc": bedroc(labels),
    }


def _labels_from(rows: list[dict]) -> list[int]:
    """점수 오름차순으로 정렬한 뒤의 라벨 수열."""
    s = sorted(rows, key=lambda r: r["score"])
    return [1 if r.get("label") == "known_active" else 0 for r in s]


def bootstrap(rows: list[dict], n_boot: int, seed: int = 42) -> dict:
    """화합물을 복원추출해 지표 분포를 얻는다."""
    rnd = random.Random(seed)
    n = len(rows)
    dist: dict[str, list[float]] = {m: [] for m in METRICS}
    for _ in range(n_boot):
        sample = [rows[rnd.randrange(n)] for _ in range(n)]
        v = compute(_labels_from(sample))
        for m in METRICS:
            if v[m] is not None:
                dist[m].append(v[m])
    return dist


def ci(values: list[float], lo: float = 2.5, hi: float = 97.5) -> dict:
    if not values:
        return {"lo": None, "hi": None, "n": 0}
    v = sorted(values)
    k = lambda p: v[min(len(v) - 1, max(0, int(round(p / 100 * (len(v) - 1)))))]
    return {"lo": round(k(lo), 3), "hi": round(k(hi), 3), "n": len(v)}


def permutation(rows: list[dict], n_perm: int, seed: int = 42) -> dict:
    """라벨을 섞어 '우연히 이만큼' 의 분포를 만든다."""
    rnd = random.Random(seed)
    obs = compute(_labels_from(rows))
    labels = [1 if r.get("label") == "known_active" else 0 for r in rows]
    counts = {m: 0 for m in METRICS}
    for _ in range(n_perm):
        rnd.shuffle(labels)
        v = compute(labels)  # 섞인 라벨은 이미 '정렬 후 순서' 로 본다
        for m in METRICS:
            if v[m] is not None and obs[m] is not None and v[m] >= obs[m]:
                counts[m] += 1
    # 가산 보정 — 관측값 자신을 분포에 포함시킨다
    return {
        m: {
            "observed": obs[m],
            "p_value": round((counts[m] + 1) / (n_perm + 1), 5),
            "n_perm": n_perm,
        }
        for m in METRICS
    }


def paired_diff(a: list[dict], b: list[dict], n_boot: int, seed: int = 42) -> dict:
    """같은 복원추출 인덱스를 두 팔에 적용해 차이의 신뢰구간을 구한다.

    두 팔은 같은 라이브러리를 쓰므로 짝지어져 있다. 따로 부트스트랩해 빼면
    상관을 무시하게 되어 구간이 쓸데없이 넓어진다.
    """
    ka = {r.get("chembl_id"): r for r in a}
    kb = {r.get("chembl_id"): r for r in b}
    common = [k for k in ka if k in kb and k is not None]
    if not common:
        return {"error": "두 팔에 공통 화합물이 없다"}

    rnd = random.Random(seed)
    n = len(common)
    diffs: dict[str, list[float]] = {m: [] for m in METRICS}
    for _ in range(n_boot):
        idx = [common[rnd.randrange(n)] for _ in range(n)]
        va = compute(_labels_from([ka[k] for k in idx]))
        vb = compute(_labels_from([kb[k] for k in idx]))
        for m in METRICS:
            if va[m] is not None and vb[m] is not None:
                diffs[m].append(va[m] - vb[m])

    obs_a = compute(_labels_from([ka[k] for k in common]))
    obs_b = compute(_labels_from([kb[k] for k in common]))
    out = {"n_common": n, "metrics": {}}
    for m in METRICS:
        d = diffs[m]
        c = ci(d)
        # 0 이 구간 밖이면 차이가 유의하다고 본다
        sig = c["lo"] is not None and (c["lo"] > 0 or c["hi"] < 0)
        out["metrics"][m] = {
            "treatment": obs_a[m],
            "control": obs_b[m],
            "diff": (round(obs_a[m] - obs_b[m], 3)
                     if obs_a[m] is not None and obs_b[m] is not None else None),
            "ci95": c,
            "significant": sig,
        }
    return out


def load(path: Path) -> list[dict]:
    rows = [json.loads(l) for l in path.read_text().splitlines() if l.strip()]
    return [r for r in rows if r.get("score") is not None]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--treatment", required=True, help="docked_all.jsonl 경로")
    ap.add_argument("--control", default=None)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    t = load(Path(args.treatment))
    print(f"[*] {args.gene} 재검증 — treatment {len(t)}개 "
          f"(활성 {sum(1 for r in t if r.get('label') == 'known_active')}개)")
    print(f"    부트스트랩 {args.n_boot}회 · 순열검정 {args.n_perm}회")
    print()

    obs = compute(_labels_from(t))
    dist = bootstrap(t, args.n_boot)
    perm = permutation(t, args.n_perm)

    result = {"gene": args.gene, "n": len(t), "treatment": {}}
    print(f"    {'지표':8s} {'관측':>9s}  {'95% 신뢰구간':>20s}  {'p 값':>8s}")
    for m in METRICS:
        c = ci(dist[m])
        result["treatment"][m] = {"observed": obs[m], "ci95": c,
                                  "p_value": perm[m]["p_value"]}
        span = f"[{c['lo']}, {c['hi']}]" if c["lo"] is not None else "-"
        print(f"    {m:8s} {obs[m] if obs[m] is not None else 0:9.3f}  {span:>20s}  "
              f"{perm[m]['p_value']:8.4f}")

    if args.control:
        c = load(Path(args.control))
        print()
        print(f"[*] 짝지은 비교 — control {len(c)}개")
        pd = paired_diff(t, c, args.n_boot)
        result["paired"] = pd
        if "error" in pd:
            print(f"    [!] {pd['error']}")
        else:
            print(f"    공통 화합물 {pd['n_common']}개")
            print(f"    {'지표':8s} {'treatment':>10s} {'control':>9s} {'차이':>9s}  "
                  f"{'95% 신뢰구간':>22s}  판정")
            for m in METRICS:
                d = pd["metrics"][m]
                span = (f"[{d['ci95']['lo']}, {d['ci95']['hi']}]"
                        if d["ci95"]["lo"] is not None else "-")
                verdict = "유의" if d["significant"] else "불확실"
                print(f"    {m:8s} {d['treatment'] or 0:10.3f} {d['control'] or 0:9.3f} "
                      f"{d['diff'] or 0:9.3f}  {span:>22s}  {verdict}")

    out = Path(args.out) if args.out else ROOT / "results" / f"stats_{args.gene}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False))
    print()
    print(f"[+] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
