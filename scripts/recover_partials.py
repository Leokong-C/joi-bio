#!/usr/bin/env python3
"""
잘린 샤드의 부분 결과를 회수한다 — **명시적으로만.**

## 왜 자동이 아닌가

`dock.py` 는 리간드 하나마다 `partial_shard{N}.jsonl` 에 append 한다. 샤드가
취소돼도 거기까지는 남는다. 2026-10·11 두 달 모두 러너 소비의 46% 가 취소된
샤드로 날아갔고, 이 파일이 그걸 막는다.

그런데 **부분 결과를 자동으로 평가셋에 넣으면 안 된다.** 불완전한 평가셋은
숫자가 없는 게 아니라 **그럴듯한 숫자가 나오고 방향까지 바뀌기** 때문이다.
2026-11-01 이 그랬다 — 치료군 7/8, 대조군 2/8 로 compare 를 돌렸더니
`BIAS_LIKELY` 가 나왔고 **대조군 AUC 가 치료군보다 높았다**(부호가 뒤집혔다).
활성 14개 표본의 잡음이었다.

그래서 이 도구는 **보여주기만 하고**, `--promote` 를 명시해야 옮긴다.
옮긴 뒤에도 게이트의 `COVERAGE_MIN` 검사가 불완전한 샤드 구성을 `INCOMPLETE`
로 덮는다. 우회하지 않는다.

## 쓸 자리

- 디버깅: 잘린 샤드가 어디까지 갔나, 어떤 리간드에서 느려졌나
- 재실행 설계: 실제 소요 시간 분포를 보고 timeout 을 정한다
  (11-01 의 실패 원인이 그 분포를 모른 채 추정한 것이었다)
- 부분 데이터로 **탐색적** 확인 — 보고서에는 올리지 않는다
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def read_rows(p: Path) -> list[dict]:
    out = []
    for line in p.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            out.append(json.loads(line))
        except json.JSONDecodeError:
            # 마지막 줄이 잘렸을 수 있다 — flush 하므로 드물지만 가능하다
            print(f"    [!] {p.name}: 줄 하나가 깨졌다 (마지막 줄로 보임). 건너뛴다")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-dir", required=True,
                    help="targets/<GENE> 또는 부분 파일이 있는 디렉터리")
    ap.add_argument("--promote", action="store_true",
                    help="부분 결과를 docked_shard 로 옮긴다. 명시해야만 동작한다")
    args = ap.parse_args()

    tdir = Path(args.target_dir)
    partials = sorted(tdir.glob("partial_shard*.jsonl"))
    complete = sorted(tdir.glob("docked_shard*.jsonl"))

    print(f"[*] {tdir}")
    print(f"    완주한 샤드 {len(complete)}개, 잘린 샤드 {len(partials)}개\n")

    if not partials:
        print("    잘린 샤드가 없다 — 회수할 것이 없다")
        return 0

    total_scored = 0
    print(f"    {'파일':28s} {'기록':>7s} {'스코어 있음':>10s} {'최고':>8s}")
    for p in partials:
        rows = read_rows(p)
        scored = [r for r in rows if r.get("score") is not None]
        total_scored += len(scored)
        best = min((r["score"] for r in scored), default=None)
        print(f"    {p.name:28s} {len(rows):7,d} {len(scored):10,d} "
              f"{best if best is not None else '-':>8}")

    print(f"\n    회수 가능한 스코어 {total_scored:,}개")

    if not args.promote:
        print("\n    --promote 없이는 옮기지 않는다.")
        print("    불완전한 평가셋은 숫자가 없는 게 아니라 그럴듯한 숫자를 낸다 —")
        print("    2026-11-01 에 치료군 7/8·대조군 2/8 로 돌렸더니 부호가 뒤집혔다.")
        return 0

    print("\n[!] --promote: 부분 결과를 docked_shard 로 옮긴다")
    print("    게이트의 COVERAGE_MIN 검사가 여전히 INCOMPLETE 로 덮을 것이다.")
    moved = 0
    for p in partials:
        shard = p.stem.replace("partial_shard", "")
        dest = tdir / f"docked_shard{shard}.jsonl"
        if dest.exists():
            print(f"    {p.name} → 건너뜀 ({dest.name} 이 이미 있다. 완주한 샤드를 덮지 않는다)")
            continue
        rows = [r for r in read_rows(p) if r.get("score") is not None]
        rows.sort(key=lambda r: r["score"])
        dest.write_text("\n".join(json.dumps(r, ensure_ascii=False) for r in rows) + "\n")
        print(f"    {p.name} → {dest.name} ({len(rows):,}개)")
        moved += 1

    print(f"\n[+] {moved}개 샤드를 옮겼다. 이 평가셋은 불완전하다 — "
          f"보고서에 올리기 전에 게이트 판정을 확인할 것")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
