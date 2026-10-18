#!/usr/bin/env python3
"""
scripts/build_open_receptor.py — 유연 도킹에서 얻은 측쇄 형태를 고정한 수용체를 만든다.

왜
  2026-10-17 유연 측쇄 1단계(실행 38973161302)에서 M252/H293/F254 를 풀었더니
  **M252 만 움직였다** — 측쇄 RMSD 2.12 A, 말단 CE 3.96 A, belzutifan 접촉
  1.53 → 3.21 A. PHE254(0.15 A)·HIS293(0.37 A)은 자유를 줬는데도 제자리였다.

  움직이는 잔기가 하나뿐이라면, 그 회전 형태를 고정한 수용체를 만들어 평범한
  리지드 도킹을 돌리면 된다. 유연 도킹은 리간드마다 측쇄 탐색을 반복하므로
  실측 2.15배가 들지만, 형태를 미리 고정하면 **1배**다.

  나머지 9개 타깃의 단가가 여기서 갈린다. 이 수용체가 유연 도킹 이득의
  대부분을 회수하면 이후 전부를 리지드 비용으로 돌릴 수 있다.

주의
  이것은 '유연 도킹의 근사' 이지 대체가 아니다. 리간드마다 다른 측쇄 형태를
  요구하는 경우에는 한 형태로 고정한 수용체가 그 리간드를 놓친다. 그래서
  11월 실행에서 두 방식을 나란히 돌려 차이를 측정한다.

입력
  --receptor   원래 구조 (6CZW 체인 A)
  --pose       유연 도킹 출력 PDBQT (리간드 + 움직인 측쇄 좌표)
  --residues   고정할 잔기. 생략하면 포즈에 들어 있는 유연 잔기 전부
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def parse_flex_residues(pdbqt: Path, model: int = 1) -> dict[str, dict[str, tuple]]:
    """PDBQT 의 MODEL n 에서 BEGIN_RES 블록별 원자 좌표를 뽑는다."""
    out: dict[str, dict[str, tuple]] = {}
    cur: str | None = None
    inside = False
    for line in pdbqt.read_text().splitlines():
        if line.startswith(f"MODEL {model}"):
            inside = True
            continue
        if inside and line.startswith("ENDMDL"):
            break
        if not inside:
            continue
        if line.startswith("BEGIN_RES"):
            p = line.split()
            cur = f"{p[1]}:{p[2]}:{p[3]}"  # RESNAME:CHAIN:RESID
            out[cur] = {}
            continue
        if line.startswith("END_RES"):
            cur = None
            continue
        if cur and line.startswith(("ATOM", "HETATM")):
            out[cur][line[12:16].strip()] = (
                float(line[30:38]), float(line[38:46]), float(line[46:54])
            )
    return out


def substitute(receptor: Path, flex: dict[str, dict[str, tuple]],
               wanted: set[str] | None) -> tuple[list[str], dict]:
    """수용체 PDB 의 해당 잔기 원자 좌표를 유연 도킹 결과로 바꾼다.

    주쇄(N, CA, C, O)는 건드리지 않는다. 유연 측쇄 도킹은 측쇄만 움직이므로
    주쇄가 달라졌다면 파싱이 잘못된 것이다.
    """
    index: dict[tuple[str, str], tuple[str, dict[str, tuple]]] = {}
    for key, atoms in flex.items():
        rn, ch, ri = key.split(":")
        if wanted and f"{rn}{ri}" not in wanted and ri not in wanted:
            continue
        index[(ch, ri)] = (rn, atoms)

    BACKBONE = {"N", "CA", "C", "O", "OXT"}
    out, stats = [], {"substituted": {}, "skipped_backbone": 0, "unmatched": []}
    seen: dict[tuple[str, str], set[str]] = {}

    for line in receptor.read_text().splitlines():
        if not line.startswith("ATOM"):
            out.append(line)
            continue
        ch, ri, name = line[21], line[22:26].strip(), line[12:16].strip()
        hit = index.get((ch, ri))
        if hit is None or name not in hit[1]:
            out.append(line)
            continue
        rn, atoms = hit
        if line[17:20].strip() != rn:
            stats["unmatched"].append(f"{ch}:{ri} 구조는 {line[17:20].strip()}, 포즈는 {rn}")
            out.append(line)
            continue
        if name in BACKBONE:
            # 주쇄는 그대로 둔다. 다만 좌표가 다르면 경고 대상이다.
            old = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            if math.dist(old, atoms[name]) > 0.05:
                stats["skipped_backbone"] += 1
            out.append(line)
            continue
        x, y, z = atoms[name]
        old = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        out.append(f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}")
        k = f"{rn}{ri}"
        stats["substituted"].setdefault(k, []).append(
            {"atom": name, "shift": round(math.dist(old, (x, y, z)), 3)}
        )
        seen.setdefault((ch, ri), set()).add(name)

    for (ch, ri), (rn, atoms) in index.items():
        missing = set(atoms) - BACKBONE - seen.get((ch, ri), set())
        if missing:
            stats["unmatched"].append(f"{rn}{ri} 미치환 원자 {sorted(missing)}")
    return out, stats


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--pose", required=True, help="유연 도킹 출력 PDBQT")
    ap.add_argument("--residues", nargs="*", default=None,
                    help="고정할 잔기 (예: MET252 또는 252). 생략하면 전부")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-shift", type=float, default=0.5,
                    help="이 값 미만으로만 움직인 잔기는 치환하지 않는다. 움직이지 "
                         "않은 잔기를 바꾸면 노이즈만 들어간다")
    args = ap.parse_args()

    rec, pose = Path(args.receptor), Path(args.pose)
    for p in (rec, pose):
        if not p.exists():
            print(f"[!] 없음: {p}", file=sys.stderr)
            return 1

    flex = parse_flex_residues(pose)
    if not flex:
        print(f"[!] 포즈에서 유연 잔기를 찾지 못했다: {pose}", file=sys.stderr)
        return 1
    print(f"[*] 포즈의 유연 잔기: {', '.join(flex)}")

    # 실제로 움직인 잔기만 고른다
    wanted = set(args.residues) if args.residues else None
    moved: set[str] = set()
    BACKBONE = {"N", "CA", "C", "O", "OXT"}
    orig: dict[str, dict[str, tuple]] = {}
    for line in rec.read_text().splitlines():
        if line.startswith("ATOM"):
            key = f"{line[17:20].strip()}:{line[21]}:{line[22:26].strip()}"
            if key in flex:
                orig.setdefault(key, {})[line[12:16].strip()] = (
                    float(line[30:38]), float(line[38:46]), float(line[46:54]))
    for key, atoms in flex.items():
        rn, _, ri = key.split(":")
        o = orig.get(key, {})
        side = [n for n in atoms if n in o and n not in BACKBONE]
        if not side:
            continue
        rms = math.sqrt(sum(math.dist(o[n], atoms[n]) ** 2 for n in side) / len(side))
        tag = "치환" if rms >= args.min_shift else "제자리 — 건너뜀"
        print(f"    {rn}{ri:<4s} 측쇄 RMSD {rms:5.2f} A   {tag}")
        if rms >= args.min_shift and (wanted is None or f"{rn}{ri}" in wanted or ri in wanted):
            moved.add(f"{rn}{ri}")

    if not moved:
        print(f"[!] {args.min_shift} A 이상 움직인 잔기가 없다. 수용체를 바꿀 이유가 없다.",
              file=sys.stderr)
        return 1

    lines, stats = substitute(rec, flex, moved)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n")

    print()
    print(f"[+] {out}")
    for k, atoms in stats["substituted"].items():
        shifts = [a["shift"] for a in atoms]
        print(f"    {k}: 원자 {len(atoms)}개 치환, 최대 이동 {max(shifts):.2f} A")
    if stats["skipped_backbone"]:
        print(f"    [!] 주쇄 좌표가 다른 원자 {stats['skipped_backbone']}개 — 그대로 두었다")
    for u in stats["unmatched"]:
        print(f"    [!] {u}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
