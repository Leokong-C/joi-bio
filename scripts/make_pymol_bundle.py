#!/usr/bin/env python3
"""
PyMOL 로 열 수 있는 묶음을 만든다.

Vina 가 내놓는 PDBQT 는 PDB 뒤에 부분전하와 AD4 원자타입 두 칼럼이 더 붙은
형식이라 PyMOL 이 깨끗하게 읽지 못한다. 그 두 칼럼을 떼고 표준 PDB 로 바꾼다.

유연 측쇄 도킹 결과는 한 파일에 리간드(UNL)와 움직인 측쇄(BEGIN_RES/END_RES)가
MODEL 단위로 같이 들어 있다. 최적 포즈(MODEL 1)만 꺼내 리간드와 측쇄를 나눈다.

출력: exp/<gene>/pymol/
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

ELEM = {
    "A": "C", "C": "C", "N": "N", "NA": "N", "OA": "O", "O": "O",
    "SA": "S", "S": "S", "HD": "H", "H": "H", "F": "F", "Cl": "CL",
    "CL": "CL", "Br": "BR", "BR": "BR", "I": "I", "P": "P",
}


def _pdb_line(line: str, resname: str | None, resseq: int | None, serial: int) -> str:
    """PDBQT 한 줄을 표준 PDB 로. 전하/타입 칼럼(67~)을 떼고 element 를 채운다."""
    ad_type = line[77:79].strip() if len(line) > 77 else ""
    elem = ELEM.get(ad_type, ad_type[:1] if ad_type else line[12:14].strip()[:1])
    rn = (resname or line[17:20]).ljust(3)[:3]
    rs = resseq if resseq is not None else int(line[22:26] or 1)
    return (
        f"{line[:6]}{serial:5d} {line[12:16]}{line[16]}{rn} {line[21]}"
        f"{rs:4d}{line[26]}   {line[30:54]}  1.00  0.00          {elem:>2s}\n"
    )


def convert_flex_pose(src: Path, out_lig: Path, out_sc: Path) -> tuple[int, int]:
    """MODEL 1 에서 리간드와 움직인 측쇄를 갈라 각각 PDB 로 쓴다."""
    lines = src.read_text().splitlines()
    try:
        start = lines.index("MODEL 1")
    except ValueError:
        start = 0
    end = next((i for i in range(start + 1, len(lines)) if lines[i].startswith("ENDMDL")), len(lines))
    block = lines[start:end]

    lig, sc = [], []
    in_res: tuple[str, int] | None = None
    for ln in block:
        if ln.startswith("BEGIN_RES"):
            p = ln.split()
            in_res = (p[1], int(p[3]))
            continue
        if ln.startswith("END_RES"):
            in_res = None
            continue
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        if in_res is None:
            lig.append(ln)
        else:
            sc.append((ln, in_res))

    out_lig.write_text(
        "".join(_pdb_line(ln, "LIG", 1, i + 1) for i, ln in enumerate(lig)) + "END\n"
    )
    out_sc.write_text(
        "".join(
            _pdb_line(ln, rn, rs, i + 1) for i, (ln, (rn, rs)) in enumerate(sc)
        )
        + "END\n"
    )
    return len(lig), len(sc)


PML = """\
# Joi_bio — EPAS1 (HIF-2 alpha) PAS-B
#
# PyMOL 에서:  pymol joi_bio.pml
# 또는 PyMOL 안에서:  @joi_bio.pml
#
# 장면 세 개가 들어 있다. F1~F3 또는 scene 메뉴로 전환.

reinitialize
bg_color white
set ray_opaque_background, 0

# ---------- 구조 ----------
load receptor_rigid.pdb,     rigid
load receptor_m252open.pdb,  m252open
load poses_rigid.sdf,        poses_rigid

# 유연 측쇄 도킹 결과 (벤치마크 3개)
load flex_belzutifan_lig.pdb,  flex_belz_lig
load flex_belzutifan_sc.pdb,   flex_belz_sc
load flex_PT2399_lig.pdb,      flex_pt2399_lig
load flex_PT2385_lig.pdb,      flex_pt2385_lig

# ---------- 기본 표현 ----------
hide everything
show cartoon, rigid or m252open
color grey80, rigid
color palecyan, m252open
set cartoon_transparency, 0.45, rigid
set cartoon_transparency, 0.45, m252open

# 문헌이 지목한 세 잔기 (Ren et al. 2022, Mol Pharmacol 102:240)
select key_res_rigid, rigid and resi 252+254+293
select key_res_open,  m252open and resi 252+254+293
show sticks, key_res_rigid or key_res_open
color grey50, key_res_rigid and elem C
color deepteal, key_res_open and elem C
util.cnc("key_res_rigid or key_res_open")

# ---------- 장면 1: 리지드 도킹 포즈 ----------
disable m252open
disable flex_belz_lig
disable flex_belz_sc
disable flex_pt2399_lig
disable flex_pt2385_lig
show sticks, poses_rigid
color yellow, poses_rigid and elem C
util.cnc("poses_rigid")
orient rigid and resi 252+254+293
zoom rigid and resi 252+254+293, 8
scene overview, store, message="리지드 도킹 — 벤치마크 3개 + 상위 히트 20개 (노란색). 회색 막대가 M252/F254/H293"

# ---------- 장면 2: M252 측쇄 뒤집힘 ----------
disable poses_rigid
enable m252open
show sticks, flex_belz_lig
color orange, flex_belz_lig and elem C
util.cnc("flex_belz_lig")
enable flex_belz_lig
distance d_rigid, rigid and resi 252 and name CE, flex_belz_lig
color red, d_rigid
hide labels, d_rigid
orient rigid and resi 252
zoom rigid and resi 252, 6
scene m252_flip, store, message="M252 뒤집힘 — 회색이 원본 6CZW, 청록이 유연 도킹 후. 측쇄 RMSD 2.37 A, CE 원자 3.96 A 이동. 접촉 거리 1.53 -> 3.21 A"

# ---------- 장면 3: 유연 측쇄 도킹 결과 ----------
disable rigid
enable m252open
show sticks, flex_belz_sc
color salmon, flex_belz_sc and elem C
util.cnc("flex_belz_sc")
enable flex_belz_sc
enable flex_pt2399_lig
enable flex_pt2385_lig
show sticks, flex_pt2399_lig or flex_pt2385_lig
color marine, flex_pt2399_lig and elem C
color purple, flex_pt2385_lig and elem C
util.cnc("flex_pt2399_lig or flex_pt2385_lig")
orient m252open and resi 252+254+293
zoom m252open and resi 252+254+293, 8
scene flex_result, store, message="유연 측쇄 도킹 — belzutifan(주황) PT2399(파랑) PT2385(보라). 연어색이 움직인 측쇄. 세 잔기를 풀었으나 실제로 움직인 것은 M252 하나"

# ---------- 마무리 ----------
scene overview, recall
set stick_radius, 0.14
set cartoon_side_chain_helper, on
print ""
print "=== Joi_bio EPAS1 ==="
print "장면: overview / m252_flip / flex_result   (F1~F3 또는 scene 메뉴)"
print ""
print "유연 측쇄 1단계 결과 (실행 38973161302, 표본 53개)"
print "  AUC      0.7667 -> 0.8478"
print "  BEDROC   0.9656 -> 0.9978"
print "  벤치마크 순위 중앙값 +9 계단"
print "  판정     INDUCED_FIT_SUPPORTED"
print ""
print "움직인 잔기는 M252 하나뿐이다 (측쇄 RMSD 2.37 A)."
print "PHE254 는 0.16 A, HIS293 은 0.40 A 로 사실상 제자리."
print ""
"""

README = """\
Joi_bio — PyMOL 묶음 (EPAS1 / HIF-2 alpha PAS-B)
생성 2026-10-23

여는 법
  pymol joi_bio.pml
  또는 PyMOL 안에서  @joi_bio.pml

파일
  receptor_rigid.pdb        6CZW chain A, 880 원자. 리지드 도킹에 쓴 원본 수용체
  receptor_m252open.pdb     같은 구조에서 M252 측쇄만 유연 도킹 결과로 바꾼 것
                            880 원자 보존, MET252 의 CG/SD/CE 만 좌표가 다르다
  poses_rigid.sdf           리지드 도킹 포즈 23개 (벤치마크 3 + 상위 히트 20)
                            각 분자에 vina_score, group, chembl_id 속성이 붙어 있다
  flex_*_lig.pdb            유연 측쇄 도킹의 최적 포즈 (리간드)
  flex_belzutifan_sc.pdb    같은 도킹에서 움직인 측쇄 3개
  joi_bio.pml               장면 3개가 든 PyMOL 스크립트

무엇을 보면 되나
  M252 측쇄 뒤집힘이 이 묶음의 핵심이다. Ren et al. 2022 (Mol Pharmacol
  102:240) 가 belzutifan 결합 시 M252 가 포켓 밖으로 뒤집히고 H293 이
  수소결합, F254 가 pi-pi 상호작용을 한다고 보고했다. M252A 변이에서
  Ki 가 23 nM -> 3.5 uM 로 약 150배 악화된다.

  우리는 이 세 잔기를 유연하게 풀고 도킹했는데, 실제로 움직인 것은
  M252 하나였다. 측쇄 RMSD 2.37 A, CE 원자 3.96 A, SD 2.43 A 이동.
  리간드와의 접촉 거리가 1.53 A 에서 3.21 A 로 벌어진다.
  PHE254 는 0.16 A, HIS293 은 0.40 A 로 사실상 제자리다.

주의
  Vina 스코어는 순위용이다. 결합 자유에너지 추정에 쓰지 않는다.
  이 포즈들은 실험으로 확인된 바 없다.
  EPAS1 치료군은 BIAS_LIKELY 판정 — 상위권 순위의 75%가 구조 선택에서 나왔다.
  자세한 내용은 docs/context-of-use.md 참조.
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", default="EPAS1")
    ap.add_argument("--root", default=".")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    exp = root / "exp" / args.gene
    out = exp / "pymol"
    out.mkdir(parents=True, exist_ok=True)

    shutil.copy(exp / "treatment" / "receptor.pdb", out / "receptor_rigid.pdb")
    shutil.copy(exp / "m252open" / "receptor.pdb", out / "receptor_m252open.pdb")
    shutil.copy(exp / "poses_v2.sdf", out / "poses_rigid.sdf")

    names = {
        "CHEMBL4585668": "belzutifan",
        "CHEMBL4165848": "PT2399",
        "CHEMBL4173075": "PT2385",
    }
    for cid, label in names.items():
        src = exp / "flex" / "poses" / f"{cid}.pdbqt"
        if not src.exists():
            print(f"    [!] 없음: {src.name}")
            continue
        nl, ns = convert_flex_pose(
            src, out / f"flex_{label}_lig.pdb", out / f"flex_{label}_sc.pdb"
        )
        print(f"    {label:12s} 리간드 {nl:3d} 원자, 측쇄 {ns:3d} 원자")

    (out / "joi_bio.pml").write_text(PML)
    (out / "README.txt").write_text(README)
    print(f"\n[+] {out}")
    for f in sorted(out.iterdir()):
        print(f"    {f.name:28s} {f.stat().st_size:>8,d} bytes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
