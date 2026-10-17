#!/usr/bin/env python3
"""
modules.small_molecule.contacts — 후보를 포즈로 판정한다.

왜 점수만으로는 안 되는가
  2026-10-17 진단. EPAS1 상위 후보 20개를 뜯어보니:
    · 전부 decoy 라벨 — 검증 지표가 오답으로 세는 것들이 곧 신규 후보다
    · 알려진 활성보다 작다 (중원자 20 대 27, MW 321 대 402)
    · 크기-점수 상관 +0.730 — 클수록 점수가 나쁘다. PAS-B 공동이 작고 닫혀 있어
      큰 분자가 아예 들어가지 못한다(라이브러리의 22.6%)
    · LE 중앙값 0.479 — 통상 유망 기준 0.3 을 크게 넘는다. 중원자 20개가
      -10 kcal/mol 이면 50 nM 급인데 조각 화합물로는 비현실적이다

  종합하면 상위권이 '포켓을 잘 채운 분자' 가 아니라 '충돌을 피한 작은 분자'일
  수 있다. Vina 점수는 이 둘을 구분하지 못한다.

무엇으로 판정하는가
  승인약이 실제로 쓰는 상호작용을 재현하는지 본다. belzutifan 포즈를 기준으로:

    hbond_h293    HIS293 과 수소결합 거리에 극성 원자가 있는가
                  (Ren et al. 2022 이 특정한 상호작용)
    m252_pocket   M252 가 비켜주는 서브포켓을 점유하는가
    overlap       belzutifan 포즈와 부피가 겹치는가 — 같은 자리를 쓰는가
    buried        포켓 잔기와 4.5 A 이내 접촉 수 — 공동을 채우는가

  네 가지를 모두 만족해야 후보로 올린다. 점수가 아무리 좋아도 승인약의
  결합 양식을 재현하지 못하면 스코어링 함수의 빈틈일 가능성이 높다.

한계
  이것은 계산 필터이지 검증이 아니다. 통과한 화합물도 실험 전에는 가설이다.
  이 모듈의 목적은 '틀린 것을 걸러내는 것' 이지 '맞다고 증명하는 것' 이 아니다.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[2]

HBOND_MAX = 3.5      # 수소결합 공여체-수용체 중원자 거리 상한
CONTACT_MAX = 4.5    # 매몰 접촉으로 셀 거리
OVERLAP_MAX = 2.0    # 기준 리간드 원자와 '같은 자리'로 볼 거리
POLAR = {"N", "O", "S"}


def load_receptor(pdb: Path) -> list[tuple[str, int, str, tuple[float, float, float], str]]:
    out = []
    for line in pdb.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        el = (line[76:78].strip() or line[12:16].strip()[0]).upper()
        out.append((
            line[17:20].strip(),
            int(line[22:26]),
            line[12:16].strip(),
            (float(line[30:38]), float(line[38:46]), float(line[46:54])),
            el,
        ))
    return out


def mol_atoms(mol: Chem.Mol) -> list[tuple[str, tuple[float, float, float]]]:
    conf = mol.GetConformer()
    return [
        (a.GetSymbol(), tuple(conf.GetAtomPosition(a.GetIdx())))
        for a in mol.GetAtoms()
        if a.GetAtomicNum() > 1
    ]


def evaluate(mol: Chem.Mol, rec: list, ref: list | None, m252_centre) -> dict:
    atoms = mol_atoms(mol)
    d = math.dist

    # HIS293 극성 원자와의 최단 거리
    h293 = [a for a in rec if a[1] == 293 and a[4] in POLAR and a[2] not in ("N", "O")]
    if not h293:  # 측쇄 극성 원자가 없으면 주쇄라도
        h293 = [a for a in rec if a[1] == 293 and a[4] in POLAR]
    hb = min((d(a[3], p) for a in h293 for s, p in atoms if s in POLAR), default=None)

    # M252 서브포켓 — M252 측쇄가 비켜나며 열리는 공간
    m252 = min((d(m252_centre, p) for _, p in atoms), default=None)

    # 포켓 매몰도
    buried = sum(
        1 for _, p in atoms
        if any(d(a[3], p) <= CONTACT_MAX for a in rec)
    )

    # 기준 리간드와의 자리 겹침
    ov = None
    if ref:
        ov = sum(1 for _, p in atoms if any(d(q, p) <= OVERLAP_MAX for _, q in ref))
        ov = round(100 * ov / len(atoms), 1)

    return {
        "n_heavy": len(atoms),
        "hbond_h293_dist": round(hb, 2) if hb is not None else None,
        "m252_pocket_dist": round(m252, 2) if m252 is not None else None,
        "buried_contacts": buried,
        "buried_frac": round(buried / len(atoms), 3) if atoms else None,
        "overlap_pct": ov,
    }


def judge(m: dict, thresholds: dict) -> tuple[bool, list[str]]:
    fails = []
    if m["hbond_h293_dist"] is None or m["hbond_h293_dist"] > thresholds["hbond"]:
        fails.append(f"H293 수소결합 없음 ({m['hbond_h293_dist']} A)")
    if m["m252_pocket_dist"] is None or m["m252_pocket_dist"] > thresholds["m252"]:
        fails.append(f"M252 서브포켓 미점유 ({m['m252_pocket_dist']} A)")
    if m["overlap_pct"] is not None and m["overlap_pct"] < thresholds["overlap"]:
        fails.append(f"기준 리간드와 자리 불일치 ({m['overlap_pct']}%)")
    if m["buried_frac"] is not None and m["buried_frac"] < thresholds["buried"]:
        fails.append(f"포켓 매몰 부족 ({m['buried_frac']})")
    return (not fails), fails


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--targetdir", default=None)
    ap.add_argument("--reference", default="BELZUTIFAN",
                    help="자리 비교 기준 리간드의 이름 또는 ChEMBL ID")
    ap.add_argument("--hbond", type=float, default=HBOND_MAX)
    ap.add_argument("--m252", type=float, default=5.0)
    ap.add_argument("--overlap", type=float, default=40.0)
    ap.add_argument("--buried", type=float, default=0.85)
    args = ap.parse_args()

    tdir = Path(args.targetdir) if args.targetdir else ROOT / "targets" / args.gene
    sdf = tdir / "poses.sdf"
    if not sdf.exists():
        print(f"[!] 포즈 파일이 없다: {sdf}")
        return 1
    rec = load_receptor(tdir / "receptor.pdb")

    mols = [m for m in Chem.SDMolSupplier(str(sdf), removeHs=False) if m is not None]
    if not mols:
        print("[!] 포즈를 읽지 못했다")
        return 1

    # 1 순위 포즈만 쓴다
    best: dict[str, Chem.Mol] = {}
    for m in mols:
        cid = m.GetProp("chembl_id") if m.HasProp("chembl_id") else m.GetProp("_Name")
        if cid not in best:
            best[cid] = m

    ref_key = next(
        (k for k, m in best.items()
         if args.reference.upper() in (
             (m.GetProp("_Name") if m.HasProp("_Name") else "") + k).upper()),
        None,
    )
    ref = mol_atoms(best[ref_key]) if ref_key else None
    if ref is None:
        print(f"[!] 기준 리간드 {args.reference} 를 포즈에서 찾지 못했다. 자리 비교를 건너뛴다.")

    # M252 서브포켓 중심 — 측쇄 말단(SD, CE)의 중점
    m252_atoms = [a for a in rec if a[1] == 252 and a[2] in ("SD", "CE")]
    if not m252_atoms:
        print("[!] M252 SD/CE 를 찾지 못했다")
        return 1
    m252_centre = tuple(
        sum(a[3][i] for a in m252_atoms) / len(m252_atoms) for i in range(3)
    )

    th = {"hbond": args.hbond, "m252": args.m252,
          "overlap": args.overlap, "buried": args.buried}

    out = []
    for cid, m in best.items():
        metrics = evaluate(m, rec, ref, m252_centre)
        score = float(m.GetProp("vina_score")) if m.HasProp("vina_score") else None
        smi = Chem.MolToSmiles(Chem.RemoveHs(m))
        rd = Chem.MolFromSmiles(smi)
        metrics.update(
            chembl_id=cid,
            name=(m.GetProp("_Name") if m.HasProp("_Name") else cid),
            group=(m.GetProp("group") if m.HasProp("group") else None),
            vina_score=score,
            ligand_efficiency=(
                round(-score / metrics["n_heavy"], 3)
                if score and metrics["n_heavy"] else None
            ),
            mw=round(Descriptors.MolWt(rd), 1) if rd else None,
            clogp=round(Crippen.MolLogP(rd), 2) if rd else None,
            smiles=smi,
        )
        ok, fails = judge(metrics, th)
        metrics["contact_pass"] = ok
        metrics["contact_fails"] = fails
        out.append(metrics)

    out.sort(key=lambda r: (not r["contact_pass"], r["vina_score"] or 0))
    (tdir / "contacts.json").write_text(
        json.dumps({"gene": args.gene, "thresholds": th,
                    "reference": ref_key, "compounds": out},
                   indent=2, ensure_ascii=False)
    )

    passed = [r for r in out if r["contact_pass"]]
    print(f"[*] {args.gene} 포즈 접촉 판정 — {len(out)}개 중 {len(passed)}개 통과")
    print(f"    기준: H293 수소결합 ≤{th['hbond']} A · M252 서브포켓 ≤{th['m252']} A · "
          f"자리 겹침 ≥{th['overlap']}% · 매몰 ≥{th['buried']}")
    print()
    hdr = f"{'화합물':22s} {'점수':>7s} {'LE':>6s} {'H293':>6s} {'M252':>6s} {'겹침':>6s} {'판정':>5s}"
    print(hdr)
    for r in out:
        print(f"{(r['name'] or r['chembl_id'])[:22]:22s} "
              f"{r['vina_score'] or 0:7.2f} {r['ligand_efficiency'] or 0:6.3f} "
              f"{r['hbond_h293_dist'] or 99:6.2f} {r['m252_pocket_dist'] or 99:6.2f} "
              f"{r['overlap_pct'] or 0:6.1f} {'통과' if r['contact_pass'] else '탈락':>5s}")
        if not r["contact_pass"]:
            print(f"    └ {' · '.join(r['contact_fails'])}")
    print()
    print(f"[+] {tdir / 'contacts.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
