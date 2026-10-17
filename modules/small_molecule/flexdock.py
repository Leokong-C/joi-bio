#!/usr/bin/env python3
"""
modules.small_molecule.flexdock — 유연 측쇄 도킹 (메커니즘 검증).

왜
  리지드 도킹은 리간드가 들어올 때의 측쇄 재배열을 모델링하지 못한다. EPAS1
  6CZW 에서 belzutifan 포즈를 재보니 MET252 측쇄가 2.05 Å 거리에 있다 —
  탄소-탄소 반데르발스 접촉이 보통 3.4~4.0 Å 이므로 이것은 충돌이다.
  Ren et al. 2022 (Mol Pharmacol 102:240) 은 belzutifan 이 결합하면 M252 가
  포켓 밖으로 뒤집힌다고 보고했고, H293 수소결합·F254 π–π 도 특정했다.
  우리 구조에서 셋 다 접촉 범위(2.05 / 3.30 / 3.91 Å)에 있다.

무엇을 묻는가
  이 단계는 농축도를 재계산하지 않는다. 유연 도킹은 자유도가 늘어 스코어가
  체계적으로 좋아지므로, 일부만 유연하게 돌린 점수를 리지드 점수와 한 줄에
  세우면 순위가 무의미해진다. 대신 같은 화합물 집합을 두 방식으로 모두 돌려
  짝지어 비교한다.

    belzutifan 의 개선폭이 다른 화합물의 개선폭보다 큰가?

  전부 똑같이 좋아지면 척도 이동일 뿐 순위는 그대로다. 벤치마크만 선택적으로
  좋아져야 유도 적합이 원인이라는 뜻이 된다.

덤으로 얻는 것
  리지드 대비 소요 배수(M). 전체 라이브러리 재실행 비용을 추정 없이 계산할 수 있다.
"""
from __future__ import annotations

import argparse
import json
import os
import random
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from meeko import MoleculePreparation, PDBQTWriterLegacy
from vina import Vina

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[2]


def _candidate_cmds(pdb: Path, outstem: Path, flex_res: list[str]) -> list[list[str]]:
    """mk_prepare_receptor 의 인자명은 meeko 버전마다 다르다.

    0.6.x  --read_pdb / -o / -p / -f
    0.5.x  --macromolecule / -o / -f
    실패 시 다음 후보로 넘어간다. 한 번의 CI 실패가 100 러너분이므로 여기서 흡수한다.
    """
    exe = ["mk_prepare_receptor.py"]
    if shutil.which("mk_prepare_receptor.py") is None:
        if shutil.which("mk_prepare_receptor"):
            exe = ["mk_prepare_receptor"]
        else:
            exe = [sys.executable, "-m", "meeko.cli.mk_prepare_receptor"]
    # -f / --flexres 는 인자를 하나만 받는다. 잔기 여러 개는 쉼표로 묶어 넘긴다.
    # (2026-10-17 실행 38971290707 에서 확인: 공백 구분은 unrecognized arguments)
    csv = ",".join(flex_res)
    rep: list[str] = []
    for r in flex_res:
        rep += ["-f", r]
    # 6CZW 의 SER292·TYR342 는 대체 위치(altloc A)를 갖는다. meeko 는 altloc 라벨이
    # 보이면 어느 변이를 쓸지 모르겠다며 거부하므로 --default_altloc 으로 지정한다.
    # (2026-10-17 실행 38972510507 에서 확인)
    head = ["--read_pdb", str(pdb), "-o", str(outstem), "-p"]
    base = [
        head + ["-f", csv, "--allow_bad_res", "--default_altloc", "A"],
        head + ["-f", csv, "--default_altloc", "A"],
        head + rep + ["--allow_bad_res", "--default_altloc", "A"],
        head + ["-f", csv, "--allow_bad_res"],
        ["-i", str(pdb), "-o", str(outstem), "-p", "-f", csv, "--allow_bad_res",
         "--default_altloc", "A"],
    ]
    return [exe + a for a in base]


def prepare_flex_receptor(pdb: Path, outstem: Path, flex_res: list[str]) -> tuple[Path, Path] | None:
    """Meeko 로 수용체를 rigid + flex 두 PDBQT 로 나눈다.

    flex_res 는 "A:252" 같은 체인:잔기번호 형식.
    """
    for i, cmd in enumerate(_candidate_cmds(pdb, outstem, flex_res), 1):
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
        except (FileNotFoundError, subprocess.TimeoutExpired) as e:
            print(f"[!] 후보 {i} 실행 불가: {e}", file=sys.stderr)
            continue
        rigid = Path(f"{outstem}_rigid.pdbqt")
        flex = Path(f"{outstem}_flex.pdbqt")
        if r.returncode == 0 and rigid.exists() and flex.exists():
            if i > 1:
                print(f"    (인자 후보 {i} 에서 성공)")
            nflex = sum(1 for l in flex.read_text().splitlines() if l.startswith("BEGIN_RES"))
            print(f"    rigid {rigid.stat().st_size // 1024} KB · "
                  f"flex {flex.stat().st_size} B · 유연 잔기 {nflex}개")
            if nflex != len(flex_res):
                print(f"[!] 요청한 유연 잔기 {len(flex_res)}개 중 {nflex}개만 반영됐다. "
                      "잔기 번호·체인을 확인할 것.", file=sys.stderr)
            return rigid, flex
        print(f"[!] 후보 {i} 실패 (exit {r.returncode})", file=sys.stderr)
        tail = (r.stderr or r.stdout or "").strip().splitlines()[-16:]
        for line in tail:
            print(f"      {line}", file=sys.stderr)

    made = sorted(p.name for p in outstem.parent.glob(f"{outstem.name}*"))
    print(f"[!] 모든 인자 후보가 실패했다. 생성된 파일: {made}", file=sys.stderr)
    return None


def ligand_pdbqt(smiles: str, seed: int = 42) -> str | None:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        p = AllChem.ETKDGv3()
        p.randomSeed = seed
        if AllChem.EmbedMolecule(mol, p) != 0:
            p.useRandomCoords = True
            if AllChem.EmbedMolecule(mol, p) != 0:
                return None
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=600)
        except Exception:
            pass
        setups = MoleculePreparation().prepare(mol)
        if not setups:
            return None
        txt, ok, _ = PDBQTWriterLegacy.write_string(setups[0])
        return txt if ok else None
    except Exception:
        return None


def dock_one(job) -> dict:
    rigid, flex, center, size, rec, exh, posedir = job
    out = {k: rec[k] for k in ("chembl_id", "pref_name", "label", "smiles", "group", "source")
           if k in rec}
    t0 = time.time()
    try:
        lig = ligand_pdbqt(rec["smiles"])
        if lig is None:
            out.update(flex_score=None, error="ligand_prep_failed")
            return out
        v = Vina(sf_name="vina", cpu=1, seed=42, verbosity=0)
        v.set_receptor(rigid_pdbqt_filename=rigid, flex_pdbqt_filename=flex)
        v.set_ligand_from_string(lig)
        v.compute_vina_maps(center=list(center), box_size=list(size))
        v.dock(exhaustiveness=exh, n_poses=3)
        e = v.energies(n_poses=3)
        out.update(flex_score=float(e[0][0]), seconds=round(time.time() - t0, 1))
        # 벤치마크는 포즈까지 남긴다. 유연 도킹 출력에는 움직인 측쇄 좌표가 함께
        # 들어 있어, 리지드 구조와 겹쳐 보면 M252 가 실제로 뒤집혔는지 눈으로 확인된다.
        if posedir:
            p = Path(posedir) / f"{rec['chembl_id']}.pdbqt"
            try:
                v.write_poses(str(p), n_poses=3, overwrite=True)
                out["pose_file"] = p.name
            except Exception as pe:  # noqa: BLE001
                out["pose_error"] = str(pe)[:120]
        return out
    except Exception as ex:  # noqa: BLE001
        out.update(flex_score=None, error=f"{type(ex).__name__}: {ex}"[:200],
                   seconds=round(time.time() - t0, 1))
        return out


def select_set(rows: list[dict], n_active: int, n_decoy: int, seed: int = 42) -> list[dict]:
    """벤치마크 전부 + 활성 일부 + decoy 일부. 짝지어 비교할 집합."""
    random.seed(seed)
    bench = [r for r in rows if r.get("source") == "benchmark_spikein"]
    bench_ids = {r["chembl_id"] for r in bench}
    acts = [r for r in rows if r.get("label") == "known_active" and r["chembl_id"] not in bench_ids]
    decs = [r for r in rows if r.get("label") == "decoy"]
    # 리지드 스코어 분포 전체를 덮도록 균등 추출
    acts.sort(key=lambda r: r["score"])
    decs.sort(key=lambda r: r["score"])
    pick = lambda xs, k: [xs[round(i * (len(xs) - 1) / max(1, k - 1))] for i in range(min(k, len(xs)))]
    return bench + pick(acts, n_active) + pick(decs, n_decoy)


def _median(xs: list[float]) -> float | None:
    if not xs:
        return None
    xs = sorted(xs)
    n = len(xs)
    return round(xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2, 3)


def _group(r: dict) -> str:
    if r.get("source") == "benchmark_spikein" or r.get("group") == "benchmark":
        return "benchmark"
    return "active" if r.get("label") == "known_active" else "decoy"


def analyse(rows: list[dict]) -> dict:
    """짝지어 비교. 유연화가 활성/decoy 분리도를 개선했는가?

    왜 '중앙값 개선폭' 으로 판정하지 않는가
      2026-10-17 실행 38973161302 에서 배운 것. 처음에는 decoy 의 중앙값 개선폭을
      기준선으로 두고 벤치마크가 그보다 더 좋아졌는지를 봤다. 그 통계는 틀렸다.

      라이브러리의 22.6% 는 리지드에서 스코어가 양수였다 — Vina 에서 양수는
      '포켓에 들어가지 못했다' 는 뜻이다. 측쇄를 풀면 포켓이 열려 이들이 들어오고,
      스코어가 +17 에서 -5 로 뛴다. 개선폭 20 이 넘는 값이 중앙값을 끌어당긴다.
      그런데 이것은 스코어링이 좋아진 것이 아니라 '입장' 이라는 다른 현상이다.
      실제로 이들은 유연 상위권에 끼지 못한다(중앙값 -5.6 대 상위 10위 경계 -10.6).

    그래서 무엇으로 판정하는가
      같은 화합물 집합을 두 방식으로 돌렸으므로, 같은 집합 안의 활성/decoy
      분리도 변화를 직접 본다. AUC 와 BEDROC 의 짝지어 변화가 1차 지표다.
      이 지표는 양수 스코어 화합물에 휘둘리지 않는다 — 순위만 쓰기 때문이다.
      벤치마크 순위 변화를 2차로 보고, 새로 입장한 화합물 수는 위험 신호로 따로 센다.

    주의
      이 집합은 리지드 스코어 분포를 따라 균등 추출한 것이라 활성 비율이 높다
      (43%). AUC·BEDROC 의 절대값은 전체 라이브러리 값과 비교할 수 없다.
      같은 집합 안에서의 '변화' 만 유효하다.
    """
    from core.validation.gate import bedroc, roc_auc

    def sep(key: str) -> tuple[float | None, float | None]:
        s = sorted(rows, key=lambda r: r[key])
        lab = [1 if r.get("label") == "known_active" else 0 for r in s]
        if not any(lab) or all(lab):
            return None, None
        return round(roc_auc(lab), 4), round(bedroc(lab), 4)

    r_auc, r_bed = sep("rigid_score")
    f_auc, f_bed = sep("flex_score")

    by: dict[str, list[float]] = {"benchmark": [], "active": [], "decoy": []}
    for r in rows:
        by[_group(r)].append(r["delta"])

    # 리지드에서 들어가지 못했던 화합물 — 개선폭 통계에서 분리한다
    admitted = [r for r in rows if r["rigid_score"] > 0]
    fitted = [r for r in rows if r["rigid_score"] <= 0]

    out: dict = {
        "n": len(rows),
        "n_by_group": {k: len(v) for k, v in by.items()},
        "separation": {
            "rigid": {"auc": r_auc, "bedroc": r_bed},
            "flex": {"auc": f_auc, "bedroc": f_bed},
            "auc_delta": round(f_auc - r_auc, 4) if (f_auc and r_auc) else None,
            "bedroc_delta": round(f_bed - r_bed, 4) if (f_bed and r_bed) else None,
        },
        "pocket_admission": {
            "n": len(admitted),
            "share_pct": round(100 * len(admitted) / len(rows), 1) if rows else None,
            "n_active": sum(1 for r in admitted if r.get("label") == "known_active"),
            "flex_score_median": _median([r["flex_score"] for r in admitted]),
            "flex_score_best": min((r["flex_score"] for r in admitted), default=None),
            "note": "리지드 스코어가 양수였던 화합물. 포켓에 들어가지 못했다가 측쇄가 "
                    "풀리며 입장한 것으로, 스코어링 개선과는 다른 현상이다.",
        },
        # 원래 들어갔던 화합물만으로 본 개선폭 — 이쪽이 '재스코어링' 효과다
        "median_delta_fitted_only": {
            k: _median([r["delta"] for r in fitted if _group(r) == k])
            for k in ("benchmark", "active", "decoy")
        },
        "median_delta_all": {k: _median(v) for k, v in by.items()},
    }

    rr = sorted(rows, key=lambda r: r["rigid_score"])
    fr = sorted(rows, key=lambda r: r["flex_score"])
    rpos = {r["chembl_id"]: i + 1 for i, r in enumerate(rr)}
    fpos = {r["chembl_id"]: i + 1 for i, r in enumerate(fr)}
    out["benchmarks"] = sorted(
        (
            {
                "name": r.get("pref_name") or r["chembl_id"],
                "rigid_score": r["rigid_score"],
                "flex_score": r["flex_score"],
                "delta": r["delta"],
                "rank_rigid": rpos[r["chembl_id"]],
                "rank_flex": fpos[r["chembl_id"]],
                "rank_change": rpos[r["chembl_id"]] - fpos[r["chembl_id"]],
            }
            for r in rows
            if _group(r) == "benchmark"
        ),
        key=lambda b: b["rank_flex"],
    )
    # 유연 상위 10위의 구성 — 새로 입장한 화합물이 상위권을 잠식했는가
    top = fr[:10]
    out["flex_top10"] = {
        "n_active": sum(1 for r in top if r.get("label") == "known_active"),
        "n_benchmark": sum(1 for r in top if _group(r) == "benchmark"),
        "n_newly_admitted": sum(1 for r in top if r["rigid_score"] > 0),
        "score_cutoff": round(top[-1]["flex_score"], 2) if top else None,
    }

    d_auc = out["separation"]["auc_delta"]
    rank_gain = _median([b["rank_change"] for b in out["benchmarks"]])
    if d_auc is None:
        out["verdict"] = "INCONCLUSIVE"
        out["verdict_reason"] = "활성/decoy 분리도를 계산할 수 없다."
    elif d_auc >= 0.03 and (rank_gain or 0) > 0:
        out["verdict"] = "INDUCED_FIT_SUPPORTED"
        out["verdict_reason"] = (
            f"같은 집합에서 AUC 가 {d_auc:+.4f} 변했고 벤치마크 순위가 중앙값 "
            f"{rank_gain:+.0f}계단 올랐다. 새로 입장한 화합물은 상위 10위에 "
            f"{out['flex_top10']['n_newly_admitted']}개만 들어왔다. "
            "측쇄 재배열이 순위를 실제로 개선한다 — 전체 라이브러리 재실행이 정당하다."
        )
    elif d_auc >= 0.0:
        out["verdict"] = "INDUCED_FIT_WEAK"
        out["verdict_reason"] = (
            f"AUC 변화가 {d_auc:+.4f} 에 그친다. 개선 방향이긴 하나 작다. "
            "전체 재실행의 기대 이득이 비용에 비해 불확실하다."
        )
    else:
        out["verdict"] = "SEPARATION_DEGRADED"
        out["verdict_reason"] = (
            f"AUC 가 {d_auc:+.4f} 로 나빠졌다. 자유도가 늘면서 decoy 가 활성보다 "
            "더 이득을 봤다는 뜻이다. 유연화는 이 포켓에서 역효과다."
        )
    return out


def report(a: dict) -> None:
    s = a["separation"]
    pa = a["pocket_admission"]
    print()
    print("[*] 짝지어 비교 — 유연 측쇄가 활성/decoy 분리도를 개선했는가?")
    print(f"    같은 {a['n']}개 집합 내 분리도 (절대값은 전체 라이브러리와 비교 불가)")
    print(f"      리지드  AUC {s['rigid']['auc']}  BEDROC {s['rigid']['bedroc']}")
    print(f"      유연    AUC {s['flex']['auc']}  BEDROC {s['flex']['bedroc']}")
    print(f"      변화    AUC {s['auc_delta']:+}  BEDROC {s['bedroc_delta']:+}")
    print()
    print(f"    포켓 입장 — 리지드에서 스코어 양수였던 화합물")
    print(f"      {pa['n']}개 ({pa['share_pct']}%), 그중 활성 {pa['n_active']}개")
    print(f"      유연 점수 중앙값 {pa['flex_score_median']}, 최고 {pa['flex_score_best']}")
    print()
    print(f"    유연 상위 10위 구성")
    t = a["flex_top10"]
    print(f"      활성 {t['n_active']}/10, 벤치마크 {t['n_benchmark']}개, "
          f"새로 입장 {t['n_newly_admitted']}개, 경계 {t['score_cutoff']}")
    print()
    print(f"    개선폭 중앙값 — 원래 들어갔던 화합물만 (입장 효과 제외)")
    for k in ("benchmark", "active", "decoy"):
        print(f"      {k:10s} Δ={a['median_delta_fitted_only'][k]}")
    print()
    print("    벤치마크별 (53개 집합 내 순위)")
    for b in a["benchmarks"]:
        print(f"      {b['name']:16s} {b['rigid_score']:7.2f} → {b['flex_score']:7.2f}  "
              f"Δ{b['delta']:+.2f}  순위 {b['rank_rigid']}위 → {b['rank_flex']}위 "
              f"({b['rank_change']:+d})")
    print()
    print(f"    판정: {a['verdict']}")
    print(f"    {a['verdict_reason']}")


def main() -> int:
    # stdout 이 버퍼링되면 stderr 의 실패 메시지가 먼저 찍혀 로그 순서가 뒤집힌다.
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except Exception:
        pass
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--targetdir", required=True, help="리지드 결과가 있는 디렉터리")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--flex-res", nargs="+", default=["A:252", "A:293", "A:254"],
                    help="유연화할 잔기 (체인:번호)")
    ap.add_argument("--n-active", type=int, default=20)
    ap.add_argument("--n-decoy", type=int, default=30)
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--prep-only", action="store_true",
                    help="수용체 분리까지만 하고 끝낸다 (CI 사전점검용)")
    args = ap.parse_args()

    tdir, outdir = Path(args.targetdir), Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    meta = json.loads((tdir / "target_meta.json").read_text())
    pocket = meta["pocket"]

    print(f"[*] {args.gene} 유연 측쇄 도킹")
    print(f"    유연화 잔기: {', '.join(args.flex_res)}")
    prep = prepare_flex_receptor(tdir / "receptor.pdb", outdir / "receptor", args.flex_res)
    if prep is None:
        return 1
    rigid, flex = prep
    if args.prep_only:
        print("[+] 사전점검 통과 — 수용체 분리 성공")
        return 0

    rows = [json.loads(l) for l in (tdir / "docked_all.jsonl").read_text().splitlines() if l.strip()]
    rows = [r for r in rows if r.get("score") is not None]
    sel = select_set(rows, args.n_active, args.n_decoy)
    print(f"    대상 {len(sel)}개 "
          f"(벤치마크 {sum(1 for r in sel if r.get('source')=='benchmark_spikein')}, "
          f"활성 {sum(1 for r in sel if r.get('label')=='known_active')}, "
          f"decoy {sum(1 for r in sel if r.get('label')=='decoy')})")

    posedir = outdir / "poses"
    posedir.mkdir(exist_ok=True)
    workers = args.workers or (os.cpu_count() or 2)
    jobs = [(str(rigid), str(flex), pocket["center"], pocket["size"], r, args.exhaustiveness,
             str(posedir) if _group(r) == "benchmark" else None)
            for r in sel]
    results, t0 = [], time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(dock_one, j) for j in jobs]
        for i, f in enumerate(as_completed(futs), 1):
            results.append(f.result())
            if i % 10 == 0 or i == len(futs):
                el = time.time() - t0
                print(f"    {i}/{len(futs)}  {el/60:.1f}분", flush=True)

    rigid_by_id = {r["chembl_id"]: r["score"] for r in rows}
    ok = []
    for r in results:
        if r.get("flex_score") is None:
            continue
        r["rigid_score"] = rigid_by_id.get(r["chembl_id"])
        if r["rigid_score"] is not None:
            r["delta"] = round(r["flex_score"] - r["rigid_score"], 3)
            ok.append(r)

    (outdir / "flex_results.jsonl").write_text(
        "\n".join(json.dumps(r, ensure_ascii=False) for r in ok)
    )
    analysis = analyse(ok)
    (outdir / "flex_analysis.json").write_text(json.dumps(analysis, indent=2, ensure_ascii=False))
    report(analysis)

    wall = (time.time() - t0) / 60
    # GitHub Actions 는 잡의 벽시계 시간으로 과금한다. 코어 수를 곱하면 안 된다
    # (초기 구현의 오류). 리지드 기준선도 같은 단위다: 8샤드 합계 1,223 러너분 /
    # 2,157 리간드 = 0.567 러너분/리간드.
    per_lig = wall / max(1, len(ok))
    rigid_per_lig = 1223 / 2157
    m = per_lig / rigid_per_lig
    full_lib_min = round(2157 * per_lig)
    stats = {
        "gene": args.gene,
        "flex_residues": args.flex_res,
        "n_docked": len(ok),
        "n_failed": len(results) - len(ok),
        "workers": workers,
        "wall_minutes": round(wall, 1),
        "runner_min_per_ligand": round(per_lig, 4),
        "rigid_min_per_ligand": round(rigid_per_lig, 4),
        "time_multiplier_M": round(m, 2),
        "full_library_runner_minutes": full_lib_min,
        "full_library_billable_minutes": max(0, full_lib_min - 2000),
        "full_library_cost_usd": round(max(0, full_lib_min - 2000) * 0.008, 2),
    }
    (outdir / "flex_stats.json").write_text(json.dumps(stats, indent=2, ensure_ascii=False))

    print()
    print(f"[+] 성공 {len(ok)}/{len(results)}")
    print(f"    소요 {stats['wall_minutes']}분 (워커 {workers}개)")
    print(f"    리간드당 {stats['runner_min_per_ligand']} 러너분 "
          f"vs 리지드 {stats['rigid_min_per_ligand']}")
    print(f"    소요 배수 M = {stats['time_multiplier_M']}")
    print(f"    전체 라이브러리 환산 {stats['full_library_runner_minutes']:,} 러너분 "
          f"→ 과금 {stats['full_library_billable_minutes']:,}분 "
          f"≈ ${stats['full_library_cost_usd']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
