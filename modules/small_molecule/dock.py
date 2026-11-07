#!/usr/bin/env python3
"""
modules.small_molecule.dock — AutoDock Vina 도킹.

리간드: SMILES → ETKDGv3 3D 임베딩 → MMFF 최적화 → meeko PDBQT
수용체: obabel로 PDBQT 변환 (강체, Gasteiger 전하, pH 7.4 수소 부가)
박스:   core.structure.fetch가 정한 center/size

GitHub Actions 6시간 제한을 고려해 --n-shards로 매트릭스 병렬화한다.

알려진 한계 (보고서에 명시할 것)
  - 리지드 수용체. 측쇄 유연성과 유도 적합(induced fit)을 반영하지 않는다.
  - Vina 스코어링은 금속 배위와 PPI 계면을 과소평가한다.
  - 리간드 프로토네이션은 ChEMBL 표기 그대로 — pH 7.4 상태를 계산하지 않는다.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

# 의존성은 모듈 로드 시점에 터뜨린다. 워커 안에서 조용히 실패하면
# 리간드 준비가 전부 실패해도 '도킹은 돌았다'처럼 보인다 (실제로 겪은 사고).
from meeko import MoleculePreparation, PDBQTWriterLegacy  # noqa: E402
from vina import Vina  # noqa: E402

RDLogger.DisableLog("rdApp.*")


# AutoDock 4 원자 타입 — Vina 가 수용체에서 받아들이는 것들.
# 이 집합 밖의 타입이 PDBQT 에 들어가면 Vina 가 파이썬 예외가 아니라 네이티브에서
# 죽는다. 스택 트레이스도 없고, 버퍼에 쌓인 출력도 함께 사라진다.
AD4_RECEPTOR_TYPES = {
    "H", "HD", "HS", "C", "A", "N", "NA", "NS", "OA", "OS", "F",
    "MG", "P", "SA", "S", "CL", "CA", "MN", "FE", "ZN", "BR", "I", "SI", "B",
}


def sanitize_receptor_pdbqt(pdbqt: Path) -> dict:
    """수용체 PDBQT 에서 Vina 가 모르는 원자 타입을 걸러낸다.

    왜 필요한가
      2026-10-19 MET 대조군(실행 39186266215)에서 8샤드가 전부 1~3초 만에
      죽었고 로그는 0바이트였다. 파이썬 예외였다면 트레이스백이 남았을 것이다.
      출력이 통째로 사라진 것은 버퍼를 비우기 전에 시그널로 죽었다는 뜻 —
      네이티브 크래시다.

      treatment(5HTI)와 control(3DKC)의 파이프라인 경로에서 다른 것은 보조인자
      하나뿐이었다. 3DKC 는 Mg2+ 와 Cl- 를 갖고 있고, fetch.py 의 화이트리스트가
      이들을 수용체에 보존한다. Open Babel 이 이온에 매긴 타입이 Vina 의 파서를
      통과하지 못한 것으로 보인다.

      어느 타입이 범인인지 단정하지 않고, 알려진 집합 밖은 전부 떨어뜨린다.
      떨어뜨린 것은 반드시 로그에 남긴다 — 조용히 지우면 다음 사람이 같은
      자리에서 또 헤맨다.

    주의
      금속 보조인자는 생물학적으로 결합 부위의 일부일 수 있다. 키나제 ATP
      부위의 Mg2+ 가 그렇다. 떨어뜨린다는 것은 그 상호작용을 모델링하지
      않겠다는 선언이므로, 결과 해석에 반영해야 한다.
    """
    kept, dropped = [], []
    hist: dict[str, int] = {}
    for line in pdbqt.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            kept.append(line)
            continue
        t = line[77:79].strip().upper() if len(line) >= 78 else ""
        hist[t or "?"] = hist.get(t or "?", 0) + 1
        if t in AD4_RECEPTOR_TYPES:
            kept.append(line)
        else:
            dropped.append((t or "?", line[17:20].strip(), line[22:26].strip()))

    if dropped:
        pdbqt.write_text("\n".join(kept) + "\n")
    return {
        "atom_types": dict(sorted(hist.items(), key=lambda kv: -kv[1])),
        "n_dropped": len(dropped),
        "dropped": sorted({f"{t} ({res}{num})" for t, res, num in dropped}),
    }


MIN_RECEPTOR_ATOMS = 200

# 단원자 이온 — Open Babel 의 Gasteiger 전하 계산이 이들에서 실패한다.
MONATOMIC_IONS = {
    "MG", "ZN", "MN", "FE", "FE2", "CA", "CU", "CU1", "NI", "CO",
    "K", "NA", "CL", "BR", "IOD", "F", "CD", "HG", "SR", "BA", "CS", "RB", "LI",
}

# 촉매 금속 — 지우면 결합 부위 자체가 사라지는 것들.
#
# 2026-10-25 발견. 이온 제거는 MET 대조군의 obabel 실패를 고치려고 넣었는데
# (가스테이거 전하 부여가 고립 단원자 이온에서 깨진다), 종류를 가리지 않는다.
# config/targets.yaml 의 CA9(탄산무수화효소 IX)는 pocket_mode=metal_site 이고
# 활성부위 Zn2+ 에 설폰아마이드 억제제가 직접 배위한다. 그 아연을 지우면
# 도킹할 부위가 없어진다 — 11월에 CA9 로 넓혔다면 조용히 무의미한 결과가
# 나왔을 것이다.
#
# 그래서 촉매 금속은 기본적으로 남긴다. 남기면 obabel 이 깨질 수 있으므로
# prep_receptor 가 전환에 실패하면 그때 해당 이온을 빼고 한 번 더 시도하되,
# 뺐다는 사실을 메타데이터에 남겨 절대 해석에서 빠지지 않게 한다.
CATALYTIC_METALS = {"ZN", "MN", "FE", "FE2", "CU", "CU1", "NI", "CO", "MG"}


def strip_monatomic_ions(
    pdb: Path, keep_catalytic: bool = True
) -> tuple[Path, list[str]]:
    """Gasteiger 전하 계산을 막는 단원자 이온을 뺀 사본을 만든다.

    keep_catalytic=True 면 촉매 금속(Zn, Mn, Fe, Cu, Ni, Co, Mg)은 남긴다.
    CA9 처럼 금속 배위가 결합의 본질인 타깃에서 그 금속을 지우면 도킹할
    부위 자체가 사라지기 때문이다. 남겨서 obabel 이 깨지면 호출부가
    keep_catalytic=False 로 한 번 더 시도하고, 그 사실을 기록한다.

    왜
      2026-10-19 MET 대조군(3DKC)에서 수용체 변환이 종료코드 0 으로 끝나면서
      빈 파일을 냈다. 변형 시험 결과:

        -xr -p 7.4 --partialcharge gasteiger   → 0 molecules converted
        -xr                                     → 1 molecule, 2,453 원자
        -xr --partialcharge gasteiger           → 0 molecules converted
        이온 제거 + 위 전체 옵션                  → 1 molecule, 2,960 원자

      -p 7.4 가 아니라 --partialcharge gasteiger 가 범인이다. 고립된 Mg2+/Cl-
      에서 전하 할당이 실패한다. 3DKC 의 HETATM 은 그 둘뿐이었다.

    왜 폴백이 아니라 이온 제거인가
      옵션을 줄여 폴백하면 두 팔의 수용체 준비 방식이 달라진다. 이 실험의
      전제는 '구조만 다르고 나머지는 동일' 이므로, 양성자화 여부가 달라지면
      측정한 차이가 구조 선택 때문인지 준비 방식 때문인지 구분할 수 없게 된다.
      편향 측정이라는 목적 자체가 무너진다.

      이온 제거는 이온이 없는 구조(MET treatment 5HTI 는 보조인자가 비어 있다)
      에서는 아무 일도 하지 않는다. 두 팔이 같은 명령으로 준비된다.

    해석에 반영할 것
      금속 보조인자는 생물학적으로 결합 부위의 일부일 수 있다. 키나제 ATP
      부위의 Mg2+ 가 그렇다. 이것을 빼면 Mg 배위는 모델링되지 않는다.
      다만 Vina 의 스코어 함수는 금속 배위를 제대로 다루지 못하고, 비교
      대상인 두 팔이 똑같이 Mg 없이 준비되므로 비교 자체는 성립한다.
      절대값 해석에만 영향이 있다.
    """
    drop = MONATOMIC_IONS - (CATALYTIC_METALS if keep_catalytic else set())
    kept, removed = [], []
    for line in pdb.read_text().splitlines():
        if line.startswith("HETATM") and line[17:20].strip().upper() in drop:
            removed.append(f"{line[17:20].strip()}{line[22:26].strip()}")
            continue
        kept.append(line)
    if not removed:
        return pdb, []
    suffix = "_noion" if keep_catalytic else "_nometal"
    out = pdb.with_name(pdb.stem + suffix + ".pdb")
    out.write_text("\n".join(kept) + "\n")
    return out, sorted(set(removed))


def _record_metal_stripped(pdbqt: Path, removed: list[str]) -> None:
    """촉매 금속을 뺐다는 사실을 수용체 옆에 남긴다. 보고서가 읽어 간다."""
    try:
        (pdbqt.parent / "metal_stripped.json").write_text(
            json.dumps({"metal_stripped": True, "removed": removed,
                        "note": "촉매 금속을 제거하고 변환했다. 금속 배위는 "
                                "모델링되지 않으므로 절대값 해석에 반영할 것."},
                       ensure_ascii=False, indent=2))
    except Exception:
        pass


def prep_receptor(pdb: Path, pdbqt: Path) -> bool:
    """수용체 PDB → PDBQT.

    조용히 실패하지 않는다
      2026-10-19 MET 대조군에서 8샤드가 1~3초 만에 죽고 로그가 0바이트였다.
      원인은 네이티브 크래시가 아니라 침묵 세 겹이었다.
        1. obabel 이 종료코드 0 으로 끝나면서 빈 파일을 썼다
        2. 이 함수가 그것을 메시지 없이 False 로 반환했다 (예외 경로에서만 출력)
        3. dock.py 가 메시지 없이 return 1 했다
      obabel 의 stdout/stderr 를 삼키고 있었던 것이 결정적이었다. 이제 남긴다.
    """
    if pdbqt.exists() and pdbqt.stat().st_size > 0:
        return True

    # 1차: 촉매 금속을 남기고 시도. CA9 처럼 금속 배위가 결합의 본질인
    # 타깃에서 아연을 지우면 도킹할 부위가 사라진다.
    # 2차: 그래도 obabel 이 깨지면 금속까지 빼고 한 번 더. 뺐다는 사실은
    # metal_stripped 로 기록해 절대값 해석에서 빠지지 않게 한다.
    r = None
    src = pdb
    for keep_catalytic in (True, False):
        src, removed = strip_monatomic_ions(pdb, keep_catalytic=keep_catalytic)
        if removed:
            print(f"    단원자 이온 {len(removed)}개 제외: {', '.join(removed)}",
                  flush=True)
        if not keep_catalytic:
            print("    [!] 촉매 금속까지 제외하고 재시도한다 — 금속 배위가 "
                  "모델링되지 않는다. 절대값 해석에 반드시 반영할 것.",
                  file=sys.stderr, flush=True)
            _record_metal_stripped(pdbqt, removed)

        cmd = [
            "obabel", str(src), "-O", str(pdbqt),
            "-xr", "-p", "7.4", "--partialcharge", "gasteiger",
        ]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            print(f"[!] 수용체 변환 실패: {e}", file=sys.stderr, flush=True)
            return False

        size = pdbqt.stat().st_size if pdbqt.exists() else 0
        if r.returncode == 0 and size > 0:
            break
        if keep_catalytic:
            print("[!] 촉매 금속을 남긴 채로는 변환이 실패했다 — 금속을 빼고 재시도",
                  file=sys.stderr, flush=True)
            continue

        print(f"[!] 수용체 변환 실패 — obabel exit={r.returncode}, 출력 {size} 바이트",
              file=sys.stderr, flush=True)
        for stream, name in ((r.stderr, "stderr"), (r.stdout, "stdout")):
            for line in (stream or "").strip().splitlines()[-25:]:
                print(f"    obabel {name}: {line}", file=sys.stderr, flush=True)
        print(f"    입력: {src} ({src.stat().st_size if src.exists() else 0} 바이트)",
              file=sys.stderr, flush=True)
        return False

    # 빈 파일이 아니어도 원자가 없으면 의미가 없다. Vina 는 원자 없는 PDBQT 를
    # 오류 없이 받아들이므로(위 진단에서 set_receptor OK 였다) 여기서 막는다.
    n = sum(1 for l in pdbqt.read_text().splitlines()
            if l.startswith(("ATOM", "HETATM")))
    if n < MIN_RECEPTOR_ATOMS:
        print(f"[!] 수용체 PDBQT 의 원자가 {n}개뿐이다 (최소 {MIN_RECEPTOR_ATOMS}). "
              "변환이 반쯤 된 파일로 도킹하면 결과가 전부 무의미하다.",
              file=sys.stderr, flush=True)
        return False
    print(f"    수용체 PDBQT 원자 {n:,}개", flush=True)
    return True


def ligand_pdbqt(smiles: str, seed: int = 42) -> str | None:
    """SMILES → 3D → PDBQT 문자열."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        params.useSmallRingTorsions = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            params.useRandomCoords = True
            if AllChem.EmbedMolecule(mol, params) != 0:
                return None
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=600)
        except Exception:
            pass
        prep = MoleculePreparation()
        setups = prep.prepare(mol)
        if not setups:
            return None
        pdbqt, ok, _err = PDBQTWriterLegacy.write_string(setups[0])
        return pdbqt if ok else None
    except Exception:
        return None


def dock_one(job: tuple) -> dict:
    rec_path, center, size, rec, exhaustiveness = job
    out = {k: rec[k] for k in rec if k != "protected"}
    try:
        lig = ligand_pdbqt(rec["smiles"])
        if lig is None:
            out.update(score=None, error="ligand_prep_failed")
            return out

        v = Vina(sf_name="vina", cpu=1, seed=42, verbosity=0)
        v.set_receptor(rec_path)
        v.set_ligand_from_string(lig)
        v.compute_vina_maps(center=list(center), box_size=list(size))
        v.dock(exhaustiveness=exhaustiveness, n_poses=5)
        energies = v.energies(n_poses=5)
        best = float(energies[0][0])
        out.update(
            score=best,
            pose_spread=round(float(energies[-1][0]) - best, 3),
            ligand_efficiency=None,
        )
        mol = Chem.MolFromSmiles(rec["smiles"])
        if mol is not None and mol.GetNumHeavyAtoms():
            out["ligand_efficiency"] = round(-best / mol.GetNumHeavyAtoms(), 4)
        return out
    except Exception as e:  # noqa: BLE001
        out.update(score=None, error=f"{type(e).__name__}: {e}"[:200])
        return out


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--targetdir", default=None)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--max-ligands", type=int, default=0, help="샤드당 상한, 0=전체")
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--workers", type=int, default=0)
    args = ap.parse_args()

    tdir = Path(args.targetdir) if args.targetdir else root / "targets" / args.gene
    meta = json.loads((tdir / "target_meta.json").read_text())
    pocket = meta["pocket"]

    if pocket["confidence"] == "low":
        print(f"[!] {args.gene}: blind docking (포켓 신뢰도 낮음). 결과 해석 주의.")

    rec_pdbqt = tdir / "receptor.pdbqt"
    if not prep_receptor(tdir / "receptor.pdb", rec_pdbqt):
        print("[!] 수용체를 준비하지 못해 중단한다.", file=sys.stderr, flush=True)
        return 1

    # 수용체 원자 타입 점검 — Vina 가 모르는 타입은 네이티브 크래시를 낸다.
    san = sanitize_receptor_pdbqt(rec_pdbqt)
    print(f"    수용체 원자 타입: {san['atom_types']}", flush=True)
    if san["n_dropped"]:
        print(f"    [!] Vina 가 모르는 원자 {san['n_dropped']}개를 제외했다: "
              f"{', '.join(san['dropped'])}", flush=True)
        print("        금속 보조인자라면 그 상호작용은 모델링되지 않는다 — "
              "결과 해석에 반영할 것.", flush=True)

    # 사전 점검 — 리간드 하나를 실제로 준비해 본다. 여기서 막히면 전체가 막히므로
    # 수백 개를 돌려 100% 실패를 확인하느라 한 시간을 쓰지 않는다.
    probe = ligand_pdbqt("CC(=O)Oc1ccccc1C(=O)O")
    if probe is None:
        print(
            "[!] 사전 점검 실패: 리간드 PDBQT 변환이 되지 않는다. "
            "meeko/rdkit 의존성(scipy 포함)을 확인할 것.",
            file=sys.stderr,
        )
        return 1
    print(f"    사전 점검 OK (프로브 리간드 {len(probe.splitlines())}줄)", flush=True)

    # 수용체 사전 점검 — 지금까지 리간드만 프로브하고 수용체는 그대로 믿었다.
    # MET 대조군이 여기서 죽었다. 격자 계산까지 실제로 해 본다.
    print("    수용체 사전 점검: set_receptor + compute_vina_maps ...", flush=True)
    try:
        _v = Vina(sf_name="vina", cpu=1, seed=42, verbosity=0)
        _v.set_receptor(str(rec_pdbqt))
        _v.set_ligand_from_string(probe)
        _v.compute_vina_maps(center=list(pocket["center"]), box_size=list(pocket["size"]))
    except Exception as e:  # noqa: BLE001
        print(f"[!] 수용체 사전 점검 실패: {type(e).__name__}: {e}", file=sys.stderr, flush=True)
        return 1
    print("    수용체 사전 점검 OK", flush=True)

    lib = [
        json.loads(l) for l in (tdir / "library.jsonl").read_text().splitlines() if l.strip()
    ]
    # 회전결합이 많은 리간드일수록 Vina 탐색이 오래 걸린다. 인덱스로 그냥
    # 자르면 무거운 리간드가 한 샤드에 몰려 그 샤드만 타임아웃된다
    # (EPAS1 1차 실행에서 4개 중 3개가 350분 타임아웃으로 취소됐다).
    # 복잡도 순으로 정렬한 뒤 라운드로빈으로 나눠 부하를 고르게 만든다.
    lib.sort(key=lambda r: -(r.get("rotb") or 0))
    shard_lib = [r for i, r in enumerate(lib) if i % args.n_shards == args.shard]
    if args.max_ligands:
        shard_lib = shard_lib[: args.max_ligands]

    workers = args.workers or (os.cpu_count() or 2)
    n_act = sum(1 for r in shard_lib if r["label"] == "known_active")
    print(
        f"[*] {args.gene} 도킹: {len(shard_lib)} 리간드 (활성 {n_act}) "
        f"샤드 {args.shard + 1}/{args.n_shards}, 워커 {workers}"
    )
    print(f"    박스 center={pocket['center']} size={pocket['size']} ({pocket['mode']})")

    jobs = [
        (str(rec_pdbqt), pocket["center"], pocket["size"], r, args.exhaustiveness)
        for r in shard_lib
    ]
    # 한 리간드가 끝날 때마다 바로 쓴다.
    #
    # 왜: 예전에는 샤드의 리간드를 **전부** 끝낸 뒤에야 출력을 썼다. 그래서
    # 취소된 샤드는 산출물이 0 이었고, 2026-10 과 2026-11 두 달 모두 전체
    # 러너 소비의 46% 가 그렇게 날아갔다 (10월 2,554분, 11월 1,194분).
    #
    # 이름을 `partial_shard*` 로 둔 이유: 하류 네 곳이
    # `docked_shard*.jsonl` 글롭으로 샤드를 모은다. 부분 파일이 그 글롭에
    # 걸리면 불완전하고 정렬도 안 된 데이터가 **조용히** 평가셋에 섞인다 —
    # 이 프로젝트가 잡아내려고 만들어진 바로 그 종류의 실패다.
    #
    # 완주하면 정렬된 docked_shard 를 쓰고 부분 파일을 지운다. 그러므로
    # **partial_shard 파일이 남아 있다 = 그 샤드는 잘렸다** 는 신호다.
    partial = tdir / f"partial_shard{args.shard}.jsonl"
    partial.unlink(missing_ok=True)

    results, t0 = [], time.time()
    pf = partial.open("w")
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = [ex.submit(dock_one, j) for j in jobs]
        for i, f in enumerate(as_completed(futures), 1):
            r = f.result()
            results.append(r)
            pf.write(json.dumps(r, ensure_ascii=False) + "\n")
            pf.flush()  # 버퍼에 남으면 취소 시 그대로 사라진다
            if i % 100 == 0 or i == len(futures):
                el = time.time() - t0
                rate = i / el if el else 0
                eta = (len(futures) - i) / rate if rate else 0
                print(
                    f"    {i}/{len(futures)}  {el / 60:.1f}분 경과, "
                    f"잔여 ~{eta / 60:.1f}분",
                    flush=True,
                )

    pf.close()

    scored = [r for r in results if r.get("score") is not None]
    scored.sort(key=lambda r: r["score"])
    out = tdir / f"docked_shard{args.shard}.jsonl"
    out.write_text("\n".join(json.dumps(r, ensure_ascii=False) for r in scored) + "\n")
    # 완주했으므로 부분 파일을 지운다 — 남아 있으면 '잘린 샤드' 신호다
    partial.unlink(missing_ok=True)

    errs: dict[str, int] = {}
    for r in results:
        if r.get("score") is None:
            e = r.get("error", "unknown").split(":")[0]
            errs[e] = errs.get(e, 0) + 1

    print(f"[+] 성공 {len(scored)}/{len(results)} → {out.name}")
    if errs:
        print(f"    실패 내역: {errs}")
    if scored:
        print(f"    최고 스코어: {scored[0]['score']:.2f} kcal/mol ({scored[0]['chembl_id']})")
    fail_rate = 1 - len(scored) / max(1, len(results))
    if fail_rate > 0.10:
        print(f"    [!] 실패율 {fail_rate:.0%} — 10% 초과. 리간드 준비 단계 점검 필요.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
