#!/usr/bin/env python3
"""
선정된 타깃의 **평가셋과 스크리닝셋**을 처음부터 갈라서 만든다.

## 이 스크립트가 막는 것

세 가지가 실제로 우리를 물었고, 전부 라이브러리 단계에서 생겼다.

1. **교란 2 — 후보/decoy 순환.** 채점에서 음성으로 센 화합물을 그대로 후보로
   내놓고 있었다. 그래서 **평가셋과 스크리닝셋을 한 번에, 서로 겹치지 않게**
   만든다. 나중에 가르지 않는다.
2. **11-01 `rotb` 누락.** 고정본 decoy 2,100개에 `rotb` 가 없어 전부 0으로
   읽혔고, `dock.py` 의 부하 분산이 조용히 실패해 무거운 리간드가 한 샤드에
   몰렸다. 여기서 **모든 행에 `rotb` 를 계산해 넣는다.**
3. **교란 5 — 기억이 도킹을 이긴다.** 비활성이 추정이면 Tanimoto 상한을 걸게
   되고 그 상한이 간격을 만든다. 그래서 **실측 비활성만 쓴다** (C3).

## PAINS 는 표시만 하고 버리지 않는다

평가셋에서 PAINS 를 버리면 측정된 사실을 우리가 고른 것으로 바꾸는 셈이다.
벤치마크는 측정된 그대로여야 한다. 플래그만 달아 두고, 후보를 고를 때
(스크리닝셋 쪽에서) 쓴다.

## 활성 비율을 왜 2.8% 로 두는가

`gate.py` 의 문턱이 활성 비율 2~3% 를 전제한다. EPAS1 평가셋이 2.78% 였으므로
그 근처로 두면 EF 천장이 비슷해져 이전 숫자와 나란히 놓을 수 있다. 비율을
바꾸려면 게이트 문턱도 함께 봐야 한다.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import random
from pathlib import Path


def canon_rows(pairs: list[tuple[str, str, int]]) -> tuple[list[dict], dict]:
    """정규화 SMILES 로 중복을 없앤다. 라벨이 충돌하면 활성을 남기고 센다."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
    RDLogger.DisableLog("rdApp.*")

    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    pains = FilterCatalog(params)

    seen: dict[str, dict] = {}
    bad = conflict = dup = 0
    for cid, smi, y in pairs:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            bad += 1
            continue
        key = Chem.MolToSmiles(m)
        if key in seen:
            dup += 1
            if seen[key]["label"] != ("known_active" if y == 1 else "decoy"):
                conflict += 1
                if y == 1:                      # 한 번이라도 활성이면 활성으로 둔다
                    seen[key]["label"] = "known_active"
                    seen[key]["label_conflict"] = True
                else:
                    seen[key]["label_conflict"] = True
            continue
        seen[key] = {
            "smiles": key,
            "chembl_id": f"TDC{cid}",       # 파이프라인이 이 키로 화합물을 가리킨다
            "pref_name": "",
            "label": "known_active" if y == 1 else "decoy",
            "source": "tdc_hts",
            "measured_inactive": y == 0,     # 추정이 아니라 실제로 음성이었다
            "protected": False,
            # 11-01 에 이게 없어서 부하 분산이 조용히 실패했다
            "rotb": Descriptors.NumRotatableBonds(m),
            "mw": round(Descriptors.MolWt(m), 2),
            "pains": bool(pains.HasMatch(m)),
        }
    return list(seen.values()), {"unparseable": bad, "duplicates": dup,
                                 "label_conflicts": conflict}


def write(path: Path, rows: list[dict]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "\n".join(json.dumps(r, ensure_ascii=False) for r in rows) + "\n"
    path.write_text(text)
    return hashlib.sha256(text.encode()).hexdigest()


def main() -> int:
    import warnings
    warnings.filterwarnings("ignore")
    from tdc.single_pred import HTS

    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, help="TDC HTS 이름")
    ap.add_argument("--gene", required=True, help="exp/<GENE>/ 에 쓴다")
    ap.add_argument("--eval-n", type=int, default=5000)
    ap.add_argument("--eval-active-fraction", type=float, default=0.028)
    ap.add_argument("--screen-n", type=int, default=5000)
    ap.add_argument("--data-dir", default="tdc_raw")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--active-max-sim", type=float, default=0.0,
                    help="평가셋 활성끼리 이 Tanimoto 이상으로 닮지 않게 고른다. "
                         "0 이면 끄고 무작위로 뽑는다")
    ap.add_argument("--seed", type=int, default=1023)
    args = ap.parse_args()

    out_dir = Path(args.out_dir or f"exp/{args.gene}")
    print(f"[1] {args.dataset} 내려받기")
    df = HTS(name=args.dataset, path=args.data_dir).get_data()
    smi_col = "Drug" if "Drug" in df.columns else df.columns[1]
    id_col = "Drug_ID" if "Drug_ID" in df.columns else df.columns[0]
    pairs = [(str(r[id_col]), str(r[smi_col]), int(r["Y"])) for _, r in df.iterrows()]
    print(f"    원본 {len(pairs):,}개")

    print("[2] 정규화·중복 제거·rotb·PAINS")
    rows, qc = canon_rows(pairs)
    acts = [r for r in rows if r["label"] == "known_active"]
    decs = [r for r in rows if r["label"] != "known_active"]
    print(f"    고유 {len(rows):,}개 — 활성 {len(acts):,}, 실측 비활성 {len(decs):,}")
    print(f"    버린 것: 파싱 실패 {qc['unparseable']:,}, 중복 {qc['duplicates']:,}, "
          f"라벨 충돌 {qc['label_conflicts']:,}")
    print(f"    PAINS 표시 {sum(1 for r in rows if r['pains']):,}개 (버리지 않는다)")

    print("[3] 평가셋·스크리닝셋을 겹치지 않게 가른다")
    rng = random.Random(args.seed)
    rng.shuffle(acts)
    rng.shuffle(decs)

    want_act = min(round(args.eval_n * args.eval_active_fraction), len(acts))
    dropped_sim = 0
    if args.active_max_sim > 0:
        # 활성끼리 닮은 쌍이 평가셋에 같이 들어가면, 단백질 없이 "활성과
        # 닮았나"만 세는 기억 베이스라인이 그 쌍을 서로 찾아낸다. 11-16 에
        # 활성 140개로 재니 바로 그 이유로 게이트를 3/3 통과했다.
        # 탈락한 활성은 버리지 않고 스크리닝셋으로 보낸다 — 거기서는
        # 숨은 양성으로 쓸모가 있다.
        from rdkit import Chem, DataStructs, RDLogger
        from rdkit.Chem import rdFingerprintGenerator
        RDLogger.DisableLog("rdApp.*")
        fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=4096)
        picked, picked_fp, rest = [], [], []
        for r in acts:
            fp = fpgen.GetFingerprint(Chem.MolFromSmiles(r["smiles"]))
            too_close = bool(picked_fp) and max(
                DataStructs.BulkTanimotoSimilarity(fp, picked_fp)) >= args.active_max_sim
            if too_close:
                dropped_sim += 1          # 닮아서 뺀 것만 센다
                rest.append(r)
            elif len(picked) < want_act:
                picked.append(r)
                picked_fp.append(fp)
            else:
                rest.append(r)            # 자리가 차서 남은 것 — 닮아서가 아니다
        print(f"    활성 다양성 필터 (Tanimoto < {args.active_max_sim}): "
              f"{len(picked)}개 선택, 닮아서 뺀 것 {dropped_sim}개, "
              f"자리가 차서 남은 것 {len(rest) - dropped_sim}개")
        acts = picked + rest          # 앞쪽 picked 만 평가셋에 들어간다
    n_eval_act = min(want_act, len(acts))
    n_eval_dec = min(args.eval_n - n_eval_act, len(decs))
    ev = acts[:n_eval_act] + decs[:n_eval_dec]

    rest_act = acts[n_eval_act:]
    rest_dec = decs[n_eval_dec:]
    n_scr_act = min(len(rest_act), args.screen_n)
    n_scr_dec = min(args.screen_n - n_scr_act, len(rest_dec))
    sc = rest_act[:n_scr_act] + rest_dec[:n_scr_dec]

    rng.shuffle(ev)
    rng.shuffle(sc)

    # 겹치면 교란 2 가 그대로 돌아온다. 가정하지 말고 확인한다.
    overlap = {r["smiles"] for r in ev} & {r["smiles"] for r in sc}
    if overlap:
        print(f"    [!] 평가셋과 스크리닝셋이 {len(overlap)}개 겹친다 — 중단한다")
        return 1
    print(f"    겹침 0개 확인")

    h_ev = write(out_dir / "library_eval.jsonl", ev)
    h_sc = write(out_dir / "library_screen.jsonl", sc)

    stats = {
        "gene": args.gene, "dataset": args.dataset, "seed": args.seed,
        "source_rows": len(pairs), "unique_rows": len(rows), "qc": qc,
        "pains_flagged": sum(1 for r in rows if r["pains"]),
        "active_max_sim": args.active_max_sim,
        "actives_moved_to_screen_by_similarity": dropped_sim,
        "eval": {"n_total": len(ev), "n_active": n_eval_act,
                 "active_fraction": round(n_eval_act / len(ev), 5),
                 "sha256": h_ev, "file": "library_eval.jsonl"},
        "screen": {"n_total": len(sc), "n_active": n_scr_act,
                   "sha256": h_sc, "file": "library_screen.jsonl"},
        "disjoint_verified": True,
    }
    (out_dir / "library_stats.json").write_text(
        json.dumps(stats, ensure_ascii=False, indent=2))

    print(f"\n[4] 산출물 — {out_dir}/")
    print(f"    평가셋      {len(ev):>6,}개  활성 {n_eval_act:>4}  "
          f"({100 * n_eval_act / len(ev):.2f}%)  sha256 {h_ev[:16]}…")
    print(f"    스크리닝셋  {len(sc):>6,}개  활성 {n_scr_act:>4}  "
          f"(남은 활성 전부)      sha256 {h_sc[:16]}…")
    print(f"    rotb 중앙값 {sorted(r['rotb'] for r in ev)[len(ev) // 2]} "
          f"(도킹 부하 분산에 쓴다)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
