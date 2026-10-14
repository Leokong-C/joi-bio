#!/usr/bin/env python3
"""
core.report.build — 타깃별 결과를 하나의 보고서로 묶는다.

검증 FAIL 또는 INSUFFICIENT 타깃은 후보를 나열하지 않고 실패 원인만 기록한다.
근거 없는 히트 주장을 막는 것이 이 모듈의 역할이다.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
VERDICT_MARK = {"PASS": "통과", "FAIL": "미통과", "INSUFFICIENT": "데이터 부족"}


def load(p: Path):
    try:
        return json.loads(p.read_text()) if p.exists() else None
    except json.JSONDecodeError:
        return None


def target_section(t: dict, tdir: Path) -> str:
    gene = t["gene"]
    title = gene + (f" ({t['alias']})" if t.get("alias") else "")
    lines = [f"## {title}", ""]
    lines.append(f"- 질환: {t['disease']} · 분류: {t['class']}")
    lines.append(f"- 선정 근거: {t['rationale']}")

    meta = load(tdir / "target_meta.json")
    if not meta:
        lines.append("- **상태: 미실행**")
        return "\n".join(lines) + "\n"
    if meta.get("error"):
        lines.append(f"- **상태: 구조 확보 실패** — {meta['error']}")
        return "\n".join(lines) + "\n"

    rc, pk = meta.get("receptor", {}), meta.get("pocket", {})
    src = f"{meta.get('source')} {meta.get('pdb_id') or ''}".strip()
    res = f", {meta['resolution']} Å" if meta.get("resolution") else ""
    lines.append(f"- 구조: {src}{res} · 체인 {rc.get('chain')} · 원자 {rc.get('n_atoms')}")
    lines.append(
        f"- 포켓: {pk.get('mode')} (신뢰도 {pk.get('confidence')}) — {pk.get('note', '')}"
    )
    if rc.get("cofactors"):
        lines.append(f"- 보조인자 유지: {', '.join(rc['cofactors'])}")
    if meta.get("warning"):
        lines.append(f"- 주의: {meta['warning']}")
    if meta.get("declared_risk"):
        lines.append(f"- 사전 식별 리스크: {meta['declared_risk']}")

    val = load(tdir / "validation.json")
    if not val:
        lines.append("- **상태: 검증 미실행**")
        return "\n".join(lines) + "\n"

    fmt = lambda k: val[k] if val.get(k) is not None else "n/a"
    lines.append("")
    lines.append(
        f"**검증: {VERDICT_MARK.get(val['verdict'], val['verdict'])}** — {val['verdict_reason']}"
    )
    lines.append("")
    lines.append("| 지표 | 값 | 임계값 |")
    lines.append("| --- | --- | --- |")
    lines.append(f"| EF1% | {fmt('ef1')} | ≥ {val['thresholds']['ef1']} |")
    lines.append(f"| EF5% | {fmt('ef5')} | — |")
    lines.append(f"| AUC | {fmt('auc')} | ≥ {val['thresholds']['auc']} |")
    lines.append(f"| BEDROC | {fmt('bedroc')} | ≥ {val['thresholds']['bedroc']} |")
    lines.append(
        f"| 화합물 수 | {val['n_compounds']} (활성 {val['n_known_active']}) | 활성 ≥ {val['min_actives_required']} |"
    )

    br = val.get("benchmark_recovery", {})
    if br.get("checked"):
        lines.append("")
        lines.append(f"**벤치마크 회수** — 승인약·임상 화합물이 상위 {br['top_pct']}%에 들어왔는가")
        lines.append("")
        lines.append("| 화합물 | 순위 | 백분위 | Vina (kcal/mol) | 상위 1% |")
        lines.append("| --- | --- | --- | --- | --- |")
        for f in br["checked"]:
            lines.append(
                f"| {f['name']} | {f['rank']} | {f['percentile']}% | "
                f"{f['score']} | {'O' if f['in_top_pct'] else 'X'} |"
            )
    if br.get("missing"):
        lines.append("")
        lines.append(f"라이브러리에 없는 벤치마크: {', '.join(br['missing'])}")

    if val["verdict"] != "PASS":
        lines.append("")
        lines.append("후보를 보고하지 않습니다. 검증을 통과하지 못한 예측은 근거가 없습니다.")
        return "\n".join(lines) + "\n"

    cands = [
        json.loads(l)
        for l in (tdir / "candidates.jsonl").read_text().splitlines()
        if l.strip()
    ][:10] if (tdir / "candidates.jsonl").exists() else []

    if cands:
        lines.append("")
        lines.append("**상위 후보 10개**")
        lines.append("")
        lines.append("| # | ChEMBL ID | Vina | LE | QED | MW | cLogP | TPSA | 신배설 |")
        lines.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for i, c in enumerate(cands, 1):
            a = c["admet"]
            le = c.get("ligand_efficiency")
            lines.append(
                f"| {i} | {c['chembl_id']} | {c['score']:.2f} | "
                f"{le if le is not None else '—'} | {a['qed']:.2f} | "
                f"{a['mw']:.0f} | {a['clogp']:.1f} | {a['tpsa']:.0f} | "
                f"{'O' if a['renal_clearance_favored'] else '—'} |"
            )
        lines.append("")
        lines.append(
            "위 화합물은 in silico 예측 결과이며 실험 검증 전 단계입니다 "
            "(candidate_for_testing). 결합이 확인된 것이 아닙니다."
        )
    return "\n".join(lines) + "\n"


def main() -> int:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--targetroot", default=str(ROOT / "targets"))
    ap.add_argument("--out", default=str(ROOT / "results" / "report.md"))
    args = ap.parse_args()

    troot = Path(args.targetroot)
    targets = yaml.safe_load((ROOT / "config" / "targets.yaml").read_text())["targets"]
    targets.sort(key=lambda t: t.get("priority", 99))

    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    head = [
        "# Joi_bio — 신장 타깃 가상 스크리닝 보고서",
        "",
        f"생성 {ts} · 타깃 {len(targets)}개",
        "",
        "## 요약",
        "",
        "| 타깃 | 질환 | 구조 | 포켓 신뢰도 | 검증 | 후보 |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    n_pass = n_run = 0
    for t in targets:
        tdir = troot / t["gene"]
        meta, val = load(tdir / "target_meta.json"), load(tdir / "validation.json")
        cand_f = tdir / "candidates.jsonl"
        n_cand = (
            len([l for l in cand_f.read_text().splitlines() if l.strip()])
            if cand_f.exists()
            else 0
        )
        if val:
            n_run += 1
            if val["verdict"] == "PASS":
                n_pass += 1
        src = (meta or {}).get("source", "—")
        conf = ((meta or {}).get("pocket") or {}).get("confidence", "—")
        verd = VERDICT_MARK.get((val or {}).get("verdict", ""), "—")
        head.append(
            f"| {t['gene']} | {t['disease']} | {src} | {conf} | {verd} | "
            f"{n_cand if verd == '통과' else '—'} |"
        )

    head += [
        "",
        f"**실행 {n_run}/{len(targets)}, 검증 통과 {n_pass}.** "
        "검증을 통과한 타깃의 후보만 실험 검증 대상으로 제안합니다.",
        "",
        "### 방법의 한계",
        "",
        "- 리지드 수용체 도킹입니다. 측쇄 유연성과 유도 적합을 반영하지 않습니다.",
        "- Vina 스코어링 함수는 금속 배위(CA9)와 단백질-단백질 계면(KEAP1, VHL)을 "
        "과소평가하는 경향이 있습니다. 이들 타깃의 낮은 농축도는 파이프라인 결함이 "
        "아니라 방법의 알려진 한계일 수 있습니다.",
        "- 리간드 프로토네이션 상태는 ChEMBL 표기를 그대로 사용했으며 pH 7.4 기준으로 "
        "재계산하지 않았습니다.",
        "- ADMET 값은 경험적 rule 기반이며 학습된 예측 모델이 아닙니다.",
        "",
    ]

    body = "\n".join(target_section(t, troot / t["gene"]) for t in targets)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(head) + "\n" + body)
    print(f"[+] {out} 생성 (검증 통과 {n_pass}/{n_run} 실행)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
