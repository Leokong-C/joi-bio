#!/usr/bin/env python3
"""게이트 자체 점검 — 알려진 정답 항목으로 판정 로직을 매번 검증한다.

왜 있는가. 구조 탐침이 두 번 조용히 틀린 수를 냈다(나트륨 이온을 리간드로
세어 부풀리고, 필드를 잘못 읽어 0으로 깎았다). 둘 다 표만 보면 그럴듯했다.
게이트는 그보다 더 조용히 틀릴 수 있다 — THRESHOLDS 의 키를 바꾸면
`m.get("nef1")` 이 없는 하류 스크립트에서 `or 0` 에 걸려 "기억이 게이트를
통과하지 못한다" 는 거짓 안심을 낸다. 그 경로를 여기서 막는다.

실패하면 0이 아닌 종료 코드를 낸다. 워크플로에서 게이트 앞에 돌린다.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from core.validation.gate import (  # noqa: E402
    ACTIVE_FRACTION_MAX,
    THRESHOLDS,
    evaluate,
    normalized_ef,
    perfect_ef,
)

FAILURES: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    mark = "OK  " if ok else "실패"
    print(f"  [{mark}] {name}" + (f" — {detail}" if detail else ""))
    if not ok:
        FAILURES.append(f"{name}: {detail}")


def rows_with_actives_on_top(n: int, n_act: int, on_top: int) -> list[dict]:
    """활성 on_top 개를 최상위에, 나머지 활성을 최하위에 둔 합성 결과.

    Vina 는 낮을수록 좋으므로 score 를 증가 순으로 준다.
    """
    out = []
    n_bottom_act = n_act - on_top
    for i in range(n):
        if i < on_top:
            label = "known_active"
        elif i >= n - n_bottom_act:
            label = "known_active"
        else:
            label = "decoy"
        out.append({"score": -20.0 + i * 0.01, "label": label})
    return out


print("[1] 문턱 정의 자체")
check(
    "THRESHOLDS 가 nef1 을 쓴다 (절대 EF1% 가 아니다)",
    "nef1" in THRESHOLDS and "ef1" not in THRESHOLDS,
    f"{THRESHOLDS}",
)
check(
    "nEF1% 문턱이 0~1 범위에 있다",
    0.0 < THRESHOLDS["nef1"] < 1.0,
    f"nef1={THRESHOLDS['nef1']}",
)

print("\n[2] 하류 스크립트가 게이트 키를 전부 낸다 — '조용한 0' 방지")
# THRESHOLDS 의 키 하나라도 빠진 지표 사전에 `or 0` 을 쓰면 미달로 조용히
# 깎인다. 각 스크립트가 실제로 내는 키를 코드에서 직접 확인한다.
for mod, fn_src in (
    ("scripts/trivial_baselines.py", "rank_metrics"),
    ("scripts/screen_targets.py", "memorisation"),
    ("scripts/leakage_audit.py", "nef1"),
):
    src = (Path(__file__).resolve().parents[1] / mod).read_text()
    missing = [k for k in THRESHOLDS if f'"{k}"' not in src]
    check(f"{mod} 가 {list(THRESHOLDS)} 를 모두 낸다", not missing, f"빠진 키 {missing}")

print("\n[3] nEF1% 계산")
top = perfect_ef(2000, 60, 0.01)
check("천장이 계산된다", top is not None and top > 0, f"perfect_ef(2000,60)={top}")
check(
    "천장에 닿으면 nEF1% = 1.0",
    normalized_ef(top, 2000, 60, 0.01) == 1.0,
    f"{normalized_ef(top, 2000, 60, 0.01)}",
)
check("EF1% 0 이면 nEF1% 0", normalized_ef(0.0, 2000, 60, 0.01) == 0.0)

print("\n[4] 판정 — 치료군")
# 활성 전부 최상위: nEF1% 1.0. 대조군이 없으므로 ABSTAIN 이어야 한다
# (문턱 통과만으로 PASS 를 주지 않는다 — 규칙 1).
r = evaluate(rows_with_actives_on_top(2000, 60, 60), [])
check(
    "문턱을 다 넘어도 대조군이 없으면 ABSTAIN",
    r["verdict"] == "ABSTAIN",
    f"{r['verdict']}",
)
check("nEF1% 가 판정 지표에 들어간다", "nef1" in r["checks"], f"{list(r['checks'])}")
check("천장에 닿은 경우 nef1 체크 통과", r["checks"]["nef1"] is True, f"nef1={r['nef1']}")

# 활성 전부 최하위: nEF1% 0.0 → 미달
r0 = evaluate(rows_with_actives_on_top(2000, 60, 0), [])
check("활성이 전부 바닥이면 nef1 미달", r0["checks"]["nef1"] is False, f"nef1={r0['nef1']}")

print("\n[5] 판정 — 대조군은 문턱의 PASS 요구에서 면제된다")
rc = evaluate(rows_with_actives_on_top(2000, 60, 0), [], arm="control")
check("대조군 판정은 CONTROL_REFERENCE", rc["verdict"] == "CONTROL_REFERENCE", f"{rc['verdict']}")
check("대조군은 FAIL 로 떨어지지 않는다", rc["verdict"] != "FAIL")
check("대조군도 지표는 그대로 낸다", rc["nef1"] is not None, f"nef1={rc['nef1']}")
check("대조군에 문턱을 적용하지 않았음이 기록된다", rc["thresholds_applied"] is False)
check("치료군에는 문턱을 적용했음이 기록된다", r["thresholds_applied"] is True)

print("\n[6] 활성 비율 과다 안전장치 — nEF1% 로 바꿔도 죽지 않았는가")
# 예전에는 'EF1% 천장이 문턱 아래' 라는 간접 조건이 이 역할을 했다.
# nEF1% 는 천장으로 나누므로 그 조건이 영구히 거짓이 된다.
over = int(2000 * (ACTIVE_FRACTION_MAX + 0.05))
ri = evaluate(rows_with_actives_on_top(2000, over, over), [])
check(
    f"활성 비율 {100 * (ACTIVE_FRACTION_MAX + 0.05):.0f}% 는 INVALID",
    ri["verdict"] == "INVALID",
    f"{ri['verdict']} (활성비율 {ri['active_fraction']})",
)
under = int(2000 * (ACTIVE_FRACTION_MAX - 0.15))
ru = evaluate(rows_with_actives_on_top(2000, under, under), [])
check(
    f"활성 비율 {100 * (ACTIVE_FRACTION_MAX - 0.15):.0f}% 는 INVALID 가 아니다",
    ru["verdict"] != "INVALID",
    f"{ru['verdict']}",
)

print("\n[7] 보고서 생성이 현재 문턱으로 죽지 않는가")
# run 43961801413 에서 core/report/build.py 가 val['thresholds']['ef1'] 을
# 손으로 조회해 KeyError 로 죽었다. 도킹 32샤드가 전부 끝난 뒤 보고서에서
# 멈춰 판정이 기록되지 않았다. [2] 의 소스 검사는 scripts/ 만 봤기 때문에
# 이 경로를 놓쳤다 — 여기서는 실제로 렌더링해 본다.
import json as _json  # noqa: E402
import tempfile  # noqa: E402

from core.report.build import target_section  # noqa: E402

_TARGET = {
    "gene": "SELFCHECK",
    "disease": "자체 점검",
    "class": "receptor",
    "rationale": "게이트 자체 점검용 합성 입력",
}
_META = {
    "source": "SELFCHECK",
    "pdb_id": "0XXX",
    "receptor": {"chain": "A", "n_atoms": 1},
    "pocket": {"mode": "cocrystal_ligand", "confidence": "high", "note": "합성"},
}


def _write_case(p: Path, verdict: dict) -> None:
    (p / "target_meta.json").write_text(_json.dumps(_META), encoding="utf-8")
    (p / "validation.json").write_text(_json.dumps(verdict), encoding="utf-8")


_probe = evaluate(rows_with_actives_on_top(2000, 60, 30), [])
with tempfile.TemporaryDirectory() as _td:
    _p = Path(_td)
    _write_case(_p, _probe)
    try:
        _out = target_section(_TARGET, _p)
        check("치료군 판정으로 보고서가 렌더링된다", bool(_out), f"{len(_out)}자")
        for _k in THRESHOLDS:
            check(f"보고서에 {_k} 문턱이 적힌다", str(THRESHOLDS[_k]) in _out)
    except Exception as _e:
        check("치료군 판정으로 보고서가 렌더링된다", False, f"{type(_e).__name__}: {_e}")

_probe_c = evaluate(rows_with_actives_on_top(2000, 60, 0), [], arm="control")
with tempfile.TemporaryDirectory() as _td:
    _p = Path(_td)
    _write_case(_p, _probe_c)
    try:
        _outc = target_section(_TARGET, _p)
        check("대조군 판정으로도 보고서가 렌더링된다", bool(_outc))
        check(
            "보고서가 대조군 면제를 밝힌다",
            "음성 대조군" in _outc,
            "면제를 적지 않으면 낮은 수치가 결함으로 읽힌다",
        )
    except Exception as _e:
        check("대조군 판정으로도 보고서가 렌더링된다", False, f"{type(_e).__name__}: {_e}")

print("\n[8] 실측 대조 — ci-results 에서 읽은 값과 맞는가")
# EPAS1 치료군 41259689466: n=1814, 활성 52, EF1% 29.071 -> nEF1% 0.833
# EPAS1 대조군 41259689466: n=519,  활성 14, EF1%  7.414 -> nEF1% 0.200
for label, n, na, ef1, want in (
    ("EPAS1 치료군(41259689466)", 1814, 52, 29.071, 0.8333),
    ("EPAS1 대조군(41259689466)", 519, 14, 7.414, 0.2000),
):
    got = normalized_ef(ef1, n, na, 0.01)
    check(f"{label} nEF1% = {want}", got is not None and abs(got - want) < 0.001, f"계산값 {got}")

thr = THRESHOLDS["nef1"]
check(
    f"치료군 0.8333 은 문턱 {thr} 통과 / 대조군 0.2000 은 미달",
    0.8333 >= thr > 0.2000,
    f"문턱 {thr}",
)

print()
if FAILURES:
    print(f"[!] 자체 점검 실패 {len(FAILURES)}건")
    for f in FAILURES:
        print(f"      - {f}")
    sys.exit(1)
print("[*] 자체 점검 전부 통과")
