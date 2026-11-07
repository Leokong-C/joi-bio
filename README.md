# ⚠️ `ci-results` — 실행 로그와 중간 산출물. 후보 추천이 아닙니다

이 브랜치는 GitHub Actions 실행이 남긴 **원시 출력**입니다. 코드는
[`main`](https://github.com/Leokong-C/joi-bio/tree/main) 에 있습니다.

## 여기 있는 화합물 목록을 후보로 읽지 마십시오

- 학습도 단백질도 구조도 쓰지 않는 **기억 베이스라인(AUC 0.9637)이
  도킹(0.7717)을 이깁니다.** 이 벤치마크로는 도킹이 작동한다는 것을 증명할
  수 없습니다
- 일부 `candidates_top20.jsonl` 의 화합물은 같은 실행에서 **decoy(음성)로
  라벨된 것들**입니다 — 후보와 음성을 같은 집합에서 뽑던 시기의 산출물입니다
- 구조 선택 편향 재측정은 **두 번 실패**했습니다(2026-10 라이브러리 결함,
  2026-11-01 timeout 으로 16샤드 중 7개 취소). 간판 숫자는 아직 예비치입니다
- **실험 확인이 하나도 없습니다**
- 전부 공개 ChEMBL 화합물이고 신규 구조가 아닙니다

## 실패한 실행도 그대로 둡니다

음성 결과를 지우지 않는 것이 이 프로젝트의 규칙입니다. `runs/` 아래에는
취소되거나 실패한 실행의 로그도 남아 있고, 그 산출물은 불완전합니다.
어떤 실행이 유효한지는 `main` 의 `docs/purpose.md` 를 보십시오.

## 이 저장소의 실제 산출물

"가상 스크리닝 숫자 중 얼마가 인공물인가" 를 재는 도구들입니다 — `main` 의
`scripts/leakage_audit.py`, `scripts/trivial_baselines.py`,
`core/validation/gate.py`.

- 공개 규칙: [`docs/disclosure-policy.md`](https://github.com/Leokong-C/joi-bio/blob/main/docs/disclosure-policy.md)
- 목적 점검 1~6차: [`docs/purpose.md`](https://github.com/Leokong-C/joi-bio/blob/main/docs/purpose.md)
