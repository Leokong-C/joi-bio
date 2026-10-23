# FDA·EMA 10원칙 대조표

출처: [Guiding Principles of Good AI Practice in Drug Development](https://www.fda.gov/about-fda/artificial-intelligence-drug-development/guiding-principles-good-ai-practice-drug-development)
(CDER·CBER + EMA, 2026-01-14). 전문 PDF를 열어 각 원칙의 본문을 대조했다.

작성 2026-10-23. 자평이므로 외부 감사가 아니다.

| # | 원칙 | 현재 | 근거 / 할 일 |
| --- | --- | --- | --- |
| 1 | Human-centric by design | 해당 적음 | 환자 대면 시스템이 아니다. CEO 승인 게이트가 사람이 끊는 지점 |
| 2 | Risk-based approach | **부분** | 검증 게이트가 위험 기반 판정에 가깝다. 다만 context of use별 차등이 없다 |
| 3 | Adherence to standards | 미흡 | GxP 대상 아님. 표준 지표(EF/AUC/BEDROC, Truchon α=20)는 따른다 |
| 4 | Clear context of use | **충족** | `docs/context-of-use.md` 작성(2026-10-23). 이전에는 없었다 |
| 5 | Multidisciplinary expertise | 미흡 | 1인 + AI. 도메인 검토자가 없다. 외부 리뷰 필요 |
| 6 | Data governance and documentation | **강점** | 모든 실행이 ID로 추적. 로그는 ci-results, 결과는 exp/. 24개 실행 전수 기록 |
| 7 | Model design and development | **부분** | 설계 근거를 문헌으로 댄다. 해석가능성은 포즈 수준까지. 일반화는 미검증 |
| 8 | Risk-based performance assessment | **부분** | 교란 감사 도입(탐지율 2/4). 시스템 전체 평가의 첫 걸음이나 미탐지 2건이 남았다 |
| 9 | Life cycle management | **강점** | 버전업 규율(v1.0→v3.0), 음성 결과와 오류 목록 유지. 데이터 드리프트 감시는 없다 |
| 10 | Clear, essential information | **부분** | 한계를 명문화했다(context-of-use). 비전문가용 설명은 아직 |

## 원칙별 메모

**2번 — 위험 기반.** 원문은 "proportionate validation, risk mitigation, and
oversight based on the context of use"를 요구한다. 지금은 모든 실행에 같은 게이트를
적용한다. 탐색용과 보고용을 가르는 차등이 없다.

**6번 — 데이터 거버넌스.** 원문은 데이터 출처·처리 단계·분석 선택이 추적·검증
가능해야 한다고 본다. 이 항목이 가장 강하다. 실행 ID, 커밋 해시, 러너분, 비용이
전부 기록되어 있고 틀린 가설과 실패도 지우지 않았다.

**8번 — 성능 평가.** 원문은 "how humans interact with the AI"를 포함한 시스템 전체
평가를 요구한다. 교란 감사가 여기 해당하는데 **SMILES 교체와 철회 레코드를 못 잡는다.**
스코어와 구조의 대응을 확인하는 장치가 없다는 뜻이다.

**9번 — 생애주기.** 원문은 데이터 드리프트를 잡는 정기 재평가를 요구한다. ChEMBL은
계속 갱신되고 PDB에는 새 구조가 들어온다. 지금은 재평가 일정이 없다. 시간 분할
재검증(Sheridan time-split)이 이 항목의 첫 수단이 될 수 있다.

## 가장 급한 두 가지

1. **8번** — 교란 감사 미탐지 2건. 입력 무결성 검사가 없다.
2. **5번** — 도메인 검토자 부재. 내가 틀린 것 8건 중 상당수는 외부 눈이 있었으면
   더 빨리 잡혔을 것이다.
