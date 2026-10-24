# decoy 풀 출처 — ZINC 확인 결과와 대체 경로

2026-10-25 작성. `docs/purpose.md` 2차 점검이 지목한 "다음 한 가지: 풀 출처를 ZINC로 바꿀 수
있는지 확인"에 대한 답입니다. 모든 수치는 직접 호출해 얻었습니다.

## 1. 문제 (다른 트랙의 측정, 인용)

`docs/purpose.md` 2차 점검 기록:
- decoy 물성 매칭을 적용해 불소 단독 AUC를 **0.9114 → 0.8164**로 내렸으나 기준 0.60에 미달.
- 원인: ChEMBL 무작위 풀 **2,097개 중 활성 서명(불소≥3 + 설폰 + 나이트릴 + 다이아릴에테르)을
  모두 가진 화합물이 0개**. 활성은 35 %가 보유. 맞출 상대가 없으면 매칭이 작동하지 않음.
- 결론: ChEMBL이 이 활성군의 decoy 출처로 너무 작음. DUD-E가 ZINC(약 2.3억)를 쓰는 이유.

## 2. ZINC — 네트워크는 열렸으나 사이트가 자동 접근을 차단합니다

| 호스트 | 허용 목록 | 실제 응답 |
| --- | --- | --- |
| `zinc20.docking.org` | **승인 완료** | HTTP 200이지만 본문이 "Please verify you're human" 인증 화면 (ZINC21 0.20210727.1) |
| `zinc22.docking.org` | **승인 완료** | 프록시 502 Bad Gateway |
| `cartblanche22.docking.org` / `files.docking.org` | — | 연결 시간 초과 |
| `zinc15/zinc.docking.org`, `tldr.docking.org` | 미승인 | 허용 목록 차단 |

`zinc20.docking.org`는 `/substances.json`, `/substances/ZINC…json`, `/tranches/…` 어느
경로를 불러도 동일한 12.3 KB 인증 페이지를 돌려줍니다. 이것은 허용 목록이 아니라 **ZINC 쪽의
봇 차단**입니다. 사이트의 자동 접근 차단은 정책이므로 User-Agent 위장이나 미러 우회는
하지 않았습니다.

**결론: ZINC로의 출처 교체는 현재 경로로는 불가능합니다.**

## 3. PubChem — 같은 역할을 하며 도달 가능합니다

PubChem PUG-REST는 봇 차단이 없고 허용 목록 안에 있습니다. 물성 매칭 decoy 선별에 필요한
두 기능이 모두 동작합니다.

| 기능 | 엔드포인트 | 확인 |
| --- | --- | --- |
| 부분구조(SMARTS) 검색 | `/compound/fastsubstructure/smarts/{smarts}/cids/JSON?MaxRecords=N` | 200 |
| 물성 일괄 조회 | `/compound/cid/{cids}/property/MolecularWeight,XLogP,TPSA,RotatableBondCount,MolecularFormula/CSV` | 200 |

### 결정적 수치

활성 서명 전체를 SMARTS `C(F)(F)F.S(=O)(=O).C#N.c-O-c` 로 질의한 결과:

| 출처 | 서명 전체 보유 화합물 |
| --- | --- |
| ChEMBL 무작위 풀 2,097개 | **0개** |
| PubChem | **3,000개 (MaxRecords 상한에 걸림 — 실제로는 그 이상)** |

점(`.`)으로 이은 SMARTS가 정말 같은 분자 안에서 모두 만족되는지 의심스러워 표본 30개의 SMILES를
받아 RDKit으로 직접 검증했습니다 — 불소 개수 ≥3, `S(=O)(=O)`, `C#N`, `c-O-c` 네 조건 동시 충족:

**27 / 30 통과.** 탈락 3건(CID 18768839, 56681835, 71711765)은 모두 다이아릴에테르 조건만
어긋났습니다 — PubChem과 RDKit의 방향족 결합 해석 차이로 보입니다. 즉 PubChem 질의는
**약 90 % 정밀도의 1차 필터**이며, 최종 판정은 RDKit으로 걸러야 합니다. 어차피 필요한 단계입니다.

예시(CID 10479060): `C1=CC(=CC(=C1)OS(=O)(=O)CCCC(F)(F)F)OC2=CC=CC(=C2C#N)C(F)(F)F` — `C18H13F6NO4S`

## 4. 권장 절차

1. 활성 집합에서 물성 서명을 SMARTS로 추출한다(현재: 불소≥3 + 설폰 + 나이트릴 + 다이아릴에테르).
2. PubChem `fastsubstructure`로 후보 CID를 `MaxRecords` 단위로 긁는다.
3. `property/...CSV`로 MW·XLogP·TPSA·회전결합을 받아 활성과 물성 구간을 맞춘다.
4. **RDKit으로 서명을 재검증한다** — 3 단계의 약 10 %는 가짜 양성이다.
5. 활성과 2D 유사도가 높은 것은 제외한다(decoy 정의상 물성은 같고 구조는 달라야 한다).
6. 그렇게 만든 풀로 단일 기술자 최대 AUC를 다시 재고, 기준 0.60과 비교한다.

## 5. 남는 한계

- PubChem은 구매 가능성(purchasability)을 ZINC만큼 보장하지 않습니다. 실제 구매가 필요한
  단계에서는 별도 확인이 필요합니다.
- `MaxRecords` 상한 때문에 전수 집계가 아니라 "최소 N개"만 말할 수 있습니다. 정확한 모집단
  크기가 필요하면 분할 질의가 필요합니다.
- ZINC-22는 프록시 502라 서버 측 문제일 가능성이 있어, 시일이 지난 뒤 재시도하면 열릴 수도
  있습니다. 다만 ZINC-20의 인증 게이트는 그와 무관하게 남습니다.

---

## 결정 실행 현황 (2026-10-24 승인분)

| 항목 | 상태 |
| --- | --- |
| 1. ToolUniverse 유지 | **완료** — 유지. 다른 트랙도 2차 점검에서 판정을 "벗어남"에서 "우회했으나 쓸모가 있었다"로 정정(`find` 질의가 ZINC를 먼저 가리켰다는 이유) |
| 2. PR #2 병합하지 않고 닫기 | **대기** — GitHub 자격증명 승인 카드가 6시간 만료되어 미수행 |
| 3. PR #1·#3 재배치 경로 선정 | **대기** — 2번과 같은 사유 |
| 4. abstain-dti 기여 | 현행 유지(읽기 전용 + 노션 기록만). 별도 추천을 내지 않았으므로 기존 제약을 유지 |
| 5. GPU Linux 장비 | 조달 결정 대기 |
| 6. 다음 실험 = decoy | **완료** — 이 문서가 그 답. ZINC 불가, PubChem 가능 |
| 7. 저장소 위생(`node_modules`) | **대기** — 2번과 같은 사유 |

### 재검증 결과 (요청하신 항목)

| 주장 | 재확인 |
| --- | --- |
| abstain-dti 라이선스 없음 | `/repos/Pseudo-Lab/abstain-dti/license` → **HTTP 404** 재확인 |
| PR 3건 상태 | #1 open(8파일 +1055) · #2 open(42파일 +7010) · #3 open(40파일 +3537) — 모두 미병합 |
| ToolUniverse 상시 토큰 비용 | 카탈로그 설명 636자 ≈ **159 토큰**. `skill()` 로드 시 818 토큰. 워크플로 1건 중앙값 3,542 토큰 |
| `main` 최신 상태 | 14:40 목적 점검 1차 → 14:50 **2차**(decoy 매칭 수행, ToolUniverse 판정 정정) |
