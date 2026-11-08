# Joi_bio — 신장 타깃 가상 스크리닝 보고서

생성 2026-11-08 08:43 UTC · 타깃 10개

## 요약

| 타깃 | 질환 | 구조 | 포켓 신뢰도 | 검증 | 후보 |
| --- | --- | --- | --- | --- | --- |
| EPAS1 | ccRCC | RCSB | high | — | — |
| MET | papillary_RCC | — | — | — | — |
| VHL | ccRCC | — | — | — | — |
| SLC5A2 | CKD_DKD | — | — | — | — |
| CA9 | ccRCC | — | — | — | — |
| AXL | RCC_TKI_resistance | — | — | — | — |
| NLRP3 | kidney_fibrosis | — | — | — | — |
| KEAP1 | oxidative_stress_CKD | — | — | — | — |
| APOL1 | FSGS_CKD | — | — | — | — |
| KCNJ1 | hypertension_diuresis | — | — | — | — |

**실행 1/10, 검증 통과 0.** 검증을 통과한 타깃의 후보만 실험 검증 대상으로 제안합니다.

### 방법의 한계

- 리지드 수용체 도킹입니다. 측쇄 유연성과 유도 적합을 반영하지 않습니다.
- Vina 스코어링 함수는 금속 배위(CA9)와 단백질-단백질 계면(KEAP1, VHL)을 과소평가하는 경향이 있습니다. 이들 타깃의 낮은 농축도는 파이프라인 결함이 아니라 방법의 알려진 한계일 수 있습니다.
- 리간드 프로토네이션 상태는 ChEMBL 표기를 그대로 사용했으며 pH 7.4 기준으로 재계산하지 않았습니다.
- ADMET 값은 경험적 rule 기반이며 학습된 예측 모델이 아닙니다.

## EPAS1 (HIF-2a)

- 질환: ccRCC · 분류: transcription_factor
- 선정 근거: VHL 결손 ccRCC의 핵심 구동자. belzutifan 승인 → 최상의 양성 대조군.
- 구조: RCSB 6CZW, 1.6 Å · 체인 A · 원자 880
- 포켓: cocrystal_ligand (신뢰도 high) — FO7 (24 heavy atoms) 기준

**검증: INCOMPLETE** — 평가셋 불완전 — 라이브러리 2160개 중 14개만 도킹됐다(적용률 0.7%, 문턱 90%). 샤드가 통째로 빠졌을 가능성이 크다. 지표를 측정값으로 쓰지 마라.

| 지표 | 값 | 임계값 |
| --- | --- | --- |
| EF1% | n/a | ≥ 5.0 |
| EF5% | n/a | — |
| AUC | n/a | ≥ 0.7 |
| BEDROC | n/a | ≥ 0.3 |
| 화합물 수 | 14 (활성 0) | 활성 ≥ 10 |

라이브러리에 없는 벤치마크: belzutifan, pt2385, pt2399

후보를 보고하지 않습니다. 검증을 통과하지 못한 예측은 근거가 없습니다.

## MET

- 질환: papillary_RCC · 분류: kinase
- 선정 근거: 유두상 신세포암 주요 구동 키나제. ATP 포켓 → 도킹 난이도 낮음.
- **상태: 미실행**

## VHL

- 질환: ccRCC · 분류: e3_ligase
- 선정 근거: PROTAC 리크루터로 결합 데이터 풍부. PPI 계면이라 난이도 상한 측정용.
- **상태: 미실행**

## SLC5A2 (SGLT2)

- 질환: CKD_DKD · 분류: transporter
- 선정 근거: 신보호 효과 입증된 최대 성공 사례. 막단백질 벤치마크.
- **상태: 미실행**

## CA9 (CAIX)

- 질환: ccRCC · 분류: enzyme
- 선정 근거: ccRCC 특이 표면 마커. 아연 활성부위 → 금속 배위 도킹 검증.
- **상태: 미실행**

## AXL

- 질환: RCC_TKI_resistance · 분류: kinase
- 선정 근거: RCC의 TKI 내성·전이 관여. 차세대 병용 타깃.
- **상태: 미실행**

## NLRP3

- 질환: kidney_fibrosis · 분류: inflammasome
- 선정 근거: 신장 염증·섬유화 핵심. NACHT 도메인 알로스테릭 포켓.
- **상태: 미실행**

## KEAP1

- 질환: oxidative_stress_CKD · 분류: ppi
- 선정 근거: Nrf2 경로. PPI 타깃이라 리지드 도킹의 한계 측정용.
- **상태: 미실행**

## APOL1

- 질환: FSGS_CKD · 분류: pore_forming
- 선정 근거: G1/G2 위험 변이 기반 CKD. inaxaplin 등 신규 영역 → 화이트스페이스.
- **상태: 미실행**

## KCNJ1 (ROMK)

- 질환: hypertension_diuresis · 분류: ion_channel
- 선정 근거: 신규 이뇨 타깃. 심혈관 부작용 적음.
- **상태: 미실행**
