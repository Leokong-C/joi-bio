# 선행연구·특허 조사 (2026-10-25)

샌드박스 네트워크가 오프라인이라 논문 본문은 받지 못했습니다. 아래는 검색으로 **표제와 URL이
확인된** 것과, 그중 수치까지 **재확인된** 것을 구분해 적었습니다. 본문 확인이 필요한 항목은
그렇게 표시했습니다.

## 1. 라벨 잡음 — 우리 측정이 선행연구를 재현합니다

| 문헌 | 확인 수준 |
| --- | --- |
| Kramer C., Kalliokoski T., Gedeck P., Vulpetti A. *The Experimental Uncertainty of Heterogeneous Public Ki Data.* J Med Chem, 2012-05. doi:10.1021/jm300131x | **수치 재확인** — 대규모 데이터셋의 최대 가능 R² = **0.81** (→ 최대 R 0.90) |
| Kalliokoski T. et al. *Comparability of Mixed IC50 Data.* PLoS ONE. doi:10.1371/journal.pone.0061007 | 표제·URL 확인, 본문 미확인 |
| *Combining IC50 or Ki Values from Different Sources Is a Source of Significant Noise.* J Chem Inf Model 2024. doi:10.1021/acs.jcim.4c00049 | 표제·URL 확인, 본문 미확인 |

**우리 결과와의 대조.** Kramer의 최대 R 0.90은 우리가 CASF-2016에서 측정한 천장 구간
**0.88–0.965** 안에 들어오고, 특히 복제 기반 추정치 **0.882**와 거의 같습니다. 서로 다른
데이터(공개 Ki 전반 vs CASF-2016 285쌍), 서로 다른 출처쌍(ChEMBL 내부 vs PDBbind×BindingDB),
14년 간격으로 같은 값에 도달했습니다.

**차별점.** Kramer는 공개 Ki 데이터 *일반*의 불확실성을 쟀습니다. 우리는 **스코어링 함수의
순위를 매기는 데 실제로 쓰이는 바로 그 벤치마크 집합**에서 재고, 같은 복합체에서 모든 예측자를
그 천장에 대어 봤습니다. 이 조합은 검색 결과에서 보이지 않았습니다.

**Kramer의 경고를 적용했습니다.** 반복 인용(한 측정이 여러 레코드로 들어오는 것)을 문헌 단위로
접었습니다 — Ki/Kd 819건 → 독립값 525건(고유 문헌 443편). 오차 SD 0.547 → **0.541**, 천장
0.964 → 0.965. 추정치는 이 함정에 강건합니다.

## 2. 사소한 기술자가 도킹 점수를 이긴다 — **이미 선행연구가 있습니다**

| 문헌 | URL |
| --- | --- |
| *Robustly interrogating machine learning-based scoring functions: what are they learning?* Bioinformatics 41(2), btaf040 | academic.oup.com/bioinformatics/article/41/2/btaf040/7985708 |
| *Machine-learning scoring functions trained on complexes dissimilar to the test set already outperform classical counterparts on a blind benchmark.* Brief Bioinform 22(6), bbab225 | academic.oup.com/bib/article/22/6/bbab225/6308682 |
| *In Need of Bias Control: Evaluating Chemical Data for Machine Learning in Structure-Based Virtual Screening.* J Chem Inf Model. doi:10.1021/acs.jcim.8b00712 | pubs.acs.org |
| *Comparative Assessment of Scoring Functions: The CASF-2016 Update* | — |

표제만으로도 방향이 분명합니다 — 특히 btaf040은 ML 스코어링 함수가 "무엇을 학습하는가"를
묻고 기준선 모델과 대조하는 논문입니다. **우리 핵심 주장(결합을 모르는 기술자가 Vina를 이긴다)은
새로운 발견이 아닙니다.** 본문 확인 후 우리 수치를 그들의 기준선과 나란히 놓아야 합니다.

## 3. 물성 매칭 음성군이 성능을 우연 수준으로 내린다 — **역시 선행연구가 있습니다**

| 문헌 | URL |
| --- | --- |
| *Property-Dependent Scoring Bias in Virtual Screening Validation: Data and Code* (Zenodo 18721291) | zenodo.org/records/18721291 |
| *Deep Reinforcement Learning Enables Better Bias Control in Benchmark for Virtual Screening* | biorxiv.org/content/10.1101/2023.11.03.565594 |
| *Benchmarking Methods and Data Sets for Ligand Enrichment Assessment in Virtual Screening* | researchgate (JCIM 2014) |
| *An Improved Metric and Benchmark for Assessing the Performance of Virtual Screening Models* | arxiv.org/pdf/2403.10478 |

Zenodo 18721291은 표제부터 EPAS1 트랙의 발견과 같은 주제입니다 — 가상 스크리닝 검증에서의
**물성 의존 점수 편향**. 본문 확인이 급합니다.

## 4. 특허 — **찾지 못했습니다**

특허 데이터베이스(Google Patents, EPO OPS, PatentsView, SureChEMBL)가 모두 샌드박스 허용
목록 밖이고, 웹 검색은 특허 문서를 한 건도 반환하지 않았습니다(유일한 USPTO 링크는 머신러닝
오차 추정 보정에 관한 것으로 우리 주제와 무관합니다).

가능한 경로 두 가지:
1. 네트워크 복구 후 Europe PMC의 `SRC:PAT` 소스로 EPO 색인 특허를 검색 — 도달 가능한 유일한 경로.
2. 직접 Google Patents에서 "protein-ligand scoring function", "virtual screening enrichment",
   "binding affinity prediction" 검색.

다만 판단을 미리 말씀드리면, 이 프로젝트의 산출물은 **측정 방법과 음성 결과**여서 특허 대상이
되기 어렵습니다. 특허가 걸릴 만한 것은 스코어링 함수나 스크리닝 파이프라인 자체인데 우리는
그걸 만들고 있지 않습니다. 특허 조사의 실익은 낮고, **선행연구 조사가 훨씬 중요합니다.**

## 5. 반드시 고쳐야 할 방법상 문제 하나

우리가 보고해 온 **top-1 ≤2 Å 66.3 %**는 CASF-2016 공식 docking power 수치와 **비교할 수 없습니다.**
공식 평가는 각 복합체에 대해 *미리 제공된 디코이 포즈 집합*에서 고르는 과제이고, 우리는
**처음부터 재도킹**했습니다. 후자가 훨씬 어렵습니다. 리뷰어가 가장 먼저 지적할 지점이므로
모든 문서에서 "re-docking from scratch, not the CASF decoy-pose selection task"를 명시해야 합니다.

## 6. 그래서 무엇이 남는가

| 주장 | 선행연구 | 우리 몫 |
| --- | --- | --- |
| 사소한 기술자가 도킹 점수를 이긴다 | **있음** (btaf040 등) | 재현 + MM/GBSA까지 확장 |
| ML은 스코어링에 강하고 포즈·스크리닝에 약하다 | **있음** (CASF-2016 update) | 재순위 66.0 → 40.1 % 정량화 |
| 물성 매칭 음성군은 성능을 우연 수준으로 내린다 | **있음** (Zenodo 18721291) | 자체 캠페인(EPAS1)에서 재현 |
| 라벨 잡음이 천장을 정한다 | **있음** (Kramer 2012, R² 0.81) | **CASF-2016에서 직접 측정하고 모든 예측자를 그 천장에 대어 봄** |
| 바닥(사소한 기술자)과 천장(라벨 잡음)을 한 틀에 | 보이지 않음 | **이것이 우리 기여** |
| 더 비싼 물리(MM/GBSA)도 널을 넘지 못한다 (ΔR² 0.001) | 보이지 않음 | 우리 기여 |
| 자기보고 신뢰도가 무작위 셔플과 구별되지 않는다 | 보이지 않음 | 우리 기여 |
| 재현성 바닥 — 시드만 바꿔도 284건 중 27건이 판정 역전 | 보이지 않음 | 우리 기여 |

**논지를 옮겨야 합니다.** "도킹 점수가 사소한 기술자에 진다"는 이미 알려진 사실이므로 헤드라인이
될 수 없습니다. 설 자리는 **"당신의 벤치마크에는 쓸 수 있는 여지가 얼마나 남아 있는가"** —
바닥과 천장을 함께 재서 그 사이의 폭을 보고하는 틀입니다. 우리는 두 값을 모두 측정한 유일한
사례이고, 그 둘을 하나의 점검표로 묶으면 재사용 가능한 방법이 됩니다.
