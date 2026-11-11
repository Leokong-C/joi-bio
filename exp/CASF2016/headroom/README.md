# HEADROOM

**벤치마크 점수를 두 경계 사이에서 읽는다.**

바닥 — 단백질과 리간드를 함께 본 적 없는 모델이 도달하는 값.
천장 — 라벨 자체의 측정 오차가 허용하는 값.
그 사이의 폭이 "이 벤치마크에 실제로 남은 여지"다.

CASF-2016(285개 복합체, 57개 타깃)에서 둘을 같은 라벨로 재고, 모든 예측자를 그 사이에 놓았다.

---

## 1. 한 장 요약

| 경계 | 값 |
| --- | --- |
| 천장 (라벨 잡음) | **0.88 – 0.90** |
| 바닥 (리간드+포켓, 타깃 미노출) | **0.534** |
| 바닥 (리간드+포켓, 타깃 노출) | **0.685** |

| 예측자 | R |
| --- | --- |
| RF-Score v2 | 0.741 |
| RF-Score v1 | 0.713 |
| BindingDB 외부 학습 리간드 전용 (타깃 노출) | 0.679 |
| BindingDB 외부 학습 리간드 전용 (타깃 분리) | 0.551 |
| **Vina** | **0.502** |
| 리간드 중원자 수 단독 | 0.482 |
| **Vinardo** | **0.428** |
| 포켓 조성 단독 (타깃 분리) | 0.241 |

**Vina(0.502)와 Vinardo(0.428)는 둘 다 바닥(0.534) 아래다.**

**학습 리간드를 256개에서 44,754개로 175배 늘려 얻는 것은 +0.04~0.05인데, 같은 타깃을
학습에서 보는 것이 주는 것은 +0.13~0.14다.** 바닥을 정하는 것은 데이터의 양이 아니라
분할이 모델에게 무엇을 외우도록 허용하느냐다.

---

## 2. 실험 로그 (단계별 · 전 수치)

### 2.1 데이터셋

CASF-2016 코어셋. 285개 복합체 · 57개 타깃 × 리간드 5개.
logKa 2.07 – 11.82 (평균 6.49, SD 2.17). 친화도 종류는 **Ki 178건 · Kd 107건**, IC50 없음.

수작업 정정 2건이 들어 있다 — `3syr`의 리간드를 LLP에서 **GPK**로 바로잡았고(top-1 RMSD가
4.41 Å → 0.27 Å로 바뀜, 즉 거짓 실패였다), `1a30`은 삼펩타이드(Glu-Asp-Leu)로 따로 처리했다.
이 둘 때문에 n이 284에서 **285**로 늘었다.

mmCIF의 HETATM 원자 수는 altloc을 이중 계산한다(`3gv9`·`2xnb`·`2yfe`·`2al5`·`3prs` 5건 확인).
중원자 수는 **SDF 기준**으로 센다 — 이 정정으로 리간드 크기 상관이 +0.469에서 **+0.482**로 바뀌었다.

### 2.2 재도킹 (AutoDock Vina 1.2.7)

수용체는 PDB ATOM → `obabel -xr -p 7.4` → PDBQT. 리간드는 RCSB ModelServer SDF →
RDKit AddHs → Meeko → PDBQT. 상자는 리간드 좌표 범위 + 각 축 10 Å. 복합체당 9개 포즈.
RMSD는 sPyRMSD `symmrmsd`(대칭 고려), 물리 타당성은 PoseBusters `redock` 설정.

| 항목 | 값 |
| --- | --- |
| top-1 ≤ 2 Å | **66.3 %** |
| 반환 9개 중 최선 ≤ 2 Å | **86.7 %** |
| top-1 RMSD 중앙값 | 0.935 Å |
| PoseBusters 물리 타당 | 95.8 % |
| 타당하면서 ≤ 2 Å | 64.4 % (n = 284) |

**이 66.3 %는 처음부터 재도킹한 값이며 CASF 공식 docking power와 비교할 수 없다.**
CASF는 "scoring 과정을 docking 과정과 분리"하도록 설계돼 미리 제공된 디코이 포즈 집합에서
고르는 과제다. 훨씬 쉬운 문제다. 이 단서를 빼면 숫자가 오해된다.

#### 탐색량과 시드

| 실행 | top-1 ≤2 Å | 최선 ≤2 Å | 중앙값 | CPU-시간 |
| --- | --- | --- | --- | --- |
| exhaustiveness 8, seed 42 | 66.20 % | 85.92 % | 0.957 Å | — |
| exhaustiveness 8, seed 43 | 65.49 % | 81.69 % | 0.972 Å | — |
| exhaustiveness 8, seed 44 | 66.20 % | 83.45 % | 0.907 Å | — |
| exhaustiveness 32, seed 42 | 66.55 % | 86.27 % | 0.910 Å | 34.81 |

탐색량을 4배로 올려도 +0.35 pt에 그친다. 시드 간 SD는 0.41 pt이고 **284건 중 27건이
시드만 바꿔도 2 Å 판정을 뒤집는다.**

#### 재현성 — 두 가지를 구분해야 한다

| 조건 | 결과 |
| --- | --- |
| **같은 시드로 재실행** | top-1 66.2 % → 66.2 %, **판정 역전 0건**, 284건 중 **250건이 소수점 셋째 자리까지 동일**, 최대 \|Δ\| 0.663 Å |
| **시드를 바꿈** | 성공률 SD 0.41 pt, **27건 판정 역전** |

같은 시드 재실행에서 34건이 셋째 자리에서 달라지는 원인은 Vina가 아니다 —
**Open Babel의 극성수소 배치가 프로세스 간 재현되지 않는다.** 보관된 수용체에 Vina를
`score_only`로 돌리면 원래 값이 정확히 재현된다.

#### 상자 크기 민감도 (층화 60개: 성공 20 · 경계 20 · 실패 20)

| 여유 | 평균 변 | top-1 ≤2 Å | 최선 ≤2 Å | 중앙값 |
| --- | --- | --- | --- | --- |
| +5 Å | 12.4 Å | 51.7 % | 85.0 % | 1.89 Å |
| **+10 Å (채택)** | 17.4 Å | 51.7 % | 80.0 % | 1.82 Å |
| +15 Å | 22.4 Å | 46.7 % | 66.7 % | 2.17 Å |

상자를 키우면 정확도가 떨어진다. +10 → +15는 −5.0 pt로 시드 잡음(±0.41 pt)의 12배다.

### 2.3 바닥 — 무지(無知) 기준선

Durant et al.(Bioinformatics 41(2) btaf040)의 설계를 재현했다.
**LigandBias** = RDKit 2D 기술자 205개. **ProteinBias** = 16 Å 포켓의 아미노산 20종
잔기 수 + 총계 21개(포켓 중앙값 152잔기). **BothBias** = 226개.
RandomForest 500그루, `random_state=0`.

| 모델 | 타깃 미노출 (GroupKFold-10) | 타깃 노출 (KFold-10 shuffle) | 차이 |
| --- | --- | --- | --- |
| LigandBias (205) | 0.502 | 0.641 | +0.139 |
| **ProteinBias (21)** | **0.241** | **0.509** | **+0.268** |
| BothBias (226) | 0.534 | 0.685 | +0.151 |

기술자를 14개에서 205개로 늘려도 바닥은 올라가지 않았다(0.547 → 0.534).
**바닥을 정하는 것은 특징 개수가 아니라 분할이다.**

ProteinBias가 두 배 이상 뛰는 것이 가장 선명한 증거다 — 포켓 조성은 결합을 설명하는 게
아니라 "이 타깃이 어느 타깃인지"를 알려주는 조회 키로 쓰인다.

#### 외부 학습 (BindingDB 44,754 리간드 · 58 타깃 · CASF 분자 210개 제외)

| | 학습 256개 (CASF 내부) | 학습 44,754개 (BindingDB) | 학습량 175배의 이득 |
| --- | --- | --- | --- |
| 타깃 미노출 | 0.502 | 0.551 (LOTO, 57 클러스터, 2,729초) | **+0.049** |
| 타깃 노출 | 0.641 | 0.679 | **+0.038** |
| **타깃을 보는 이득** | **+0.139** | **+0.128** | |

### 2.4 천장 — 라벨 잡음

285개를 RCSB GraphQL로 UniProt에 매핑(**285/285 성공**, 78개 접근번호).
BindingDB UniProt 엔드포인트로 **258,749건** 수집(78개 파일, 75 MB, 실패 0건).
분자 157,481개 중 분자식 사전 필터로 후보 1,092개(43초)를 추린 뒤 InChIKey 계산 →
**210개 분자 일치**(정확 154 · 골격 54 · 모호 2) → **1,864건**의 친화도 연결.

추정량은 Kalliokoski et al. 2013과 동일하다 — 쌍 차이의 SD를 √2로 나눠 측정당 오차를 구한다.

| 항목 | 값 |
| --- | --- |
| 독립 Ki/Kd를 가진 복합체 | 133 / 285 (46.7 %) |
| 두 출처 상관 R | 0.928 |
| 계통 편향 | −0.012 로그 단위 (사실상 없음) |
| 차이의 SD | 0.773 |
| \|Δ\| > 1 로그 단위 | 12.0 % |
| \|Δ\| > 2 로그 단위 | 4.5 % |
| **천장 보고값** | **0.88 – 0.90** |

| 추정 방식 | n | σ_오차 | 천장 | 보고 여부 |
| --- | --- | --- | --- | --- |
| 교차출처, 필터 없음 | 133 | 0.547 | 0.964 | 낙관적 상한으로만 |
| 주석 오류(\|Δ\|≥2.5) 제거 | 133 | 0.351 | 0.985 | 참고 |
| **복제 기반 (BindingDB 내부, n≥3)** | **45** | **0.972** | **0.882** | **헤드라인** |
| 일치 레코드 제거 | 66 | 1.062 | 0.842 | **편향 — 보고하지 않음** |

**왜 0.964를 내렸나.** 일치한 BindingDB Ki/Kd 레코드의 **20.1 %가 PDBbind 라벨과
0.02 로그 이내**다(133개 중 97개 복합체가 그런 레코드를 하나 이상 가짐). 독립 측정이 아니라
같은 1차 측정이 두 DB에 들어간 것이라 일치도가 부풀려졌다. Kalliokoski는 BindingDB가
ChEMBL에서 값을 복사해 온다고 경고한다.

**왜 그것을 걷어내지 않았나.** 라벨과 일치하는 레코드를 골라 버리면 남는 것은 불일치뿐이라
산포가 기계적으로 부풀려진다(0.842가 그 결과다). 편향된 연산이다. 대신 라벨을 전혀 건드리지
않는 복제 기반 추정치를 헤드라인으로 삼았다.

**독립 검증.** Kramer et al. 2012는 대규모 공개 Ki 데이터의 최대 가능 R²를 0.81(= R 0.90)로
추정했다. 우리 0.882와 거의 일치한다 — 다른 데이터, 다른 출처쌍, 14년 간격.

**중복 인용 점검.** Kramer는 ChEMBL 측정 쌍의 90 % 이상이 반복 인용이라고 경고한다.
문헌 단위로 접어 보니 Ki/Kd 819건 → 독립값 525건(고유 문헌 443편), 오차 SD 0.547 → **0.541**.
추정치는 이 함정에 강건하다.

### 2.5 MM/GBSA (285개 중 무작위 47개)

암시적 용매 OBC2, 에너지 최소화. OpenMM + OpenFF + PDBFixer + NAGL.
무작위로 뽑은 것은 중간 결과가 편향되지 않게 하기 위해서다.

| 예측자 | R (같은 47개) |
| --- | --- |
| RF-Score v2 | 0.756 |
| 널 모델 (forest) | 0.549 |
| 리간드 중원자 수 | 0.528 |
| **MM/GBSA** | **0.385** |
| Vina | 0.368 |
| Vinardo | 0.329 |

널을 제거한 뒤 남는 고유 기여:

| 예측자 | ΔR² | p |
| --- | --- | --- |
| **MM/GBSA** | **0.001** | **0.80** |
| Vina | 0.003 | 0.69 |
| Vinardo | 0.001 | 0.80 |
| RF-Score v2 | 0.271 | 3e-6 |

비용은 복합체당 MM/GBSA 251.6초 : Vina 도킹 94.4초 : RF-Score 1.26초 : 널 모델 0.004초.
**2.7배를 더 써도 사소한 기술자를 넘어서는 기여가 통계적으로 0이다.**

성공 41/47(87 %), 실패 7건. OpenMM이 커널을 STATUS_ACCESS_VIOLATION으로 죽이므로
복합체마다 별도 프로세스로 돌렸다.

### 2.6 포즈 재순위 (n = 282)

Vina가 낸 9개 포즈를 RF-Score v2로 다시 줄 세웠다.

| 컷오프 | Vina 1위 | RF-Score 1위 | 반환 9개 중 최선 (상한) |
| --- | --- | --- | --- |
| ≤ 1 Å | 51.8 % | 29.8 % | 59.6 % |
| **≤ 2 Å** | **66.0 %** | **40.1 %** | 86.5 % |
| ≤ 3 Å | 75.9 % | 53.5 % | 94.3 % |

중앙값 RMSD 0.97 → 2.68 Å. RF-Score가 Vina의 1위를 그대로 유지한 경우는 40.4 %.
복합체별로 **15건 구제, 88건 파괴(순 −73)**.

**친화도 예측 1위 함수가 포즈 선택에서는 최악이다.** 복합체 사이에서 세기를 가르는 능력과
복합체 안에서 포즈를 가르는 능력은 서로 다른 문제이고 전이되지 않는다.

### 2.7 언어모델 축 (855건 구조화 예측)

claude-sonnet-5, thinking 비활성, 오류 0건, 2.3분, 652k in / 137k out 토큰.

| 조건 | R | ρ | 잔차 SD | MAE | 타깃 내 ρ | 기권율 | "안다" 보고 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A: 리간드만 | 0.384 | 0.506 | 2.01 | 1.96 | 0.366 | 96.8 % | 4.9 % |
| B: +단백질명 | 0.655 | 0.663 | 1.65 | 1.38 | 0.521 | 20.4 % | 26.0 % |
| C: +PDB ID | 0.689 | 0.692 | 1.58 | 1.30 | 0.554 | 8.8 % | 44.2 % |

**기권 지표** (abstain-dti 프로토콜: binds ⇔ logKa ≥ 6.0, 167 binds / 118 not)

| 조건 | 정확도 | AURC | 신뢰도 무작위 셔플 | 80 % 커버리지 선택 정확도 |
| --- | --- | --- | --- | --- |
| A | 0.495 | 0.459 | 0.508 | 0.499 |
| B | 0.744 | 0.249 | 0.256 | 0.758 |
| C | 0.761 | **0.236** | **0.244** | 0.761 |

셔플 대조와의 차이가 순열 SD(0.012) 수준이다. **자기보고 신뢰도가 자기 오답을 무작위보다
낫게 정렬하지 못한다.** 80 % 커버리지 선택 정확도가 전체 정확도와 같아 기권으로 얻는 것이 없다.

**오염 검사는 음성.** "안다"고 주장한 복합체가 오히려 더 못 맞혔다(R 0.641, n=126 vs
0.730, n=159). PDB ID는 56/285의 인식을 뒤집었지만 MAE 개선은 −0.11 대 −0.08에 그쳤다.

**난이도 축의 함정.** 이진 정확도가 리간드 크기에 따라 0.709 → 0.758 → 0.828로 오르지만,
라벨 자체가 크기와 상관(R +0.482)되어 있어 난이도가 아니라 **유병률**을 재고 있었다.
회귀 MAE는 단조롭지 않다(1.28 → 1.16 → 1.46).

---

## 3. 스스로 내린 정정 기록

| 무엇 | 전 | 후 | 이유 |
| --- | --- | --- | --- |
| 복합체 수 | 284 | **285** | `1a30` 추가, `3syr` 리간드 정정 |
| `3syr` top-1 RMSD | 4.41 Å | **0.27 Å** | 리간드가 LLP가 아니라 GPK였다 — 거짓 실패 |
| 리간드 크기 상관 | +0.469 | **+0.482** | mmCIF HETATM이 altloc을 이중 계산 |
| Vina 크기 통제 편상관 | +0.355 / +0.339 | **+0.343** | n=285 재계산 |
| PoseBusters 타당성 | 96.1 % | **95.8 %** | n 정정 |
| exhaustiveness 32 | 추정치 | **66.5 %** | 실측으로 교체 |
| 기권 지표 AURC (조건 C) | 0.233 | **0.236** | 아티팩트 재확인 |
| **천장** | **0.88 – 0.96** | **0.88 – 0.90** | 공유 측정치 오염 발견 |
| **학습량 배수** | **157배** | **175배** | 44,754/256 = 174.8 — 산술 오류 |
| 바닥 해석 | "널 모델 0.547이 바닥" | **분할마다 다름** | Durant et al. 본문 확인 후 |

**2026-11-07 전면 재검증:** 보고한 33개 수치를 저장된 아티팩트에서 다시 계산해 대조했다.
**불일치 0건.** 결과는 `results/audit.csv`.

---

## 4. 환경과 버전

수치를 좌우하는 조합 (`VERSIONS.json`):

| 패키지 | 버전 |
| --- | --- |
| Python | 3.11.17 |
| RDKit | 2026.09.1 |
| NumPy | 2.4.6 |
| SciPy | 1.17.1 |
| pandas | 3.0.6 |
| scikit-learn | 1.9.1 |
| matplotlib | 3.11.2 |

`environment.yml`은 버전을 고정하지 않는다. 위 조합과 다르면 RandomForest 동점 처리나
RDKit 기술자 정의가 달라져 **소수점 셋째 자리가 움직일 수 있다.**

재도킹까지 하려면 추가로 AutoDock Vina 1.2.7 · Meeko 0.8.0 · Open Babel 3.2.1 ·
sPyRMSD 0.9.0 · PoseBusters 0.6.5 · ODDT 0.4 · gemmi 0.7.5가 필요하다.

### 외부 DB 상태 (2026-11-07 확인)

| 출처 | 버전 | 비고 |
| --- | --- | --- |
| CASF-2016 | 고정 (2018 발표) | 갱신본 없음 |
| PDBbind | v2016 | 최신판은 등록 필요 · 허용 목록 밖 |
| ChEMBL | **ChEMBL_37** (2026-05-01) | 활성 24,527,044 · 화합물 2,921,148 · 타깃 18,552 |
| BindingDB | REST에 버전 미노출 | **2026-10-25 수집** · 월 단위 갱신 |
| RCSB PDB | 상시 갱신 | 친화도 주석 출처가 BindingDB라 독립 출처가 아님 |
| TDC | 데이터 호스트 접근 가능 (2026-11-07 승인) | `dataverse.harvard.edu` HEAD 202 · `tdcommons.ai`는 차단이나 불필요 |

---

## 5. 저장소 구성

```
exp/CASF2016/headroom/
├── README.md               이 문서
├── BACKEND.md              GitHub 을 백그라운드로 쓰는 설계
├── SPEC.md                 틀의 정의
├── MIGRATION.md            다른 PC 이관 체크리스트
├── PARTS_DRAFT.md          게시용 파트 1–8 초안
├── PRIOR_ART.md            선행연구 조사
├── VERSIONS.json           수치를 좌우하는 패키지 버전
├── environment.yml
├── featurize.py            RDKit 2D 기술자 (프로세스 풀 실패 시 직렬 폴백)
├── 01_pull.py              RCSB PDB→UniProt · BindingDB 수집
├── 02_ceiling.py           리간드 대조 · 라벨 잡음 천장
├── 03_floor.py             편향 모델 3종 × 분할 2종
├── 04_external_floor.py    BindingDB 외부 학습 바닥
├── data/
│   ├── casf_inputs.tar.gz      0.2 MB  SMILES · UniProt · 리간드 키 · SDF
│   ├── pockets.tar.gz          7.1 MB  16 Å 포켓 PDB 285 개
│   └── bindingdb_raw.tar.gz    3.9 MB  BindingDB 원자료 (2026-10-25)
├── results/                커밋된 결과 CSV/JSON
└── figures/                파이프라인 도식 · 바닥–천장 띠
```

`casf_index.csv`(CoreSet.dat 유래)는 **담지 않았다** — PDBbind 재배포 가능 여부를
확인하기 전까지 손으로 옮긴다.

`headroom_bindingdb_features.csv`(78.5 MB)도 담지 않았다 — GitHub 경고선 초과,
`featurize.py`로 13분이면 재생성된다.

---

## 6. 실행

```bash
conda env create -f environment.yml && conda activate headroom
python 01_pull.py   --index casf_index.csv --work work          # 2분
python featurize.py --in casf_smiles.csv --out casf_features.csv # 7초
python 02_ceiling.py --index casf_index.csv --smiles casf_smiles.csv --work work   # 2분
python 03_floor.py   --index casf_index.csv --smiles casf_smiles.csv \
                     --ligand-features casf_features.csv --pockets pockets/ --work work  # 1분
python featurize.py --in work/bindingdb_train.csv --out work/bindingdb_features.csv  # 13분
python 04_external_floor.py --index casf_index.csv --smiles casf_smiles.csv \
       --casf-uniprot work/casf_uniprot.csv --matched work/matched_records.csv \
       --train-features work/bindingdb_features.csv --test-features casf_features.csv \
       --work work                                              # 중첩 2분 + 분리 45분
```

네트워크 없이 분석만 하려면 `data/bindingdb_raw.tar.gz`를 풀어 `work/bindingdb/`에 두고
01단계를 건너뛴다.

### 재현 성공 판정

| 값 | 기대 |
| --- | --- |
| 복합체 / 타깃 | 285 / 57 |
| 천장 (복제 기반) | 0.88 ± 0.01 |
| 바닥 리간드+포켓, 타깃 분리 | 0.534 ± 0.01 |
| 바닥 리간드+포켓, 무작위 분할 | 0.685 ± 0.01 |
| 외부 학습, 타깃 중첩 / 분리 | 0.679 / 0.551 ± 0.01 |
| Vina / Vinardo / RF-Score v2 | 0.502 / 0.428 / 0.741 |

2026-11-07에 저장된 특징 파일만으로 바닥 6개 수치를 재계산해 **소수점 셋째 자리까지 일치**함을
확인했다.

---

## 7. 알려진 함정

1. **PDBbind General은 스크립트로 못 받는다.** 등록이 필요하고 미러는 차단. RCSB도 대체재가
   아니다 — 무작위 24개 중 친화도 주석 9건, Ki/Kd는 2건뿐.
2. **잠긴 컨테이너는 프로세스 풀을 거부한다**(Windows error 5). `featurize.py`가 감지해
   직렬로 떨어지고, RDKit 내부 스레딩 덕에 39,000 분자가 약 13분.
3. **BindingDB의 PDB 키 엔드포인트는 HTTP 500**을 돌려준다. UniProt 키를 쓰고 RCSB로 매핑한다.
4. **응답 최상위 키에 오타가 있다** — `getLindsByUniprotsResponse`. 이름으로 찾지 말고
   단일 값을 읽는다.
5. **`Ipc` 기술자가 float32 범위를 넘는다.** 학습 전에 비유한·거대값 컬럼을 걸러야 한다.
6. **한글 그림은 `Malgun Gothic`을 명시**해야 한다(Windows 전용). 위첨자 글리프가 없어
   지수는 mathtext로 쓴다.
7. **ODDT는 Windows에서 세 가지 우회가 필요하다** — `pip install oddt --no-deps`,
   import 전 `np.in1d = np.isin`, RDKit 백엔드. RF-Score v3는 NumPy 2에서 사용 불가.
8. **수용체를 리간드 16 Å 이내로 자르면** ODDT 채점이 약 40배 빨라지고 예측은 비트 단위로
   동일하다(최대 \|차이\| 1.8e-15).
9. **ZINC은 허용 목록이 열려도 사이트가 봇 차단**을 한다("Please verify you're human").
   우회하지 않는다.

---

## 8. 미해결

- 포즈 파싱 실패 2건(`3utu`, `4k18`), MM/GBSA 실패 7건.
- PoseBusters는 n = 284 기준(정정된 복합체 2건 미재검).
- **스크리닝 능력 축**이 통째로 빠져 있다 — CASF 네 축 중 하나.
- **데이터셋 일반화**가 없다. CASF-2016 단일 셋.
- **공동접힘 비교자**(Boltz-2 등)와 **수용체 유연성**(bioemu2)은 GPU 리눅스 장비가 필요하다.
- **TDC 교차검증** — 접근은 열렸으나 아직 돌리지 않았다. BindingDB 단일 출처의
  공유 측정치 오염을 푸는 가장 직접적인 길이다.

---

## 9. 선행연구

이 프로젝트의 발견 중 셋은 **우리가 처음이 아니다.** 게시·투고 시 반드시 함께 표기한다.

- **Durant G. et al.** *Robustly interrogating machine learning-based scoring functions:
  what are they learning?* Bioinformatics 2025, 41(2) btaf040 — 편향만 쓰는 기준선이
  CASF-2016에서 최신 수준에 도달하고, 널리 쓰이는 MLBSF 다섯이 그에 동등하거나 추월당한다.
- **Boyles F., Deane C. M., Morris G. M.** *Learning from the ligand.* Bioinformatics
  2020, 36, 758 — 리간드 전용 RF가 고전 스코어링 함수와 비슷한 수준이며, 리간드의
  평균 친화도를 예측하는 것으로 보인다.
- **Su M. et al.** *Comparative Assessment of Scoring Functions: The CASF-2016 Update.*
  J Chem Inf Model 2018. doi:10.1021/acs.jcim.8b00545 — docking power의 정의, 그리고
  스코어링에 강한 ML 함수가 도킹·스크리닝에서 뒤진다는 서술.
- **Kramer C., Kalliokoski T., Gedeck P., Vulpetti A.** *The Experimental Uncertainty of
  Heterogeneous Public Ki Data.* J Med Chem 2012. doi:10.1021/jm300131x — 대규모 공개
  데이터의 최대 가능 R² = 0.81.
- **Kalliokoski T. et al.** *Comparability of Mixed IC50 Data.* PLoS ONE 2013.
  doi:10.1371/journal.pone.0061007 — 우리가 쓴 σ 추정량과 필터 설계.
- Trott O., Olson A. J. J Comput Chem 2010. doi:10.1002/jcc.21334 (Vina)
- Quiroga R., Villarreal M. A. PLOS ONE 2016. doi:10.1371/journal.pone.0155183 (Vinardo)
- Ballester P. J., Mitchell J. B. O. Bioinformatics 2010. doi:10.1093/bioinformatics/btq112 (RF-Score)
- Eastman P. et al. PLOS Comput Biol 2017. doi:10.1371/journal.pcbi.1005659 (OpenMM 7)
- Buttenschoen M. et al. Chem Sci 2024. doi:10.1039/D3SC04185A (PoseBusters)

**남는 기여**는 바닥과 천장을 **같은 데이터셋·같은 라벨**에서 함께 재고 모든 예측자를 그
사이에 놓은 것, 그리고 MM/GBSA가 널을 넘지 못한다는 것(ΔR² 0.001), 자기보고 신뢰도가
무작위 셔플과 구별되지 않는다는 것, 재현성 바닥선이다.
