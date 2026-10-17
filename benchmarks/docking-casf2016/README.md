# CASF-2016 도킹·스코어링 벤치마크 — 1단계 보고

작성: 2026-10-16 · 데이터셋 확정 및 기준선 수립 단계

## 1. 목적과 평가 축

단백질–저분자 복합체 285개로 구성된 CASF-2016 코어셋에서 도킹·스코어링 성능을 측정한다.
CASF 설계의 핵심은 **스코어링 과정을 도킹 과정에서 분리**해 네 가지 능력을 따로 측정하는 데 있다 —
scoring power, ranking power, docking power, screening power.[^casf2016]
본 작업은 이 중 처음 세 가지를 대상으로 한다.

| 평가 축 | 측정 대상 | 결정 구조만으로 가능한가 |
|---|---|---|
| Scoring power | 예측 점수와 실험 logKa의 상관 | 가능 (결정 포즈 사용) |
| Ranking power | 같은 타깃의 리간드 5개 순위 재현 | 가능 |
| Docking power | 재도킹 포즈가 결정 포즈를 2 Å 이내 복원 | **불가 — 도킹 엔진 필요** |

평가 지표의 조작적 정의는 CASF-2013 프로토콜 논문의 절차를 따른다.[^nprot]
벤치마크 계보는 2009년 초판[^casf2009] → CASF-2013[^casf2013] → CASF-2016[^casf2016]이며,
기반 데이터는 PDBbind 데이터베이스다.[^pdbbind]

## 2. 데이터셋 (확정)

| 항목 | 값 |
|---|---|
| 복합체 수 | 285 |
| 타깃 클러스터 | 57개 × 리간드 5개 (모든 클러스터가 정확히 5개) |
| 실험 결합세기 logKa | 2.07 – 11.82 (평균 6.49, 표준편차 2.17) |
| 측정 유형 | Ki 178건 / Kd 107건 |
| 결정 해상도 | 1.10 – 2.50 Å |
| 등록 연도 | 1989 – 2015 |

![데이터셋 개요](./casf2016_dataset_overview.png)

왼쪽 패널은 57개 클러스터의 logKa 분포(파란 눈금 = 클러스터 중앙값), 오른쪽은 리간드 크기 기준선이다.

## 3. 구조 확보와 리간드 식별

285개 mmCIF 전부를 RCSB에서 내려받았다(189.3 MB, 실패 0건).
wwPDB 공식 아카이브(`files.wwpdb.org`)는 이 환경의 허용 목록 밖이지만 RCSB가 동일 아카이브를 제공하므로 대체했다.

PDBbind의 리간드 지정 파일이 없으므로 리간드를 구조에서 직접 식별했다. 절차는 다음과 같다.

1. mmCIF의 `_atom_site`에서 HETATM 그룹을 (체인, 잔기번호, 화합물 ID)로 묶는다.
2. 물·이온·결정화 첨가제를 제외하고, 중원자 6개 미만 그룹을 버린다.
3. `_chem_comp.type`이 펩타이드/당/핵산 결합형인 성분(즉 고분자 사슬의 수식 잔기)을 제외한다.
4. 남은 후보 중 중원자 수가 가장 많은 것을 선택하고, 동수면 단백질 접촉 원자 수로 가린다.
5. RCSB가 결합친화도 주석으로 리간드 화합물 ID를 제공하는 엔트리는 그 값을 우선한다.

**검증:** RCSB 주석이 존재하는 50개 엔트리에서 휴리스틱 단독 결과와 비교해 **49/50 일치**했다.
3단계(`_chem_comp.type` 필터)를 넣기 전에는 48/50이었고, 불일치는 공유결합된 수식 아미노산(LLP, TYS)을
리간드로 잘못 고른 경우였다.

| 리간드 출처 | 엔트리 수 |
|---|---|
| RCSB 결합친화도 주석 | 49 |
| `_chem_comp.type` 필터 휴리스틱 | 229 |
| 최대 HETATM 휴리스틱 (지방산 리간드 등) | 5 |
| 주석-휴리스틱 불일치, 휴리스틱 채택 (3syr) | 1 |
| **미해결** (1a30, 리간드가 펩타이드 사슬) | 1 |

결합부위 탐색 상자는 선택된 리간드의 좌표 범위에 각 축 +10 Å 여유를 준 크기로 정의했다 —
중앙값 17.2 × 17.1 × 17.0 Å, 최대 변 30.1 Å. 285개 전체 표는 [casf2016_core_index.csv](./casf2016_core_index.csv).

육안 확인용 대표 구조: [3o9i (최강, Ki 1.5 pM)](https://www.rcsb.org/structure/3O9I) · [4kzu (중앙값, Ki 313 nM)](https://www.rcsb.org/structure/4KZU) · [3ao4 (최약, Kd 8.5 mM)](https://www.rcsb.org/structure/3AO4).

## 4. 실행 환경 (확정)

이 머신은 네이티브 Windows(16코어 / 28 GiB / GPU 없음)이고 연결된 원격 컴퓨트 대상은 없다.
**재도킹은 이 환경에서 실행 가능하다.** 제약은 엔진이 아니라 Python 바인딩에 있다 —
PyPI의 `vina` 1.2.7은 manylinux·musllinux 휠만 배포하므로 `pip install vina`는
`No matching distribution found for vina`로 실패하고, sdist 빌드는 Boost·SWIG·MSVC를 요구하는데
이 샌드박스에는 컴파일러가 없다(`cl`·`gcc`·`clang`·`cmake`·`swig` 모두 부재). conda-forge·bioconda에도
win-64 빌드가 없다. 그러나 AutoDock Vina[^vina][^vina12] 공식 릴리스는 Windows 실행파일
`vina_1.2.7_win.exe`(1,233,920 B)를 제공하며, 내려받아 실행한 결과 `AutoDock Vina v1.2.7`을
정상 출력했다(exit 0). 따라서 파이프라인을 `vina.Vina()` API 대신 CLI 호출로 구성하면 된다.
보조 도구인 Meeko(PDBQT 준비), sPyRMSD(대칭 보정 RMSD), PoseBusters도 순수 파이썬 패키지라
Windows에서 설치된다.

도킹 포즈를 얻은 뒤에는 RMSD 성공률만으로 충분하지 않다 — 물리적으로 불가능한 포즈가 성공으로
집계되는 문제가 보고되어 있으므로 PoseBusters 류의 타당성 검사를 함께 적용한다.[^posebusters]


## 5. 지금 확보한 수치: 기준선 (널 모델)

결정 구조 포즈만으로 계산 가능한 자명한 특징이 실험 결합세기를 얼마나 설명하는지 측정했다.
스코어링 함수가 이 값을 넘지 못하면 "상호작용을 이해한다"고 말할 수 없는 바닥선이다 (n = 284).

| 특징 | Pearson R | Spearman ρ | 타깃 내 평균 Spearman |
|---|---|---|---|
| 리간드 중원자 수 | **+0.469** (p = 5.8×10⁻¹⁷) | +0.493 | +0.557 (56개 클러스터) |
| 단백질 접촉 원자 수 (4.5 Å) | +0.464 (p = 1.3×10⁻¹⁶) | +0.470 | +0.521 (57개 클러스터) |
| 접촉 수 / 중원자 수 | −0.092 (p = 0.12) | −0.103 | −0.136 |
| 결정 해상도 (대조군) | +0.173 | — | — |

해석: 리간드 크기 하나로 scoring power R ≈ 0.47, ranking power(타깃 내 평균 ρ) ≈ 0.56이 나온다.
접촉 밀도(접촉 수를 크기로 나눈 값)는 유의하지 않아, 설명력이 크기 자체에서 오고 접촉의 질에서 오지 않음을 보여준다.
결정 해상도가 logKa와 약하게 양의 상관을 보이는 것은 강한 결합자가 더 좋은 결정을 주는 경향에서 오는
교란 변수이며, 성능 수치 해석 시 유의해야 한다.

## 6. 2단계 결과 — AutoDock Vina 1.2.7 재도킹 (n = 284)

실행 환경은 이 Windows 머신이다. 공식 릴리스 실행파일 `vina_1.2.7_win.exe`를 CLI로 호출했고,
설정은 네이티브 리간드 중심 상자(좌표 범위 + 각 축 10 Å), `exhaustiveness 8`, `num_modes 9`,
`seed 42`, 작업당 `--cpu 1`을 14개 병렬로 돌렸다. 284개 전부 성공, 실패 0건, 총 44분(복합체당 중앙값 47.9초).
수용체는 ATOM 레코드만 남겨 Open Babel로 pH 7.4 PDBQT로 변환했고, 리간드는 RCSB ModelServer에서
결합차수가 포함된 SDF를 받아 Meeko로 PDBQT를 만들었다. RMSD는 sPyRMSD로 대칭을 보정해 계산했다.

![재도킹 결과](./casf2016_vina_benchmark.png)

### Docking power (포즈 복원)

| 기준 | top-1 포즈 | 9개 모드 중 최선 |
|---|---|---|
| RMSD ≤ 1 Å | 52.1 % | 59.5 % |
| **RMSD ≤ 2 Å** | **66.2 %** | **86.6 %** |
| RMSD ≤ 3 Å | 76.1 % | 94.4 % |

top-1 RMSD 중앙값 0.95 Å, 평균 2.13 Å. top-1과 best-of-9의 간격(66.2 % → 86.6 %)은
**탐색이 아니라 순위 매기기가 병목**임을 뜻한다 — 올바른 포즈는 대개 9개 안에 들어 있지만
점수가 그것을 1위로 올리지 못한다.

### 포즈의 물리적 타당성 (PoseBusters)

RMSD를 뺀 27개 검사를 모두 통과한 비율은 96.1 %다. 실패는 내부 에너지(98.2 % 통과),
사면체 키랄성(99.3 %), 비방향족 고리 평면성(99.3 %), 내부 입체 충돌(99.6 %),
이중결합 평면성(99.6 %)에서 나왔다. 타당성과 정확도를 함께 요구하면
**PB 통과 AND RMSD ≤ 2 Å = 64.4 %**다. PoseBusters가 자체 계산한 RMSD ≤ 2 Å 판정은
우리 sPyRMSD 결과와 284개 전부 일치했다.

### Scoring power (결정 포즈 점수 vs 실험 logKa)

| 예측자 | Pearson R | Spearman ρ |
|---|---|---|
| Vina 점수 (결정 포즈) | **+0.495** (p = 5.9×10⁻¹⁹) | +0.551 |
| Vina 점수 (도킹 최고 모드) | +0.526 | — |
| 리간드 중원자 수 (널 모델) | +0.469 | +0.493 |

잔차 표준편차는 1.89 log 단위다. **Vina의 스코어링은 리간드 크기만 쓰는 널 모델을 거의 넘지 못한다**
(R 0.495 대 0.469). 결합세기 예측의 설명력 대부분이 크기 효과라는 뜻이다.

### Ranking power (타깃 내 리간드 5개 순위)

| 예측자 | 타깃 내 평균 Spearman ρ | ρ ≥ 0.9 비율 |
|---|---|---|
| Vina 점수 | +0.430 (57개 타깃) | 15.8 % |
| 리간드 중원자 수 (널 모델) | +0.557 (56개 타깃) | 14.0 % |

두 예측자를 같이 계산할 수 있는 56개 타깃에서는 Vina +0.432 대 널 모델 +0.557이다.
**같은 단백질에 붙는 리간드들의 순위에서는 Vina가 널 모델보다 못하다.** 크기 차이가 작아지는
동일 타깃 내부에서는 Vina가 기댈 신호가 사라진다는 뜻으로 읽힌다.

### 종합

포즈를 찾는 능력(66.2 %)과 세기를 맞추는 능력(R = 0.495)은 성격이 다른 문제이고,
이 벤치마크에서 Vina는 전자에서는 쓸 만하지만 후자에서는 자명한 기준선과 구별되지 않는다.
단적인 예가 3k5v(Abl 키나아제–imatinib)로, 점수는 −10.98 kcal/mol로 데이터셋 상위권인데
포즈는 결정 구조에서 11.92 Å 벗어나 있다 — 점수가 좋다는 것이 자리를 맞혔다는 뜻이 아니다.

## 7. 다음 단계

1. Vinardo 등 재채점 함수로 같은 포즈를 다시 매겨 스코어링 성능이 개선되는지 확인.
2. exhaustiveness를 32까지 올려 top-1 ↔ best-of-9 간격이 탐색 부족 때문인지 점수 때문인지 분리.
3. 1a30(펩타이드 리간드) 별도 처리.
4. 결합부위 상자 크기 민감도 점검(현재 리간드 범위 + 10 Å).

## 참고문헌

모든 DOI는 Crossref에서 확인했고(10/10 해석 성공) 제목은 PubMed 기록과 정규화 후 일치했다.
검증 표: [reference_check.csv](./reference_check.csv).

[^casf2016]: Su M. et al. *Comparative Assessment of Scoring Functions: The CASF-2016 Update.* Journal of Chemical Information and Modeling (2018), DOI [10.1021/acs.jcim.8b00545](https://doi.org/10.1021/acs.jcim.8b00545)
[^casf2013]: Li Y. et al. *Comparative assessment of scoring functions on an updated benchmark: 1. Compilation of the test set.* Journal of Chemical Information and Modeling (2014), DOI [10.1021/ci500080q](https://doi.org/10.1021/ci500080q)
[^casf2009]: Cheng T. et al. *Comparative assessment of scoring functions on a diverse test set.* Journal of Chemical Information and Modeling (2009), DOI [10.1021/ci9000053](https://doi.org/10.1021/ci9000053)
[^nprot]: Su M. et al. *Assessing protein–ligand interaction scoring functions with the CASF-2013 benchmark.* Nature Protocols (2018), DOI [10.1038/nprot.2017.114](https://doi.org/10.1038/nprot.2017.114)
[^pdbbind]: Wang R. et al. *The PDBbind Database: Collection of Binding Affinities for Protein–Ligand Complexes with Known Three-Dimensional Structures.* Journal of Medicinal Chemistry (2004), DOI [10.1021/jm030580l](https://doi.org/10.1021/jm030580l)
[^vina]: Trott O., Olson A. J. *AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading.* Journal of Computational Chemistry (2009), DOI [10.1002/jcc.21334](https://doi.org/10.1002/jcc.21334)
[^vina12]: Eberhardt J. et al. *AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings.* Journal of Chemical Information and Modeling (2021), DOI [10.1021/acs.jcim.1c00203](https://doi.org/10.1021/acs.jcim.1c00203)
[^vinardo]: Quiroga R., Villarreal M. A. *Vinardo: A Scoring Function Based on AutoDock Vina Improves Scoring, Docking, and Virtual Screening.* PLOS ONE (2016), DOI [10.1371/journal.pone.0155183](https://doi.org/10.1371/journal.pone.0155183)
[^posebusters]: Buttenschoen M. et al. *PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel sequences.* Chemical Science (2024), DOI [10.1039/d3sc04185a](https://doi.org/10.1039/d3sc04185a)

문헌 메타데이터는 PubMed에서 조회했다. 추가 배경: Liu Z. et al. *Forging the Basis for Developing Protein–Ligand Interaction Scoring Functions*, Accounts of Chemical Research (2017), DOI [10.1021/acs.accounts.6b00491](https://doi.org/10.1021/acs.accounts.6b00491).

### 2단계 산출물

- `casf2016_docking_results.csv` — 284개 복합체의 결정 포즈 점수, 도킹 점수, 9개 모드 RMSD, PoseBusters 판정
- `casf2016_posebusters_top1.csv` — top-1 포즈의 PoseBusters 28개 검사 원자료
- `casf2016_pymol_bundle.zip` — PyMOL 렌더링용 구조 3종과 `render_casf.pml`
