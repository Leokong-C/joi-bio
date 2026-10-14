# Joi_bio

AI 효소 설계와 정제 공정 설계를 함께 제공하는 회사를 목표로, 공통 코어 위에
세 모듈을 쌓는다. 지금 구현된 것은 코어와 저분자 결합 모듈이다.

```
core/                    세 모듈이 공유하는 자산
  structure/             구조 조회, 체인 선택, 포켓 정의
  validation/            ★ 검증 게이트 (EF1% / AUC / BEDROC / 벤치마크 회수)
  report/                보고서 생성
modules/
  small_molecule/        ① 저분자 결합 — 구현 완료
  enzyme_eng/            ② 효소 공학·공정 — 스펙만
  protein_design/        ③ 항체·바인더 — 스펙만
```

`core/`가 핵심 자산이다. 특히 검증 게이트는 저분자든 효소든 항체든 "우리 예측이
실제로 맞는가"를 묻는 같은 질문이라 모듈이 바뀌어도 그대로 재사용된다.

## 검증 게이트가 이 저장소의 전부다

알려진 활성 화합물과 승인약을 라이브러리에 섞어 넣고, 도킹이 이들을 상위로
올리는지 측정한다. **통과하지 못한 타깃의 화합물은 보고서에 히트로 실리지
않는다.** 검증 없이 나온 후보 리스트는 파는 순간 신뢰를 잃는다.

| 판정 | 조건 |
| --- | --- |
| INSUFFICIENT | 활성 10개 미만 — 지표가 나와도 통계적으로 무의미 |
| PASS | EF1% ≥ 5.0, AUC ≥ 0.70, BEDROC ≥ 0.30 중 2개 이상 |
| FAIL | 그 외 |

여기에 더해 **벤치마크 회수**를 따로 본다. EPAS1이라면 belzutifan이 상위 1%에
들어오는지다. 지표가 좋아도 승인약을 못 찾으면 그 파이프라인은 믿을 수 없다.

## 실행

GitHub → Actions → **Joi_bio pipeline** → Run workflow → gene에 `EPAS1` 입력.
`all`을 넣으면 10개 타깃 전체를 돈다.

로컬에서 돌리려면 conda 환경이 필요하다. `pip install vina`는 Boost 1.83+의
deprecated 헤더 때문에 빌드가 깨지므로 conda-forge 바이너리를 쓴다.

```bash
micromamba env create -f environment.yml && micromamba activate joi-bio

python -m core.structure.fetch        --gene EPAS1
python -m modules.small_molecule.library --gene EPAS1 --n-decoys 2000
python -m modules.small_molecule.dock    --gene EPAS1
python -m core.validation.gate        --gene EPAS1
python -m modules.small_molecule.admet   --gene EPAS1
python -m core.report.build
```

## 타깃 10개

**신장암(RCC)** EPAS1(HIF-2α), MET, VHL, CA9, AXL
**기전·CKD** SLC5A2(SGLT2), NLRP3, KEAP1, APOL1, KCNJ1(ROMK)

선정 기준은 ① PDB 실험구조 존재 ② ChEMBL 활성 데이터 존재 ③ 신장 질환 연관성.
근거와 리스크는 `config/targets.yaml`에 적어 두었다.

## 미리 예상되는 실패

10개 전부 통과를 기대하면 안 된다. 3~5개가 통과하면 1차로는 성공이다.

- **APOL1** 실험구조가 부족해 AlphaFold 모델 + blind docking이 된다. FAIL이
  정상적인 결과다.
- **KEAP1, VHL** PPI 타깃이라 리지드 도킹으로는 농축도가 낮게 나오기 쉽다.
- **CA9** Zn²⁺ 배위를 Vina 스코어링이 과소평가한다.

어느 타깃이 왜 실패했는지가 다음 단계 설계의 자산이다.

## 방법의 한계

리지드 수용체 도킹이며 측쇄 유연성과 유도 적합을 반영하지 않는다. 리간드
프로토네이션은 ChEMBL 표기 그대로이고 pH 7.4 기준으로 재계산하지 않는다.
ADMET은 경험적 rule이며 학습된 예측 모델이 아니다. 출력 화합물은 전부
`candidate_for_testing`이고, 실험 검증 전에는 히트가 아니다.
