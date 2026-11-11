# 다른 PC로 이관하기 — HEADROOM 체크리스트

이 문서 하나로 새 PC에서 처음부터 돌릴 수 있게 썼습니다. 저장소 경로는
`exp/FLOOR/headroom/` (PR #6).

## 1. 옮겨야 하는 것 / 다시 받으면 되는 것

| 구분 | 항목 | 비고 |
| --- | --- | --- |
| **반드시 옮김** | `casf_index.csv` (code·logKa·target) | CASF-2016 `CoreSet.dat` 유래. 재배포 불가라 직접 옮겨야 함 |
| **반드시 옮김** | `casf_smiles.csv` (code·smiles) | 리간드 식별 결과. 3syr는 GPK, 1a30은 펩타이드로 수작업 정정된 값 포함 |
| **반드시 옮김** | `pockets/<code>.pdb` 285개 | 리간드 16 Å 이내 잔기. 재생성 가능하나 수작업 정정이 반영돼 있음 |
| 다시 받음 | BindingDB JSON 78개 (75 MB) | `01_pull.py`가 2분이면 재수집. **단 월 단위 갱신이라 수치가 달라질 수 있음** |
| 다시 받음 | RDKit 기술자 CSV | `featurize.py`로 재계산 (39,000 분자 ≈ 13분) |
| 다시 받음 | PDB→UniProt 매핑 | `01_pull.py`가 RCSB에서 재생성 |
| 백업 보관 | `backup_labelnoise_20261025.tar.gz` (4.0 MB) | 원자료 75 MB의 압축본. 수치 재현을 보장하려면 이걸 쓰는 편이 안전 |

**중요** — BindingDB를 새로 받으면 천장 수치가 소수점 셋째 자리에서 달라질 수 있습니다.
논문·게시물 수치를 그대로 쓰려면 백업 tar.gz를 풀어서 쓰고, 갱신 효과를 보려면 새로 받으십시오.

## 2. 환경

```bash
conda env create -f exp/FLOOR/headroom/environment.yml
conda activate headroom
```

재도킹까지 다시 하려면 별도로 AutoDock Vina 1.2.7, Meeko, Open Babel, sPyRMSD,
PoseBusters, ODDT가 필요합니다. 이 저장소의 `repro/REPRODUCE.md`에 그 환경이 적혀 있습니다.

## 3. 네트워크 허용 목록

새 PC가 같은 샌드박스를 쓴다면 다음이 열려 있어야 합니다.

| 도메인 | 쓰임 | 상태 |
| --- | --- | --- |
| `data.rcsb.org` | PDB→UniProt 매핑, 구조·리간드 | 기본 허용 |
| `www.bindingdb.org` | 친화도 수집 | **요청해서 승인받음** |
| `www.ebi.ac.uk` | ChEMBL, Europe PMC | 기본 허용 |
| `pubchem.ncbi.nlm.nih.gov` | 물성·부분구조 검색 | 기본 허용 |
| `pharos-api.ncats.io` | 타깃 성숙도(TDL) | **요청해서 승인받음** |
| `zinc20/22.docking.org` | 승인은 받았으나 **사이트가 봇 차단** | 사실상 사용 불가 |
| `pdbbind-cn.org`, `zenodo.org`, `huggingface.co`, `dataverse.harvard.edu`, `tdcommons.ai` | PDBbind·TDC 데이터 | **차단** |

## 4. 알려진 함정 (여기서 겪은 것들)

1. **PDBbind General은 스크립트로 못 받습니다.** 등록이 필요하고 미러는 차단돼 있습니다.
   RCSB도 대체재가 아닙니다 — 무작위 24개 중 친화도 주석 9건, Ki/Kd는 2건뿐입니다.
   그래서 외부 학습 바닥은 BindingDB로 대체했습니다(리간드 전용이라 구조가 불필요).
2. **프로세스 풀이 막힐 수 있습니다.** 잠긴 컨테이너는 `ProcessPoolExecutor`가 쓰는
   네임드 파이프를 거부합니다(Windows error 5). `featurize.py`가 감지해 직렬로 떨어지고,
   RDKit이 내부적으로 스레딩하므로 39,000 분자가 약 13분입니다.
3. **BindingDB의 PDB 키 엔드포인트는 HTTP 500**을 돌려줍니다. UniProt 키를 쓰고
   RCSB로 PDB→UniProt을 매핑하십시오.
4. **응답 최상위 키에 오타가 있습니다** — `getLindsByUniprotsResponse`. 이름으로 찾지 말고
   단일 값을 읽으십시오.
5. **`Ipc` 기술자가 float32 범위를 넘칩니다.** 학습 전에 비유한·거대값 컬럼을 걸러내야
   scikit-learn이 터지지 않습니다.
6. **한글 그림은 `Malgun Gothic`을 명시**해야 합니다. 안 하면 글리프가 모두 네모로 나옵니다.
   이 폰트에는 위첨자 글리프가 없어 지수는 mathtext로 써야 합니다.
7. **ODDT는 Windows에서 세 가지 우회가 필요합니다** — `pip install oddt --no-deps`,
   import 전에 `np.in1d = np.isin`, RDKit 백엔드 사용.

## 5. 실행 순서

```bash
cd exp/FLOOR/headroom
python 01_pull.py      --index casf_index.csv --work work
python featurize.py    --in casf_smiles.csv --out casf_features.csv
python 02_ceiling.py   --index casf_index.csv --smiles casf_smiles.csv --work work
python 03_floor.py     --index casf_index.csv --smiles casf_smiles.csv \
                       --ligand-features casf_features.csv --pockets pockets/ --work work
python featurize.py    --in work/bindingdb_train.csv --out work/bindingdb_features.csv
python 04_external_floor.py --index casf_index.csv --smiles casf_smiles.csv \
                       --casf-uniprot work/casf_uniprot.csv --matched work/matched_records.csv \
                       --train-features work/bindingdb_features.csv \
                       --test-features casf_features.csv --work work
```

대략 01은 2분, 02는 2분, 03은 1분, featurize는 13분, 04는 중첩 체제 2분 + 타깃 분리 체제
40분 정도입니다.

## 6. 재현 성공 판정

아래 값이 나오면 이관이 제대로 된 것입니다(BindingDB를 새로 받으면 천장 쪽은
소수점 셋째 자리가 흔들릴 수 있습니다).

| 값 | 기대 |
| --- | --- |
| 복합체 수 / 타깃 수 | 285 / 57 |
| 천장 (복제 기반) | 0.88 ± 0.01 |
| 바닥 리간드+포켓, 타깃 분리 | 0.534 ± 0.01 |
| 바닥 리간드+포켓, 무작위 분할 | 0.685 ± 0.01 |
| 외부 학습(BindingDB), 타깃 중첩 | 0.679 ± 0.01 |
| Vina / Vinardo / RF-Score v2 | 0.502 / 0.428 / 0.741 |

## 7. 아직 열려 있는 것

- **TDC(Therapeutics Data Commons)** — 데이터 호스트가 차단돼 있어 아직 못 봤습니다.
  허용되면 독립 라벨 출처로 천장 추정을 교차검증할 수 있습니다.
- **PDBbind General** — 확보되면 외부 학습 바닥을 포켓 기준선까지 확장할 수 있습니다.
- **GPU 리눅스 장비** — bioemu2 앙상블 재채점, Boltz-2 공동접힘 비교가 여기 묶여 있습니다.
