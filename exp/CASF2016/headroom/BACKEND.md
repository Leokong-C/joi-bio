# GitHub을 백그라운드로 쓰는 설계

"파이프라인은 그대로 두고 구동 PC만 바꾼다"를 성립시키려면, 저장소가 **무엇을 담고 무엇을
담지 않는지**가 분명해야 합니다. 용량을 재고 정한 결과입니다.

## 담는 것 (총 약 12 MB)

| 경로 | 내용 | 크기 |
| --- | --- | --- |
| `headroom/*.py` | 4단계 스크립트 + featurize | 20 KB |
| `headroom/*.md` | README · SPEC · MIGRATION · PARTS_DRAFT · PRIOR_ART · BACKEND | 50 KB |
| `headroom/VERSIONS.json` | 수치를 좌우하는 패키지 버전 고정 | 1 KB |
| `headroom/results/*.csv,json` | 바닥·천장·분할효과·외부학습·검증·DB버전 | 60 KB |
| `headroom/figures/*.png` | 파이프라인 도식, 바닥–천장 띠 | 330 KB |
| `headroom/data/casf_inputs.tar.gz` | SMILES · UniProt 매핑 · 리간드 키 · SDF | 0.2 MB |
| `headroom/data/pockets.tar.gz` | 16 Å 포켓 PDB 285개 | 7.1 MB |
| `headroom/data/bindingdb_raw.tar.gz` | BindingDB 원자료 78개 (수집일 2026-10-25) | 3.9 MB |

## 담지 않는 것 — 그리고 그 이유

| 항목 | 원래 크기 | 왜 안 담나 |
| --- | --- | --- |
| `headroom_bindingdb_features.csv` | **78.5 MB** | GitHub 경고선(50 MB) 초과. `featurize.py`로 13분이면 재생성 |
| 도킹 번들 (PDBQT, 포즈) | 28 MB+ | 재도킹은 어차피 비결정적이라 **결과 CSV만** 담는 것이 맞음 |
| `casf_index.csv` (CoreSet.dat 유래) | 작음 | **라이선스 확인 필요** — PDBbind 유래 데이터라 저장소를 공개로 돌리기 전에 재배포 가능 여부를 확인해야 함. 현재 private이므로 내부 사용은 무방 |

## 이 설계로 해결되는 것

1. **도킹 비결정성이 사라집니다.** 재도킹 결과 CSV를 저장소에 고정해 두면 새 PC는 재도킹하지
   않고 그 결과를 읽습니다. 시드마다 달라지던 66.3 %가 더 이상 흔들리지 않습니다.
   재도킹을 *다시 하고 싶을 때만* 하고, 그때는 값이 달라진다는 것을 알고 하는 것입니다.
2. **BindingDB 갱신에 휘둘리지 않습니다.** `bindingdb_raw.tar.gz`(수집일 2026-10-25)를
   풀어 쓰면 천장 수치가 고정됩니다. 갱신 효과를 보고 싶으면 `01_pull.py`를 돌리면 되고,
   그때는 수치가 달라진다는 것을 알고 하는 것입니다.
3. **수작업 정정이 보존됩니다.** 3syr→GPK, 1a30 펩타이드 같은 정정이 `casf_inputs.tar.gz`
   안에 들어 있어, 재생성으로 날아가지 않습니다.

## 이 설계로도 해결되지 않는 것 — 새 PC에서 반드시 해야 할 네 가지

1. **네트워크 허용 목록은 설치마다 다릅니다.** `www.bindingdb.org`, `pharos-api.ncats.io`,
   `dataverse.harvard.edu` 승인은 저장소에 들어가지 않습니다. 새 환경에서 다시 요청해야 합니다.
   단, 위 설계대로 원자료를 저장소에 담아 두면 **분석만 할 때는 네트워크가 아예 불필요합니다.**
2. **패키지 버전을 맞춰야 합니다.** `environment.yml`은 버전을 고정하지 않습니다.
   `VERSIONS.json`에 적힌 조합(RDKit 2026.09.1 · scikit-learn 1.9.1 · NumPy 2.4.6)과
   다르면 랜덤포레스트의 동점 처리나 기술자 정의가 달라져 소수점 셋째 자리가 움직일 수 있습니다.
3. **OS 차이가 남습니다.** 한글 그림은 Windows의 `Malgun Gothic`을 쓰고 있어 Linux·macOS에서는
   폰트를 바꿔야 합니다. 프로세스 풀 가용성도 환경마다 다르고(여기서는 거부됨), ODDT는
   Windows에서 세 가지 우회가 필요합니다.
4. **재도킹까지 할 거라면** Vina 1.2.7 · Meeko · Open Babel · sPyRMSD · PoseBusters가 따로
   필요하고, 그때는 결과가 ±1 pt 흔들립니다.

## 결론

**분석 파이프라인에 대해서는 맞습니다** — 저장소를 백그라운드로 두면 PC가 바뀌어도 같은 수치가
나옵니다. 방금 저장된 특징 파일만으로 바닥 6개 수치를 재계산해 소수점 셋째 자리까지 일치하는 것을
확인했습니다.

**재도킹에 대해서는 아닙니다** — 다만 결과를 저장소에 고정하면 그 문제를 *우회*할 수 있고,
이 설계가 바로 그것입니다.
