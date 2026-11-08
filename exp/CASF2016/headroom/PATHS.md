# 파일 위치 루트

<!-- exp/CASF2016/headroom/PATHS.md · 2026-11-08 -->

HEADROOM의 파일은 **네 군데**에 있습니다. 이 중 **영속적인 것은 GitHub 하나뿐**입니다.
나머지 셋은 이 PC·이 대화에 묶여 있어서, 다른 공간에서 이어받을 때는 ①만 보면 됩니다.

---

## ① GitHub — 유일한 영속 루트 ★

```
github.com/Leokong-C/joi-bio
  └─ branch: headroom/pipeline        (PR #6, 25 files)
       └─ exp/CASF2016/headroom/      ← 모든 것의 루트
```

| 경로 | 내용 |
| --- | --- |
| `README.md` | 전체 기록과 로그 (21.1 KB) |
| `SPEC.md` / `BACKEND.md` / `MIGRATION.md` | 설계 · 백엔드 규약 · 이전 계획 |
| `PARTS_DRAFT.md` / `PRIOR_ART.md` | 게시물 초안 · 선행연구 |
| `VERSIONS.json` / `environment.yml` | 고정 버전 · 환경 |
| `01_pull.py` `02_ceiling.py` `03_floor.py` `04_external_floor.py` `featurize.py` | 파이프라인 5단계 |
| `results/*.csv` `results/summary.json` | 바닥 · 천장 · 외부학습 · 분할효과 · 감사 · DB버전 |
| `figures/band.png` `figures/pipeline.png` | 도식 2장 |
| `data/casf_inputs.tar.gz` (151 KB) | SMILES · UniProt 매핑 · 리간드 키 · SDF |
| `data/pockets.tar.gz` (6.8 MB) | 16 Å 포켓 PDB 285개 (원본 29.3 MB) |
| `data/bindingdb_raw.tar.gz` (3.7 MB) | BindingDB 원자료 78개, 2026-10-25 수집 (원본 75 MB) |

**다른 PC에서 시작하는 방법**

```bash
git clone https://github.com/Leokong-C/joi-bio
cd joi-bio && git checkout headroom/pipeline
cd exp/CASF2016/headroom
conda env create -f environment.yml && conda activate headroom
python -c "import tarfile,glob; [tarfile.open(f).extractall('.') for f in glob.glob('data/*.tar.gz')]"
python featurize.py && python 03_floor.py
```

`casf_index.csv`(CoreSet.dat 유래)만 저장소에 없습니다 — **재배포 가능 여부 확인 전까지 손으로 옮깁니다.**

---

## ② 로컬 작업공간 — 임시, 쓸려 나감 ⚠

```
C:\Users\sso03\.claude-science\orgs\03c00401-07e6-4964-a511-acfb42a28cca
   \workspaces\b27db907-6f1e-40d9-a4a3-8ee6b95582f7\
```

현재 **634.6 MB**. 유휴 상태가 길어지면 **삭제됩니다.** 여기만 있고 ①에 없는 것은 잃습니다.

| 하위 경로 | 크기 | ①에 있나 |
| --- | --- | --- |
| `casf/` (863) | 137.2 MB | 포켓만 번들로 있음 |
| `dock/` (1,264) | 136.9 MB | **없음** — 결과 CSV만 보관 |
| `headroom_bindingdb_features.csv` | 78.5 MB | **없음** — 13분이면 재생성 |
| `bindingdb/` (78) | 75.0 MB | 압축본으로 있음 |
| `mmgbsa/` (56) | 39.5 MB | 결과 CSV만 |
| `tu_staging/` (1,635) | 39.2 MB | ToolUniverse 스킬 — 발행 완료 |
| `RFScore_v{1,2,3}_pdbbind2016.pickle` | 58.3 MB | **없음** — 재학습 필요 |
| `backup_headroom_20261107.tar.gz` | 27.6 MB | 아티팩트로 보관 |
| `data/` | 11.2 MB | ①에 그대로 |

도킹 번들과 RF-Score 피클을 뺀 것은 의도입니다 — 전자는 재도킹이 시드에 따라 달라져 결과 CSV를 고정하는 편이 낫고, 후자는 GitHub 경고선을 넘습니다.

---

## ③ 아티팩트 저장소 — 대화에 묶임

프로젝트 `proj_90d9aa821ba1`. 그림·CSV·JSON·백업 압축본이 버전 단위로 들어 있고,
코드에서는 `host.artifact_path(version_id)`로 경로처럼 씁니다. 파일 시스템에는 보이지 않습니다.
**핵심 산출물은 전부 ①에도 넣어 두었습니다.**

## ④ Notion — 서사 기록

「AI 도킹·스코어링 벤치마크 (CASF-2016)」 · page_id `3faed537-978c-814b-82bc-eb6b2cb982b4` (44,269자).
파트 1–8이 "무엇을 발견했나", 「📒 운영 기록」이 "어떻게 돌리고 무엇을 조심하나"로 README와 대응합니다.

---

## 대화 기록 원문

```
C:\Users\sso03\.claude-science\orgs\03c00401-07e6-4964-a511-acfb42a28cca
   \tool-results\b27db907-6f1e-40d9-a4a3-8ee6b95582f7\compaction-history\
   b27db907-6f1e-40d9-a4a3-8ee6b95582f7\compaction_*.txt
```

## 한 줄 요약

**이어받을 때 보는 곳은 `github.com/Leokong-C/joi-bio` → `headroom/pipeline` 브랜치 → `exp/CASF2016/headroom/` 하나입니다.** 나머지 셋은 이 PC·이 대화의 사본입니다.
