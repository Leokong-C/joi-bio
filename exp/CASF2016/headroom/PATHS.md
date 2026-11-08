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

---

## 번들 — 작업공간 전체 (2026-11-08)

`exp/CASF2016/headroom/bundles/`. 작업공간 전체를 압축해 넣었습니다. **634.6 MB → 203.1 MB.**
이제 로컬에만 있는 것은 없습니다.

| 묶음 | 압축 | 원본 | 멤버 | 내용 |
| --- | ---: | ---: | ---: | --- |
| `casf.tar.gz` | 33.9 MB | 137.2 MB | 869 | 구조 · 포켓 · 트리밍 수용체 |
| `rfscore.tar.gz.part{00,01}` | 58.5 MB | 58.4 MB | 3 | RF-Score v1/v2/v3 학습 모델 |
| `dock.tar.gz` | 28.8 MB | 136.9 MB | 1,271 | Vina 도킹 포즈 전량 |
| `tables.tar.gz` | 23.9 MB | 85.0 MB | 62 | 결과 CSV · JSON 전부 |
| `tu_staging.tar.gz` | 23.5 MB | 39.2 MB | 2,300 | ToolUniverse 스킬 작업물 |
| `bdb_features.tar.gz` | 23.2 MB | 83.0 MB | 2 | BindingDB 특징 행렬 · 학습셋 |
| `mmgbsa.tar.gz` | 8.3 MB | 39.5 MB | 58 | MM/GBSA 궤적 |
| `push.tar.gz` | 2.3 MB | 2.7 MB | 32 | 게시용 그림 · 문서 |
| `rerank.tar.gz` | 0.2 MB | 1.6 MB | 65 | 포즈 재순위 산출물 |
| `bin.tar.gz` | 0.5 MB | 1.2 MB | 2 | AutoDock Vina 1.2.7 실행파일 (Apache-2.0) |

`MANIFEST.json`에 파일별 SHA-256이 있습니다. **복원됐는지를 말이 아니라 해시로 확인하십시오.**

### 복원

```bash
cd exp/CASF2016/headroom/bundles
cat rfscore.tar.gz.part00 rfscore.tar.gz.part01 > rfscore.tar.gz   # 분할 합치기
python - <<'EOF'
import hashlib, json, os, tarfile, glob
man = {m["file"]: m for m in json.load(open("MANIFEST.json"))}
for f in sorted(man):
    if os.path.exists(f):
        h = hashlib.sha256(open(f,"rb").read()).hexdigest()
        print(("OK  " if h == man[f]["sha256"] else "FAIL"), f)
for f in glob.glob("*.tar.gz"):
    tarfile.open(f).extractall("../../../..")      # 저장소 루트로 풀림
EOF
```

Windows PowerShell에서 분할을 합칠 때는 `cmd /c copy /b rfscore.tar.gz.part00+rfscore.tar.gz.part01 rfscore.tar.gz`를 씁니다 — PowerShell의 `Get-Content`는 바이너리를 망칩니다.

**무결성 확인(2026-11-08, 커밋 `5f58203f`)** — 두 단계로 검사했습니다.

1. 묶음 11개를 전부 열어 모든 멤버를 끝까지 읽어 CRC 검사 — **4,664개 멤버, 손상 0건**, 분할본도 합쳐서 정상 해제.
2. GitHub에서 12개 블롭을 **다시 내려받아 SHA-256과 바이트 길이를 대조 — 11/11 일치, 불일치 0건.**

2단계가 핵심입니다. 올렸다는 말이 아니라, 내려받아 복원되는 것을 해시로 확인한 기록입니다.

### 저장소에 없는 것 — 하나뿐입니다

`articles/` 7.0 MB — 출판된 논문 PDF 4편(`10.1093_bioinformatics_btaf040`, `10.1371_journal.pone.0061007`,
`10.1038_s41746-026-02477-w`, `10.1039_d6sc01189f`). 비공개 저장소라도 **출판사 저작물을 재배포하는 형태**가
되므로 넣지 않습니다. DOI는 `PRIOR_ART.md`에 있으니 각자 내려받으면 됩니다.

`casf.tar.gz`에는 CASF-2016 구조가 들어 있습니다. **이 저장소를 비공개로 유지하는 것이 전제**이며,
`casf_index.csv`(CoreSet.dat 유래)의 재배포 가능 여부 확인은 여전히 미결입니다.
