# 다른 컴퓨터로 옮기기

<!-- exp/FLOOR/headroom/PORTING.md · 2026-11-08 -->

HEADROOM은 **GitHub 한 곳에서 전부 복원됩니다.** 로컬 작업공간은 사본이고 유휴 시 삭제됩니다.
이 문서는 새 PC에서 0부터 같은 수치를 내는 절차입니다.

## 0. 필요한 것

| 항목 | 최소 | 비고 |
| --- | --- | --- |
| OS | Windows / Linux / macOS | 수치 재현은 Windows에서 검증 |
| 코어 | 8 | 16이면 ESM-2 임베딩 191초 |
| RAM | 16 GB | ESM-2 650M 로드에 ~4 GB |
| 디스크 | 15 GB | 번들 해제 585 MB + 가중치 2.6 GB + 환경 |
| GPU | 불필요 | 바닥·천장·ESM-2 전부 CPU |
| 저장소 접근 | `Leokong-C/joi-bio` 읽기 | **비공개 저장소** |

## 1. 받기

```bash
git clone https://github.com/Leokong-C/joi-bio
cd joi-bio && git checkout headroom/pipeline
cd exp/FLOOR/headroom
```

## 2. 번들 복원 — 해시로 확인하고 풀기

```bash
cd bundles
# 분할 합치기 (RF-Score 묶음만 2분할)
#   Linux/macOS
cat rfscore.tar.gz.part00 rfscore.tar.gz.part01 > rfscore.tar.gz
#   Windows — PowerShell의 Get-Content 는 바이너리를 망칩니다
cmd /c copy /b rfscore.tar.gz.part00+rfscore.tar.gz.part01 rfscore.tar.gz
```

```python
import hashlib, json, os, tarfile, glob
man = {m["file"]: m for m in json.load(open("MANIFEST.json"))}
for f in sorted(man):
    if os.path.exists(f):
        h = hashlib.sha256(open(f, "rb").read()).hexdigest()
        print(("OK  " if h == man[f]["sha256"] else "FAIL"), f)
for f in glob.glob("*.tar.gz"):
    tarfile.open(f).extractall("../../../..")   # 저장소 루트로 풀림
```

11개 전부 OK 여야 합니다. 해제하면 4,664개 멤버, 585 MB.

## 3. 환경

```bash
conda env create -f environment.yml && conda activate headroom
```

`environment.yml`은 버전을 고정하지 않습니다. 수치를 소수점 셋째 자리까지 맞추려면
`VERSIONS.json`의 조합을 쓰십시오 — Python 3.11.17 · RDKit 2026.09.1 · NumPy 2.4.6 ·
SciPy 1.17.1 · pandas 3.0.6 · scikit-learn 1.9.1 · matplotlib 3.11.2.
다르면 RandomForest 동점 처리와 RDKit 기술자 정의가 달라져 셋째 자리가 움직입니다.

****실측 (2026-11-11).** 고정 버전을 지키지 **않고** 바닥 단계를 돌려도 세 수치가
소수점 셋째 자리까지 그대로 나왔습니다(scikit-learn 1.9.1 · NumPy 2.5.3 ·
pandas 3.0.5 · SciPy 1.18.1 · Python 3.13). 맞춘 것은 scikit-learn 하나뿐입니다.
RandomForest 동점 처리가 거기 있으므로 **1.9.1만 먼저 맞추고**, 수치가 어긋나면
그때 나머지를 맞추는 순서가 빠릅니다. RDKit은 바닥 단계에 불필요합니다 — 특징이
번들에 이미 계산돼 있고, `featurize.py`를 새로 돌릴 때만 필요합니다.

pip이 Windows ACL 오류(`WinError 5 ... source_docs`)로 깨지면** 환경을 바꿔도 소용없습니다.
순수 파이썬 패키지는 PyPI wheel을 직접 받아 `vendor/`에 풀고 `sys.path`에 얹으십시오
(`PLM.md` §1에 코드). wheel의 SHA-256은 PyPI가 같이 줍니다 — 반드시 대조하십시오.

## 4. 손으로 옮겨야 하는 것 — 한 가지뿐

| 파일 | 왜 저장소에 없나 | 어떻게 |
| --- | --- | --- |
| `articles/` 논문 PDF 4편 | 출판사 저작권 | DOI는 `PRIOR_ART.md`에 있음, 각자 내려받기 |

**정정 (2026-11-11).** 이 절은 `casf_index.csv`도 손으로 옮기라고 적고 있었으나
**틀렸습니다.** 그 파일은 번들 안에 들어 있습니다 — `bundles/casf.tar.gz`를 풀면
`casf/index.csv`로 나옵니다(285행 · 57 타깃). 같은 데이터가
`stage1-4/casf2016_core_index.csv`에도 있고, 두 벌을 대조한
결과 `code`·`logKa`·`target`·`resolution`이 전부 일치합니다(부동소수점 표기만 다름).
`03_floor.py`가 요구하는 컬럼이 정확히 이 셋입니다.

`casf_smiles.csv`, `casf2016_predictor_matrix.csv`, `pockets/` 285개도 마찬가지로
번들에 들어 있습니다(각 285행). 별도 기록에서 "반드시 옮겨야 하는 것"으로
지목됐던 네 파일 전부 해당 없습니다.

## 5. 돌리기 — 비용과 기대값

| 단계 | 명령 | 시간 | 결과 |
| --- | --- | ---: | --- |
| 특징 생성 | `python featurize.py` | 13분 | 리간드 205 · 포켓 21 |
| **바닥** | `python 03_floor.py` | 1분 | 리간드+포켓 분리 **0.534** · 노출 **0.685** |
| 천장 | `python 02_ceiling.py` | 2분 | 복제 기반 **0.88** |
| BindingDB 수집 | `python 01_pull.py` | 5분 | 78 파일 · 258,749 레코드 |
| 외부 학습 | `python 04_external_floor.py` | 46분 | 중첩 **0.679** · 분리 **0.551** |
| ESM-2 임베딩 | `python plm_embed.py` | 3분 + 가중치 2.6 GB 내려받기 | 285 × 1280 |
| ESM-2 바닥 | `python plm_floor.py` | 6분 | 개수 **0.241** vs ESM-2 **0.179** |

**번들 특징 파일은 키가 다릅니다.** `headroom_features_ligand.csv`는 `code` 키인데
`03_floor.py`는 `smiles` 키를 기대하므로 그대로 넘기면
`KeyError: "None of ['smiles'] are in the columns"`로 죽습니다. `casf_smiles.csv`로
join해 다시 키를 잡으십시오 — 285행이 고유 SMILES 279개로 줄어드는 것이 정상입니다.

```python
import pandas as pd
f = pd.read_csv("headroom_features_ligand.csv")
s = pd.read_csv("casf_smiles.csv")
m = f.merge(s, on="code").drop(columns=["code"]).drop_duplicates(subset=["smiles"])
m[["smiles"] + [c for c in m.columns if c != "smiles"]].to_csv("features_by_smiles.csv", index=False)
```

재도킹(`dock/`)은 **다시 돌리지 마십시오.** 34.8 CPU-시간이 들고 시드가 바뀌면
성공률이 SD 0.41 pt로 흔들리며 27건이 뒤집힙니다. 번들의 결과 CSV가 고정된 기준입니다.

## 6. 재현 성공 판정

| 값 | 기대 |
| --- | --- |
| 복합체 / 타깃 | 285 / 57 |
| 천장(복제 기반) | 0.88 ± 0.01 |
| 바닥 리간드+포켓, 타깃 분리 / 노출 | 0.534 / 0.685 ± 0.01 |
| 외부 학습, 중첩 / 분리 | 0.679 / 0.551 ± 0.01 |
| Vina / Vinardo / RF-Score v2 | 0.502 / 0.428 / 0.741 |
| ESM-2 포켓 / 아미노산 개수 (분리) | 0.179 / 0.241 |

하나라도 ±0.01을 넘게 벗어나면 **버전부터 의심하십시오.** `VERSIONS.json`과 대조가 먼저입니다.

## 7. 네트워크 — 허용 목록이 필요한 호스트

| 호스트 | 쓰는 단계 |
| --- | --- |
| `data.rcsb.org` | PDB → UniProt 매핑 |
| `www.bindingdb.org` | 천장·외부 학습 |
| `dl.fbaipublicfiles.com` | ESM-2 가중치 |
| `dataverse.harvard.edu` + `dvn-cloud-iqss.s3.amazonaws.com` | TDC (미실행) |
| `api.github.com` | 저장소 |

막혀 있는 채로 남은 것 — `huggingface.co`(ESMFold2 · Boltz-2 가중치),
`api.colabfold.com`(MSA 서버), `www.pdbbind-cn.org`(PDBbind General).
앞의 둘은 **GPU가 없으면 열어도 쓸 수 없습니다.**

## 7-1. 이 브랜치와 benchmark 브랜치의 관계 (2026-11-11 확인)

CASF-2016 작업은 **한 줄기**입니다. git 이력만 끊겨 있습니다.

| 브랜치 | PR | 관계 |
| --- | --- | --- |
| `docking-benchmark-casf2016` | #1 (열림) | 1~2단계 |
| `stage4-extended-benchmark` | #3 (열림, base가 `main`이 아니라 위 브랜치) | 3~4단계, 위의 직계 자손 |
| `headroom/pipeline` | #6 (열림) | 5단계 이후. **위 둘의 자손이 아니라** `main`의 다른 지점에서 갈라짐 |

**산출물은 그대로 넘어왔습니다.** stage4 브랜치의
`stage1-4/casf2016_docking_results.csv`와 이 브랜치
`bundles/dock.tar.gz` 속 같은 파일이 **sha256까지 바이트 단위로 동일**합니다
(`e3e82e73…`, 284행). MM/GBSA · 포즈 재순위 · RF-Score · LLM 결과도 `mmgbsa` ·
`rerank` · `rfscore` · `tables` 번들에 실려 있습니다. 따라서 **이 번들만 복원하면
1~4단계 결과까지 따라옵니다** — benchmark 브랜치를 따로 받을 필요가 없습니다.
재도킹을 다시 돌리지 말라는 §5의 경고가 성립하는 이유도 이것입니다.

**다만 문서는 넘어오지 않았습니다.** `STAGE4.md`, `PROJECT_DIRECTION.md`,
`repro/REPRODUCE.md`, 그림 8장은 `stage4-extended-benchmark`에만 있습니다. 특히
재현 키트의 경고 — **Open Babel의 극성 수소 배치가 프로세스마다 달라 비트 수준
재현에는 준비된 PDBQT 보관이 필요하다** — 는 이 문서에 없는 내용입니다.

`main`에는 셋 중 어느 것도 병합돼 있지 않습니다. #1·#3·#4·#5·#6 전부 열린 상태입니다.

## 7-2. 재현성 — 두 가지를 구분하십시오 (stage1-4/repro/REPRODUCE.md 흡수)

| 무엇 | 재현되는가 |
| --- | --- |
| 바닥·천장·ESM-2 수치 | **예** — 보관된 특징·결과에서 재계산하면 셋째 자리까지 일치 |
| 도킹 점수 비트 수준 | **아니오** — 아래 |

**Open Babel은 극성 수소를 프로세스마다 다르게 배치합니다.** 한 수용체 PDBQT
4,190줄 중 **26줄**이 두 번의 준비 사이에서 달라졌습니다(Ser/Cys에 붙는 H).
보관된 수용체로 재채점하면 −8.479 kcal/mol이 나오지만, 새로 준비한 수용체는
−8.355를 냅니다. **Vina 자체는 결정적입니다** — 비결정성은 전처리에 있습니다.

따라서 **비트 수준 재현에는 파이프라인이 아니라 준비된 PDBQT 파일의 보관이
필요합니다.** 번들의 `dock/prep`이 그 보관본입니다. 재도킹을 다시 돌리지 말라는
§5의 경고와 같은 뿌리입니다.

시드를 바꾸면 성공률이 SD 0.41 pt로 흔들리고 27건의 판정이 뒤집힙니다.

전체 재현 키트(정확한 패키지 버전, Vina 실행파일 해시, 실측 소요시간)는
`stage1-4/repro/REPRODUCE.md`에 있습니다.

## 8. 무엇이 이전되지 않는가

- **대화 기록과 아티팩트 저장소** — 이 프로젝트에 묶여 있습니다. 산출물은 전부 저장소에도 있습니다.
- **Notion · Slack · Gmail 커넥터** — 새 환경에서 다시 연결해야 합니다.
- **ToolUniverse 스킬** — 발행되어 있으므로 다시 만들 필요는 없습니다.
