# 다른 컴퓨터로 옮기기

<!-- exp/CASF2016/headroom/PORTING.md · 2026-11-08 -->

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
cd exp/CASF2016/headroom
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

**pip이 Windows ACL 오류(`WinError 5 ... source_docs`)로 깨지면** 환경을 바꿔도 소용없습니다.
순수 파이썬 패키지는 PyPI wheel을 직접 받아 `vendor/`에 풀고 `sys.path`에 얹으십시오
(`PLM.md` §1에 코드). wheel의 SHA-256은 PyPI가 같이 줍니다 — 반드시 대조하십시오.

## 4. 손으로 옮겨야 하는 것 — 두 가지뿐

| 파일 | 왜 저장소에 없나 | 어떻게 |
| --- | --- | --- |
| `casf_index.csv` | CASF-2016 `CoreSet.dat` 유래, **재배포 가능 여부 미확인** | 기존 PC에서 복사 |
| `articles/` 논문 PDF 4편 | 출판사 저작권 | DOI는 `PRIOR_ART.md`에 있음, 각자 내려받기 |

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

## 8. 무엇이 이전되지 않는가

- **대화 기록과 아티팩트 저장소** — 이 프로젝트에 묶여 있습니다. 산출물은 전부 저장소에도 있습니다.
- **Notion · Slack · Gmail 커넥터** — 새 환경에서 다시 연결해야 합니다.
- **ToolUniverse 스킬** — 발행되어 있으므로 다시 만들 필요는 없습니다.
