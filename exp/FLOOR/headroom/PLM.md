# 단백질 모델 축 — 설치 재검토와 ESM-2 결과

<!-- exp/FLOOR/headroom/PLM.md · 2026-11-08 -->

AlphaFold2 · ESMFold2 · Boltz-2 · ESM-2 네 개를 검토해 **실제로 돌아가는 하나만** 돌렸습니다.

## 1. 무엇이 왜 막혔나

막은 것은 GPU가 아니라 **가중치 호스트**였습니다. 코드는 전부 받아집니다.

| 호스트 | 쓰임 | 상태 |
| --- | --- | --- |
| `pypi.org` | boltz · fair-esm 코드 | 200 |
| `dl.fbaipublicfiles.com` | ESM-2 가중치 | 차단 → **허용 받음** |
| `huggingface.co` | ESMFold2 · Boltz-2 가중치 | 차단 |
| `api.colabfold.com` | MSA 서버 (AF2 · Boltz 필수) | 차단 |

여기에 더해 이 PC에는 GPU가 없고 원격 compute 대상이 0개입니다.
AF2 · ESMFold2 · Boltz-2는 셋 다 GPU 전제라 **가중치가 열려도 돌릴 하드웨어가 없습니다.**
GPU 리눅스 확보가 세 개를 한꺼번에 여는 단일 관문입니다.

### 설치기 자체가 고장나 있습니다

pip 설치가 **모든 conda 환경에서** 같은 Windows ACL 오류로 실패합니다:

```
ERROR: Could not install packages due to an OSError:
[WinError 5] 액세스가 거부되었습니다: 'C:\...\conda\envs\<env>\source_docs'
```

`plm`(새로 만든 환경) · `docking` · `md` 세 곳에서 동일하게 재현됩니다. 환경을 바꿔서 될 일이 아닙니다.

**우회** — `md` 환경에 torch 2.14.1이 이미 있었고, `fair-esm`은 순수 파이썬이라
PyPI wheel을 받아 작업공간에 풀어 `sys.path`에 얹었습니다. wheel의 SHA-256을 PyPI가
게시한 값과 대조해 확인했습니다(fair-esm 2.0.0, 93,125 B).

```python
import urllib.request, json, zipfile, io, hashlib, sys, os
meta = json.loads(urllib.request.urlopen("https://pypi.org/pypi/fair-esm/json").read())
whl  = [u for u in meta["urls"] if u["packagetype"] == "bdist_wheel"][0]
data = urllib.request.urlopen(whl["url"]).read()
assert hashlib.sha256(data).hexdigest() == whl["digests"]["sha256"]
zipfile.ZipFile(io.BytesIO(data)).extractall("vendor")
sys.path.insert(0, os.path.abspath("vendor"))
```

ESM-2 650M 가중치(2,604.5 MB)와 contact-regression 파일은 스킬이 고정한 SHA-256으로
검증해 받았습니다. 둘 다 통과.

## 2. 질문 — 단백질 언어모델이 바닥을 올리는가

HEADROOM의 바닥 모델에서 단백질 쪽 표현만 바꿨습니다.
**아미노산 개수 21개 → ESM-2 650M 평균풀링 임베딩 1280차원.**
나머지는 전부 동일합니다: 285개 복합체, RandomForest 500그루, GroupKFold-10(타깃별),
포켓은 리간드 16 Å 이내(중앙값 152잔기). 임베딩은 CPU 16스레드로 191초.

## 3. 답 — 아니오. 오히려 내려갑니다

| 단백질 표현 | 차원 | 타깃 분리 | 타깃 노출 |
| --- | ---: | ---: | ---: |
| **아미노산 개수** | 21 | **0.241** | **0.509** |
| ESM-2 650M 임베딩 | 1280 | **0.179** | 0.371 |
| ESM-2 + 아미노산 개수 | 1301 | 0.191 | 0.386 |

리간드와 합쳐도 같습니다.

| 조합 | 차원 | 타깃 분리 | 타깃 노출 |
| --- | ---: | ---: | ---: |
| 리간드만 | 205 | 0.502 | 0.641 |
| **리간드 + 아미노산 개수** | 226 | **0.534** | **0.685** |
| 리간드 + ESM-2 | 1485 | **0.439** | 0.608 |

### 차원 때문이 아닙니다

"285샘플에 1280차원을 넣었으니 당연하다"는 반론을 먼저 막았습니다. 줄여도 더 떨어집니다.

| 축소 | 타깃 분리 | 타깃 노출 |
| --- | ---: | ---: |
| PCA 64 | 0.129 | 0.370 |
| PCA 32 | 0.129 | 0.361 |
| PCA 16 | 0.090 | 0.332 |
| PCA 8 | 0.117 | 0.324 |
| Ridge(선형, 1280차원) | 0.112 | 0.330 |
| 참고 — Ridge(아미노산 개수 21) | 0.195 | 0.359 |

ESM-2의 최고값은 **0.179**(축소 없는 RF)이고, 아미노산을 세기만 한 0.241에 못 미칩니다.

## 4. 읽는 법 — 그리고 과장하지 않기

이 결과가 말하는 것은 하나입니다. **포켓 서열의 표현을 더 좋게 만들어도 바닥은 올라가지 않습니다.**
리간드 기술자를 14개에서 205개로 늘렸을 때 바닥이 0.547 → 0.534로 움직이지 않았던 것과 같은 방향이며,
표현력이 아니라 **문제 구조**가 바닥을 정한다는 증거가 하나 더 늘었습니다.

이 결과가 말하지 **않는** 것:

- "ESM-2에 정보가 없다" — 아닙니다. 측정한 것은 *평균풀링 임베딩 + RandomForest + n=285 + 타깃 분리 분할*이라는
  특정 조합입니다. 잔기별 임베딩, 접촉 예측, 미세조정은 다른 이야기이며 여기서 재지 않았습니다.
- "단백질 언어모델이 도킹에 쓸모없다" — 재지 않았습니다. 바닥 모델은 상호작용을 보지 않는 모델이고,
  ESM-2를 상호작용 모델로 쓴 것이 아닙니다.
- 평균풀링이 불리한 이유는 짐작은 되지만(패밀리 정체성을 담는데 그것이 타깃 분리 분할에서
  전이되지 않는 바로 그 정보) **검증하지 않았습니다.**

## 5. 재현

```bash
python plm_pocket_seq.py      # 포켓 PDB 285개 -> casf_pocket_seq.csv (중앙값 152잔기)
python plm_embed.py           # ESM-2 650M 임베딩, CPU 16스레드 191초 -> casf_esm2_650M.npy (285x1280)
python plm_floor.py           # 바닥 6개 + 축소 대조 6개 -> headroom_esm2_{floor,control}.csv
```

기대값: 아미노산 개수 0.241 / ESM-2 0.179 / 리간드+개수 0.534 / 리간드+ESM-2 0.439 (±0.01).
