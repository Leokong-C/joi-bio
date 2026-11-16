"""리간드로 세면 안 되는 화학 성분 — **한 곳에서만 관리한다.**

## 왜 모듈로 뺐는가

같은 결함을 두 곳에서 따로 당했다.

- 2026-11-16, `scripts/probe_structures.py`: 당사슬(NAG)을 리간드로 세어
  SPAM1 의 홀로 구조를 1개로 보고했다. 실제로는 0개다.
- 2026-11-16, `core/structure/fetch.py`: 오렉신1 구조 14개 중 11개에서
  인지질 PGW(중원자 51개)가 실제 길항제(32개)를 제치고 `ligands[0]` 이 됐다.
  **`define_pocket` 이 그 중심으로 도킹 박스를 잡으므로, 고치지 않으면
  포켓이 아니라 막 쪽 지질 자리에 도킹한다.**

목록이 두 벌이면 한쪽만 고치게 된다. 그래서 여기 하나만 둔다.

## 기준

- **당사슬** — 단백질에 붙은 글리코실화다. 결합 리간드가 아니다
- **지질·계면활성제** — 막단백질을 녹이거나 결정을 띄우려고 넣은 것이다
- 분자량이나 원자 수로는 못 거른다. PGW 는 중원자 51개로 웬만한 약보다 크다

## 한 가지 주의

지질이 **진짜 리간드**인 경우가 있다(지질 감지 수용체, 지질 게이팅 채널).
그런 타깃을 다룰 때는 이 목록을 그대로 쓰면 안 된다 — 타깃별로 판단한다.
"""

GLYCANS = {
    "NAG", "NDG", "A2G", "BMA", "MAN", "BGC", "GLC", "GAL", "GLA", "FUC",
    "FUL", "SIA", "NGA", "XYS", "RIB", "RAM", "GCU", "IDS", "SGN", "MAL",
    "SUC", "TRE", "LAT", "GLP", "NAN",
}

LIPIDS = {
    "PGW", "3PE", "PEE", "PCW", "PC1", "POV", "PX4", "PGT", "PEF", "PSC",
    "CDL", "DGA", "LPP", "PEV", "PIO", "LHG", "PLC", "SPH", "CLR", "Y01",
    "OLC", "OLA", "PLM", "STE", "MYR", "HEX", "D10", "DD9", "TRD",
}

DETERGENTS = {
    "SOG", "BOG", "BNG", "LMT", "LMU", "DAO", "UND", "OCT", "HEZ", "C8E",
    "LDA", "12P", "7PE", "P33", "PG0", "1PG", "2PG", "F09", "JEF", "MC3",
    "SDS", "BTB", "NHE",
}

#: 리간드 후보에서 빼는 전체 집합
NON_LIGAND = GLYCANS | LIPIDS | DETERGENTS


def is_non_ligand(resname: str) -> bool:
    return resname.upper() in NON_LIGAND
