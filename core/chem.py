#!/usr/bin/env python3
"""
core.chem — 화합물 조회와 유사도. 구조 선택과 라이브러리 구성이 함께 쓴다.

구조 선택이 이걸 쓰는 이유: 도킹에서는 해상도보다 '어떤 리간드가 결합한
구조인가'가 중요하다. 벤치마크 화합물과 닮은 리간드가 들어있는 구조를 쓰면
포켓이 그 화합물을 수용하는 형태로 이미 열려 있다. 리지드 도킹은 측쇄
재배열을 모델링하지 못하므로, 이 선택이 결과를 좌우한다.
"""
from __future__ import annotations

import time

import requests
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
RCSB_CHEMCOMP = "https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id}"
TIMEOUT = 90

_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def _get(url: str, params: dict | None = None, retries: int = 3):
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params or {}, timeout=TIMEOUT)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 404:
                return None
        except requests.RequestException:
            pass
        time.sleep(2 * (attempt + 1))
    return None


def chembl_by_name(name: str) -> dict | None:
    """화합물명으로 ChEMBL에서 SMILES를 찾는다. 정확 일치 → 동의어 순."""
    for params in (
        {"pref_name__iexact": name, "limit": 5},
        {"molecule_synonyms__molecule_synonym__icontains": name, "limit": 5},
    ):
        data = _get(f"{CHEMBL}/molecule.json", params)
        for m in (data or {}).get("molecules", []):
            smi = (m.get("molecule_structures") or {}).get("canonical_smiles")
            if smi:
                return {
                    "smiles": smi,
                    "chembl_id": m["molecule_chembl_id"],
                    "pref_name": m.get("pref_name") or name,
                    "query_name": name,
                }
    return None


def rcsb_ligand_smiles(comp_id: str) -> str | None:
    """PDB 화학성분 ID(예: 020)의 SMILES를 RCSB에서 가져온다."""
    data = _get(RCSB_CHEMCOMP.format(comp_id=comp_id))
    if not data:
        return None
    # 우선순위: OpenEye/CACTVS canonical SMILES
    best = None
    for d in data.get("pdbx_chem_comp_descriptor", []) or []:
        if (d.get("type") or "").upper() != "SMILES_CANONICAL":
            continue
        val = d.get("descriptor")
        if not val:
            continue
        if Chem.MolFromSmiles(val) is not None:
            return val
        best = best or val
    for d in data.get("pdbx_chem_comp_descriptor", []) or []:
        if (d.get("type") or "").upper() == "SMILES":
            val = d.get("descriptor")
            if val and Chem.MolFromSmiles(val) is not None:
                return val
    return best


def tanimoto(smiles_a: str, smiles_b: str) -> float | None:
    """Morgan(r=2, 2048bit) Tanimoto 유사도."""
    ma, mb = Chem.MolFromSmiles(smiles_a), Chem.MolFromSmiles(smiles_b)
    if ma is None or mb is None:
        return None
    fa, fb = _MORGAN.GetFingerprint(ma), _MORGAN.GetFingerprint(mb)
    return round(DataStructs.TanimotoSimilarity(fa, fb), 4)


def max_similarity(smiles: str, references: list[str]) -> tuple[float, int]:
    """references 중 가장 닮은 것과의 유사도, 그리고 그 인덱스."""
    best, idx = 0.0, -1
    for i, ref in enumerate(references):
        s = tanimoto(smiles, ref)
        if s is not None and s > best:
            best, idx = s, i
    return best, idx


def heavy_atoms(smiles: str) -> int | None:
    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumHeavyAtoms() if mol else None


def rotatable_bonds(smiles: str) -> int | None:
    from rdkit.Chem import rdMolDescriptors

    mol = Chem.MolFromSmiles(smiles)
    return rdMolDescriptors.CalcNumRotatableBonds(mol) if mol else None
