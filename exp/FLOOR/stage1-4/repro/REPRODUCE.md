# Reproducing the CASF-2016 docking / scoring benchmark

Written for a second machine (e.g. `HP_Joo`). Nothing here depends on Linux.
Everything below was executed and timed on the primary machine:
**Windows, 16 CPU cores, 28 GiB RAM, no GPU.**

---

## 0. What does and does not transfer between machines

A Claude Science installation on another PC gets its own workspace and its own
conda environments. Opening the same project there gives you the artifact store
and the project memory; it does **not** give you the software or the
intermediate files.

| Travels with the account | Must be rebuilt on the new machine |
|---|---|
| 44 project artifacts (results, figures, reports) | 6 conda environments |
| Project memory | `bin/vina.exe` |
| Connectors (GitHub, Notion, Gmail, Drive) | `casf/` inputs and `dock/prep/` PDBQT files |

The fastest path is to restore `casf2016_dock_bundle.tar.gz` (28.1 MB,
1,424 entries) rather than regenerate inputs — see §4.

---

## 1. Environments

Two environments carry the work. Installs are additive-only, so create them
fresh rather than adding to an existing one.

```
docking : rdkit openbabel meeko spyrmsd posebusters gemmi scikit-learn
md      : openmm openff-toolkit openff-interchange mdanalysis pdbfixer
```

Exact versions used for every number in this benchmark:

| package | version |  | package | version |
|---|---|---|---|---|
| python | 3.11.17 | | rdkit | 2026.09.1 |
| numpy | 2.4.6 | | openbabel | 3.2.1 |
| scipy | 1.17.1 | | meeko | 0.8.0 |
| pandas | 3.0.6 | | spyrmsd | 0.9.0 |
| scikit-learn | 1.9.1 | | posebusters | 0.6.5 |
| gemmi | 0.7.5 | | oddt | 0.4 |

### Packages that do NOT install, and why

| package | outcome | cause |
|---|---|---|
| `vina` (Python API) | fails | wheels are manylinux/musllinux only; no conda win-64 build; no C++ compiler present (`cl`, `gcc`, `clang`, `cmake`, `swig`, `msbuild` all absent). **Use the official Windows CLI executable instead — it works.** |
| `oddt` (normal install) | refused | its doc dependency `sphinxcontrib-jsmath` ships a `.pth` autorun file. Install with `pip install oddt --no-deps`; the real runtime deps are already in `docking`. |
| `pymol` | refused | conda's `qt-main` ships a link-phase hook script; the transaction is rejected. |
| `openmmforcefields` | unavailable | no win-64 build. Bypassed via OpenFF Interchange; AM1-BCC charges replaced by NAGL (`openff-gnn-am1bcc-1.0.0.pt`, 1.9 s for a 50-atom ligand). |

---

## 2. Binary

AutoDock **Vina 1.2.7**, official Windows CLI release:

```
bin/vina.exe
size    1,233,920 bytes
sha256  e0c4b2715e0c1a74...
```

Verify the size and hash prefix before trusting any numbers; a different build
will not reproduce them.

---

## 3. Data

* Target list: the CASF-2016 core set, 285 complexes = 57 target clusters x 5
  ligands, from `CoreSet.dat`. Experimental logKa 2.07-11.82 (mean 6.49,
  SD 2.17).
* Receptors: `https://files.rcsb.org/download/{code}.pdb`, ATOM records only.
* Ligands: RCSB ModelServer, which preserves bond orders —
  `https://models.rcsb.org/v1/{code}/ligand?auth_asym_id=..&auth_seq_id=..&encoding=sdf`

Ligand identity is resolved from the mmCIF: drop water, ions and crystallisation
additives, drop any `_chem_comp.type` that is peptide/saccharide/nucleic-linking,
then take the remaining candidate with the most heavy atoms. This resolves 284
of 285 automatically and was validated against RCSB affinity annotations on a
50-complex sample (49/50).

### Two corrected complexes

| code | problem | fix |
|---|---|---|
| `3syr` | the automatic pick chose the cofactor `LLP` (24 heavy atoms); the real ligand is `GPK` (20) | `casf/lig_fix/3syr.sdf`. Scoring the wrong ligand gave RMSD 4.41 A — a **false docking failure**; with GPK the top-1 RMSD is **0.27 A**. |
| `1a30` | the ligand is a tripeptide (Glu-Asp-Leu, chain C res 506-508), excluded by the peptide-linking filter | `casf/lig_fix/1a30_pep.pdb`, 26 heavy atoms, SMILES `CC(C)C[C@H](NC(=O)[C@H](CC(=O)O)NC(=O)[C@@H](N)CCC(=O)O)C(=O)O` |

### One counting pitfall

Heavy-atom counts taken from mmCIF HETATM groups double-count **alternate
locations**. Five complexes are affected — `3gv9`, `2xnb`, `2yfe`, `2al5`
(exactly 2x) and `3prs` (50 vs 53). Example: ligand `GV9` at A/364 in `3gv9` is
stored as altloc A and altloc B, 12 atoms each, and counts as 24.
**Use the heavy-atom count from the SDF.** Getting this wrong shifts the
ligand-size null model from its true R = +0.484 down to +0.469.

---

## 4. Pipeline

```
receptor : PDB ATOM records  ->  obabel -xr -p 7.4  ->  .pdbqt
ligand   : ModelServer SDF   ->  RDKit AddHs        ->  Meeko  ->  .pdbqt
dock     : vina.exe, box centre (cx,cy,cz) and side = ligand extent + 10 A
RMSD     : sPyRMSD symmrmsd on heavy atoms (symmetry-corrected)
validity : PoseBusters config="redock"
```

Baseline settings: `--exhaustiveness 8 --num_modes 9 --seed 42`, 16 cores,
**44 min** for 284 complexes, 0 failures.

Results are appended to JSONL and flushed after every complex, so any run can be
killed and resumed. Do this — the workspace is swept after idle periods, and one
completed 120-run experiment was lost before this was in place.

---

## 5. Reproducibility caveat you will hit

Rerunning the pipeline reproduces every headline number exactly (top-1 <= 2 A:
66.2 % both times, zero success-verdict flips) but only **250 of 284** scores
match to three decimals.

Cause, pinned precisely: **Open Babel places polar hydrogens
non-deterministically across processes.** 26 of 4,190 lines in one receptor
PDBQT differ between two preparations (added H on Ser/Cys). Rescoring with the
archived receptor reproduced -8.479 kcal/mol; a freshly prepared receptor gave
-8.355. **Vina itself is deterministic.**

> Bit-level reproduction therefore requires archiving the prepared PDBQT files,
> not just the pipeline. That is what `casf2016_dock_bundle.tar.gz` is for.

---

## 6. Scripts in this directory

| file | what it does | runtime |
|---|---|---|
| `null_models.py` | 14 trivial descriptors -> ridge and random forest, GroupKFold(10) by target | < 1 min |
| `score_mlsf.py` | RF-Score v1/v2 via ODDT, trained on PDBbind v2016 refined minus core | ~25 s training + ~7 min for 285 |

Both encode non-obvious workarounds in their module docstrings. Read those
before editing.

---

## 7. Measured runtimes (16 cores)

| experiment | n | wall time |
|---|---|---|
| baseline dock, exhaustiveness 8 | 284 | 44 min |
| exhaustiveness 32 | 284 | 215.6 min (34.8 CPU-hours) |
| seed repeats (43, 44) | 2 x 284 | ~88 min |
| Vinardo rescoring | 284 | ~45 min |
| box-size sensitivity, 3 paddings | 60 x 3 | 24 min |
| RF-Score v1+v2 | 285 | ~7 min |
| MD, implicit solvent, smallest complex (2,121 atoms) | 1 | **1.63 ns/day** |

The MD number is the reason the MD axis is restricted to minimisation-based
MM/GBSA: 200 minimisation iterations take 72 s, and trajectories long enough to
mean anything are out of reach on CPU.
