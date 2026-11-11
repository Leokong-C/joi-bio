# Stage 4 — extended benchmark: null models, machine learning, and a language model

All numbers here are for the **CASF-2016 core set, n = 285** (57 targets x 5 ligands).
Stage 1-3 reported n = 284; stage 4 adds `1a30` (peptide ligand) and corrects `3syr`
(see *Corrections* below), so stage-4 numbers supersede the earlier ones wherever they differ.

## The question

Stage 1-3 established that AutoDock Vina's weakness on CASF-2016 is **ranking, not sampling**:
the correct pose is generated for 85.9 % of complexes but ranked first for only 66.2 %, and
quadrupling the search effort (34.8 CPU-hours) moves that by +0.3 pt — less than three times the
seed-to-seed noise. Stage 4 asks the next question: **is a scoring function's correlation with
affinity evidence that it models binding at all?**

The test is adversarial. We build a null model from 14 descriptors that contain **no model of
binding** — ligand bulk and physicochemistry plus raw counts of how many protein heavy atoms sit
near the ligand — and ask every scoring function to beat it.

## Result 1 — physics-based scoring is almost entirely subsumed by the null

| predictor | Pearson R | residual SD (log) | within-target rho | rho >= 0.9 |
|---|---|---|---|---|
| RF-Score v2 (ML) | **0.741** | **1.46** | 0.561 | 12.3 % |
| RF-Score v1 (ML) | 0.713 | 1.53 | 0.514 | 8.8 % |
| trivial null, random forest | 0.547 | 1.82 | 0.512 | 8.8 % |
| trivial null, ridge | 0.526 | 1.85 | 0.535 | 8.8 % |
| Vina | 0.502 | 1.88 | 0.437 | 15.8 % |
| **ligand heavy atoms alone** | 0.482 | 1.91 | **0.582** | **17.9 %** |
| Vinardo | 0.428 | 1.97 | 0.411 | 7.0 % |

Cross-validation for both null models is `GroupKFold(10)` **grouped by CASF target cluster**, so the
four other ligands of a target are never in the training fold with it. Without that grouping the
null's scores are inflated by target-level leakage.

The decisive statistic is what each function adds **on top of** the null
(`results/casf2016_null_decomposition_corrected.csv`):

| predictor | partial R after removing the null | unique R^2 added |
|---|---|---|
| Vina | 0.242 | **0.041** |
| Vinardo | 0.195 | **0.027** |
| RF-Score v2 | 0.598 | **0.251** |

Vina contributes 4 % of the affinity variance that the trivial descriptors do not already explain.
RF-Score v2 contributes 25 %. **The machine-learned function carries real information; the two
physics-based functions are close to redundant with counting atoms.**

The single most important feature of the null is the number of **protein heavy atoms within 6 A of
the ligand** (32.0 % of random-forest importance) — pocket density, not chemistry.

## Result 2 — ranking is a different problem, and nothing has solved it

On within-target ordering — putting the five ligands of one protein in the right order, which is
what matters in a lead-optimisation campaign — **counting the ligand's heavy atoms is still the best
predictor on this set** (mean rho 0.582, and rho >= 0.9 for 17.9 % of targets, versus 12.3 % for
RF-Score v2 and 15.8 % for Vina). Improving the global affinity correlation does not carry over to
ranking.

## Result 3 — a language model, zero-shot, with an abstention option

285 complexes x 3 information conditions = 855 structured predictions (`claude-sonnet-5`, thinking
disabled, forced tool call returning `abstain` / `pKd` / `confidence` / `recognized_complex` /
`basis`). 0 failures, 2.3 min, 652k input and 137k output tokens.

| condition | R | within-target rho | MAE | abstains | claims recognition |
|---|---|---|---|---|---|
| A: ligand SMILES only | 0.384 | 0.366 | 1.96 | **96.8 %** | 4.9 % |
| B: + protein name | 0.655 | 0.521 | 1.38 | 20.4 % | 26.0 % |
| C: + PDB ID | **0.689** | **0.554** | **1.30** | 8.8 % | 44.2 % |

The language model beats Vina and both null models on affinity correlation and sits below
RF-Score v2. Denied the target identity (condition A) it abstains 97 % of the time — it recognises
that the question is not answerable from a SMILES string alone.

### Contamination check: negative

CASF-2016 affinities are public, so the obvious objection is memorisation. Three independent
checks say no:

1. Complexes the model **claims** to recognise are predicted **worse**, not better — in condition C,
   R = 0.641 (n = 126) for claimed recognition versus R = 0.730 (n = 159) without.
   Memorisation would produce the opposite ordering.
2. Supplying the PDB ID flipped `recognized_complex` from false to true for **56 of 285** complexes,
   yet their mean absolute error improved by only **-0.11** log units — indistinguishable from the
   **-0.08** seen across the other 229.
3. Predictions landing within 0.1 log of the experimental value: **4.2 %** without the PDB ID,
   **3.9 %** with it. The identifier does not unlock stored values.

What the PDB ID actually buys is **confidence, not accuracy**: abstention falls from 20.4 % to
8.8 % while MAE moves 1.38 -> 1.30. That is the case for building abstention into the protocol.

## Result 4 - MM/GBSA: the costliest physics method adds nothing

Vina is an empirical scoring function, so the obvious objection to Result 1 is that the problem is
Vina's approximations rather than physics-based scoring as such. We therefore rescored the crystal
poses with minimisation-based, single-trajectory MM/GBSA - a far more expensive calculation with a
real force field and an implicit-solvent model.

Method: OpenMM 8.6.1 on CPU, receptor repaired with PDBFixer, `ff14sb_off_impropers_0.0.4` for
protein and `openff-2.2.0` for ligand, NAGL (`openff-gnn-am1bcc-1.0.0.pt`) partial charges,
`GBSAOBCForce` (OBC2, solvent dielectric 78.5, SA energy 2.25936 kJ/mol/nm^2), 200 minimisation
iterations on the complex, then component energies from that one minimised geometry. Run on a
**random 47 of the 285** (seed 42, 87 % success, 252 s per complex); the sample is random rather
than stratified so the subset is unbiased, and all comparators below are recomputed on the same 47.

| predictor | Pearson R | cost per complex |
|---|---|---|
| RF-Score v2 (ML) | **0.756** | 1.3 s |
| trivial null, random forest | 0.549 | < 0.01 s |
| ligand heavy atoms / trivial null, ridge | 0.528 | < 0.01 s |
| **MM/GBSA** | **0.385** | **252 s** |
| Vina | 0.368 | 94 s |
| Vinardo | 0.329 | 94 s |

MM/GBSA edges past Vina and loses heavily to the trivial null. Its unique contribution beyond the
null is **dR^2 = 0.001 (p = 0.80)** - statistically nothing; Vina gives 0.003 (p = 0.69) and
Vinardo 0.001 (p = 0.80) on this subset, while RF-Score v2 adds **0.271 (p = 3e-6)** at 1/200 of
the cost.

**This strengthens Result 1 rather than qualifying it.** The failure of physics-based scoring on
CASF-2016 is not an artefact of Vina's empirical functional form: computing the energy properly,
with a force field and implicit solvent, lands in the same place.

Four implementation traps, all documented in `repro/mmgbsa_one.py`:

1. `Topology.from_pdb` returns one molecule **per chain**; taking `molecules[0]` silently truncates
   the receptor and crashes `GBSAOBCForce` with "Index out of range" (11 of the first 17 complexes).
2. PDBFixer's `removeHeterogens` keeps free standard residues that sit in ATOM records as their own
   chain - a zwitterionic proline in `3bv9` - and `create_interchange` then raises
   `NonIntegralMoleculeChargeError`. Drop molecules under 60 atoms.
3. Several deposited ligand SDFs have unspecified stereocentres, so `Molecule.from_file` needs
   `allow_undefined_stereo=True`. MM parameters do not depend on the assignment.
4. OpenMM kills the whole interpreter with `STATUS_ACCESS_VIOLATION` on some complexes, so each
   complex runs in its own subprocess with a resumable JSONL log.

Seven of 54 attempts failed and are listed with their errors in `results/casf2016_mmgbsa_failures.csv`
(two `UnassignedChemistryInPDBError` on nonstandard residues, two pre-fix stereochemistry errors,
one access violation, one timeout, one charge error).

## Corrections to earlier stages

* **n = 284 -> 285.** `1a30` (tripeptide ligand Glu-Asp-Leu, excluded earlier by the
  peptide-linking filter) is now included, and `3syr`'s ligand is corrected from the cofactor `LLP`
  to the real ligand `GPK`. With the correct ligand `3syr`'s top-1 RMSD is **0.27 A**, not 4.41 A —
  the earlier entry was a **false docking failure**. Including `1a30`, which is an outlier, lowers
  the null models slightly (ridge 0.553 -> 0.526, forest 0.556 -> 0.547).
* **Ligand-size null R = +0.484, not +0.469.** Heavy-atom counts taken from mmCIF HETATM groups
  double-count **alternate locations**. Five complexes are affected: `3gv9`, `2xnb`, `2yfe`, `2al5`
  (exactly 2x) and `3prs` (50 vs 53). Verified directly: ligand `GV9` at A/364 in `3gv9` is stored
  as altloc A and altloc B, 12 atoms each, and was counted as 24. Use the SDF-derived count.
  Vina's size-controlled partial R moves from +0.355 to **+0.343** (+0.339 was the n = 284
  value and appears in the first version of this file, the Notion page and the summary email).

### Figure provenance

Both figures that carried pre-correction numbers have been rebuilt at n = 285:
`figures/casf2016_poster.png` (v2, eight stages) and `figures/casf2016_size_decomposition.png`
(partial R = +0.343, variance split 0.071 / 0.162 / 0.090, total R^2 = 0.323). No figure in this
directory now shows a superseded value. `repro/poster_v3.py` is the poster script and
`results/poster_numbers.json` holds the 49 values it drew.

## Does ML scoring need a Linux machine?

No. ODDT 0.4 runs on Windows and trained RF-Score v1/v2/v3 from PDBbind v2016 in 22-31 s each.
Three workarounds are required, none of them OS-specific, and all three are documented in
`repro/score_mlsf.py`:

1. install with `pip install oddt --no-deps` — its documentation dependency `sphinxcontrib-jsmath`
   ships a `.pth` autorun file that the sandbox refuses;
2. shim `np.in1d = np.isin` before importing oddt — removed in NumPy 2.0;
3. use the RDKit backend — the Open Babel backend imports the legacy top-level `pybel` module.

RF-Score v3 is excluded: its AutoDock Vina term calls `np.hstack` on a generator, which NumPy 2
rejects. v2 scores highest of the three anyway.

**No leakage.** ODDT trains on `train_set='refined' & ~core`, i.e. the PDBbind refined set with the
core set removed, and tests on the core set. Verified by reading `oddt/scoring/functions/RFScore.py`.

Only GPU deep-learning co-folding (DiffDock, Boltz, AlphaFold3) genuinely needs different hardware.

## Layout

```
results/    22 tables - the numbers behind every claim above
figures/     8 panels - see STAGE4.md section they belong to
repro/       REPRODUCE.md + 2 scripts, for running this on another machine
```

`repro/REPRODUCE.md` lists exact package versions, the Vina binary checksum, what does and does not
transfer between machines, measured runtimes, and the one reproducibility caveat that matters:
**Open Babel places polar hydrogens non-deterministically across processes**, so bit-level
reproduction requires archiving the prepared PDBQT files, not just the pipeline. Vina itself is
deterministic.
