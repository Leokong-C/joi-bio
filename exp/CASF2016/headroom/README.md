# HEADROOM

**How much room is actually left on a protein–ligand scoring benchmark?**

A benchmark score means nothing on its own. It has to be read between two bounds:

- a **floor** — what a model reaches without ever looking at the protein–ligand
  interaction (ligand descriptors, pocket composition, or just target identity);
- a **ceiling** — what the labels permit, given that binding constants measured
  in different laboratories disagree.

HEADROOM measures both on the same complexes with the same labels, and places
every predictor between them. On CASF-2016 the band turns out to be narrow, and
its width depends less on the model than on how the data were split.

## What it reproduces

| Quantity | Value on CASF-2016 (n = 285) |
| --- | --- |
| Ceiling, replicate-only estimate | **0.882** |
| Ceiling, independent literature value (Kramer et al. 2012) | 0.90 |
| Ceiling, cross-source estimate (optimistic — see below) | 0.964 |
| Floor, BothBias, target-disjoint split | **0.534** |
| Floor, BothBias, random split | **0.685** |
| Vina / Vinardo | 0.502 / 0.428 |
| RF-Score v1 / v2 | 0.713 / 0.741 |

Two results drive the conclusions:

1. **Vina (0.502) and Vinardo (0.428) sit below BothBias (0.534)** — below a
   model that never sees the two molecules together.
2. **The floor moves with the split, not the model.** On identical data with
   identical features, ProteinBias goes from 0.241 (target unseen) to 0.509
   (target seen). Pocket composition is acting as a lookup key for target
   identity, not as a physical description of binding.

## Setup

```bash
conda env create -f environment.yml
conda activate headroom
```

Network access is required to `data.rcsb.org` and `www.bindingdb.org`. On a
sandboxed host both may need to be allowlisted first.

## Inputs you must supply

| File | Columns | Where it comes from |
| --- | --- | --- |
| `casf_index.csv` | `code`, `logKa`, `target` | CASF-2016 core set (`CoreSet.dat`): PDB id, log Ka, target cluster |
| `casf_smiles.csv` | `code`, `smiles` | ligand SMILES per complex (RCSB ModelServer SDF → RDKit) |
| `pockets/` | `<code>.pdb` | receptor trimmed to residues within 16 Å of the ligand |

## Pipeline

```bash
# 1. map the core set to UniProt and pull BindingDB affinities  (~2 min, 75 MB)
python 01_pull.py --index casf_index.csv --work work

# 2. label-noise ceiling                                        (~2 min)
python 02_ceiling.py --index casf_index.csv --smiles casf_smiles.csv --work work

# 3. the floor, both splits                                     (~1 min)
python 03_floor.py --index casf_index.csv --smiles casf_smiles.csv \
    --ligand-features casf_features.csv --pockets pockets/ --work work

# descriptors for any SMILES table (step 3 and 4 inputs)
python featurize.py --in casf_smiles.csv --out casf_features.csv
python featurize.py --in work/bindingdb_train.csv --out work/bindingdb_features.csv

# 4. the external-training floor, from BindingDB                 (~20 min)
python 04_external_floor.py --index casf_index.csv --smiles casf_smiles.csv \
    --casf-uniprot work/casf_uniprot.csv --matched work/matched_records.csv \
    --train-features work/bindingdb_features.csv --test-features casf_features.csv \
    --work work
```

## How the ceiling is estimated

Two independent measurements of the same quantity differ by an amount that
carries both measurements' error, so

```
sigma_error = SD(differences) / sqrt(2)
ceiling_R   = sqrt(1 - sigma_error^2 / SD(label)^2)
```

This is the estimator of Kalliokoski et al. 2013 (doi:10.1371/journal.pone.0061007).

**The cross-source number is optimistic.** 20.1 % of matched BindingDB Ki/Kd
records sit within 0.02 log units of the PDBbind label: they re-report the same
primary measurement rather than an independent one, so the two sources agree
more than two laboratories would. Dropping those records is *not* a fix — it
keeps only the disagreeing ones and inflates the spread. The replicate-only
estimate compares BindingDB measurements to each other and never touches the
label, so report that one.

## Constraints found the hard way

- **PDBbind General is not scriptable.** `pdbbind-cn.org` requires registration;
  Zenodo and Hugging Face mirrors may be unreachable on a sandboxed host. RCSB
  is not a substitute — in a random sample of 24 PDBbind codes it carried a
  binding affinity for 9 and a Ki/Kd for only 2. Step 4 uses BindingDB instead,
  which needs no structures for a ligand-only model.
- **Process pools may be unavailable.** A locked-down container can refuse the
  named pipes `ProcessPoolExecutor` needs (Windows error 5). `featurize.py`
  detects this and falls back to serial; RDKit already threads internally, so
  ~39,000 molecules take about 13 minutes.
- **BindingDB's PDB-keyed endpoint returns HTTP 500.** Use the UniProt-keyed one
  and map PDB → UniProt through RCSB, as step 1 does.
- **The response key is misspelled upstream** (`getLindsByUniprotsResponse`).
  Read the single top-level value rather than keying on the name.

## Prior art this builds on

- Durant et al., *Robustly interrogating machine learning-based scoring
  functions: what are they learning?*, Bioinformatics 41(2) btaf040 — bias-only
  baselines reach state-of-the-art on CASF-2016 when trained on PDBbind General.
- Boyles, Deane & Morris, *Learning from the ligand*, Bioinformatics 36, 758 —
  a ligand-only model predicts about as well as a classical scoring function.
- Kramer, Kalliokoski, Gedeck & Vulpetti, *The experimental uncertainty of
  heterogeneous public Ki data*, J Med Chem 2012, doi:10.1021/jm300131x —
  maximum attainable R² on large public sets ≈ 0.81.

## One number to be careful with

Our top-1 RMSD ≤ 2 Å of 66.3 % comes from **re-docking from scratch**. It is not
CASF-2016's docking-power figure, which scores a *provided* set of decoy poses.
The two are not comparable; label ours as re-docking wherever it appears.
