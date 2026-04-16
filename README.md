# BOLANTS

**BOLANTS** (BOLtz-2 + plANTS) is a virtual screening pipeline that combines
[PLANTS](https://github.com/purnawanpp/plants) molecular docking with
[Boltz-2](https://github.com/jwohlwend/boltz) AI-based binding affinity prediction.

Given a receptor structure, a set of ligand candidates, and a binding site location,
BOLANTS automatically docks each ligand into the binding pocket and predicts how
strongly it binds — producing a ranked CSV you can use to prioritize compounds for
experimental follow-up.

This project is a fork of [Boltzina](https://github.com/ohuelab/boltzina),
replacing AutoDock Vina with PLANTS as the docking engine.

---

## How It Works

```mermaid
flowchart TD
    subgraph IN["📥  Inputs"]
        direction LR
        A["📄 protein.yaml\n(sequence + ligand SMILES)"]
        B["🧬 protein.pdb\n(receptor structure)"]
        C["🔬 Ligand SMILES\n(one per compound)"]
    end

    subgraph PREP["⚙️  Preparation"]
        direction TB
        A -->|"boltz predict\n--use_msa_server"| D["📋 manifest.json\n+ constraints.npz"]
        B -->|"SPORES / obabel"| E["🗂️ receptor.mol2"]
        C -->|"ligand_preparation.py"| F["📁 ligand.pdb\n(3D conformer, unique atom names)"]
        F -->|"obabel"| G["🗂️ ligand.mol2"]
    end

    subgraph DOCK["🌿  PLANTS Docking"]
        direction TB
        E --> PL(["⚡ PLANTS\nchemplp scoring"])
        G --> PL
        PL --> I["📦 docked_ligands.mol2\nranking.csv  →  TOTAL_SCORE"]
    end

    subgraph SCORE["🤖  Boltz-2 Affinity Scoring"]
        direction TB
        I -->|"obabel + maxit"| J["📄 complex.cif\n(protein–ligand complex)"]
        D --> BZ(["🧠 Boltz-2\naffinity model"])
        J --> BZ
    end

    BZ --> OUT["📊 boltzina_results.csv\naffinity_pred_value · docking_score"]

    IN --> PREP
    PREP --> DOCK
    DOCK --> SCORE
```

### Stage 0 — Input preparation (`prepare.py`)

`prepare.py` takes your raw inputs and generates everything the pipeline needs:

| Input | What it produces |
|---|---|
| `protein.pdb` | `protein_clean.pdb` (single chain, ATOM records only) |
| protein sequence | `protein.yaml` (Boltz-2 format: sequence + reference SMILES) |
| binding site X Y Z | `plants_input.conf` (binding pocket definition) |
| ligand SMILES / SDF / MOL2 / PDB | `ligands/input_pdbs/*.pdb` (3D conformers via RDKit) |
| — | `ligands/prepared_mols.pkl` (RDKit mol objects for later use) |
| all of the above | `config.json` (master config pointing to all files) |

Ligand 3D conformers are generated using RDKit + Boltz-2's `compute_3d_conformer`,
which produces geometries compatible with the rest of the pipeline.

---

### Stage 1 — Boltz-2 structure prediction (`boltz predict`)

Before any docking, BOLANTS automatically runs:

```
boltz predict protein.yaml --use_msa_server
```

This step runs **once per protein** and produces:

- **`manifest.json`** — token-level description of the system: which tokens belong
  to the protein, which to the ligand, and affinity metadata (chain ID, molecular weight).
  Required by the Boltz-2 affinity model.
- **`constraints.npz`** — ligand chemical graph (bond orders, atom types, chirality).
  Boltz-2 uses this to understand the ligand's covalent structure.
- **MSA files** — multiple sequence alignment pulled from an external server.
  Encodes the protein's evolutionary history, which Boltz-2 uses to infer structural context.

The output is cached — subsequent runs with the same protein reuse it without re-running prediction.

---

### Stage 2 — Receptor conversion

```
protein.pdb  →  receptor.mol2
```

PLANTS works exclusively with MOL2 format. BOLANTS converts the receptor using
**SPORES** (preferred, better atom-type assignment) or **obabel** as a fallback.
Accurate atom types in the MOL2 file directly affect docking quality.

---

### Stage 3 — PLANTS docking (per ligand)

For each ligand:

```
ligand.pdb  →  ligand.mol2  →  PLANTS  →  docked_ligands.mol2 + ranking.csv
```

PLANTS searches the binding pocket defined by the center coordinates and radius in
`plants_input.conf`. It uses the **ChemPLP** scoring function, which accounts for:
- Steric clashes and shape complementarity
- Hydrogen bond geometry
- Hydrophobic contacts
- Metal coordination (if applicable)

With `cluster_structures 5`, PLANTS generates the 5 best distinct poses
(clustered by RMSD ≥ 2.0 Å) and reports their scores in `ranking.csv`.
The top-ranked pose (`TOTAL_SCORE`) is selected for Boltz-2 scoring.

Multiple ligands can be docked in parallel using `--num_workers N`.

---

### Stage 4 — Complex assembly

The docked ligand pose is merged with the receptor into a single structure:

```
receptor.pdb + docked_ligand.mol2  →  complex.pdb  →  complex.cif
```

The PDB→CIF conversion is done by **maxit**, which produces a standards-compliant
mmCIF file that Boltz-2 can read. This step is necessary because Boltz-2's affinity
model reads mmCIF input, not PDB.

---

### Stage 5 — Boltz-2 affinity scoring (per ligand)

```
complex.cif + manifest.json + constraints.npz  →  Boltz-2 affinity model  →  affinity_pred_value
```

Boltz-2 does **not** predict a new structure here. It takes the docked complex as
input and predicts the binding affinity directly — this is much faster than a full
structure prediction run.

The model outputs a continuous value (approximately in pKd/pKi units):
a higher value means stronger predicted binding.

---

### Stage 6 — Results

All per-ligand scores are collected into a single CSV:

```
boltzina_results.csv
```

| Column | Source | Interpretation |
|---|---|---|
| `ligand` | ligand file name | compound identifier |
| `affinity_pred_value` | Boltz-2 | predicted binding strength (higher = better) |
| `docking_score` | PLANTS `TOTAL_SCORE` | physical fit in the pocket (more negative = better) |

Compounds with **both** a good docking score and a high affinity prediction are the
strongest hits. Using both scores together reduces false positives compared to either
metric alone.

---

## What changed from Boltzina

| | Boltzina (original) | BOLANTS |
|---|---|---|
| **Docking engine** | AutoDock Vina | PLANTS |
| **Ligand format** | PDBQT | MOL2 |
| **Receptor format** | PDBQT (Meeko) | MOL2 (SPORES/obabel) |
| **Score source** | `REMARK VINA RESULT:` | `ranking.csv → TOTAL_SCORE` |
| **Config** | Vina grid txt | PLANTS `.conf` file |
| **Input preparation** | manual | `prepare.py` (fully automated) |

---

## Installation

```bash
git clone https://github.com/itsraiko/BOLANTS.git
cd BOLANTS
bash setup.sh
```

`setup.sh` does the following automatically:
1. `pip install .` — installs all Python dependencies (boltz, openbabel, pdb-tools, …)
2. Downloads Boltz-2 model weights to `~/.boltz/`
3. Downloads and builds **maxit** from source

**PLANTS** must be installed separately (platform-specific binary):
```bash
# Download from: https://github.com/purnawanpp/plants
chmod +x /path/to/PLANTS
```
Then set `plants_bin` in your `config.json`, or add PLANTS to your `PATH`.

> **SPORES** (optional): if installed and in PATH, BOLANTS uses it instead of obabel
> for better MOL2 atom-type assignment. Download from the same PLANTS repository.

---

## Quick Start

### Step 1 — Prepare input files

**Option A — SMILES file:**

```bash
python prepare.py \
    --protein  protein.pdb \
    --smiles   ligands.smi \
    --center   X Y Z \
    --output   my_run \
    --plants_bin /path/to/PLANTS
```

SMILES file format (`ligands.smi`):
```
O=C(N...)c1n...   ZINC000343638897
CC(=O)Nc1...      ZINC000012345678
```

**Option B — Structure files (SDF / MOL2 / MOL / PDB):**

```bash
python prepare.py \
    --protein  protein.pdb \
    --ligands  ligand1.sdf ligand2.mol2 ligand3.pdb \
    --center   X Y Z \
    --output   my_run \
    --plants_bin /path/to/PLANTS
```

You can mix formats freely; multi-molecule SDF files are split automatically.
Both options generate `my_run/config.json` and all required input files.

### Step 2 — Run

```bash
python run.py my_run/config.json --float32_matmul_precision medium
```

> For GPUs with limited VRAM (e.g. RTX 4060), `--float32_matmul_precision medium`
> significantly reduces memory usage with minimal accuracy impact.

Results are saved to `my_run/results/boltzina_results.csv`.

---

## `prepare.py` Options

| Option | Default | Description |
|---|---|---|
| `--protein` | *(required)* | Receptor PDB file |
| `--smiles` | *(required\*)* | SMILES file (`SMILES NAME` per line) |
| `--ligands` | *(required\*)* | One or more structure files (SDF/MOL2/MOL/PDB) |
| `--center X Y Z` | *(required)* | Binding site center coordinates (Å) |
| `--radius` | `10.0` | Binding site radius (Å) |
| `--output` | `bolants_run` | Output directory |
| `--fname` | PDB filename stem | Project name used in output filenames |
| `--chain` | first chain | Protein chain to use |
| `--plants_bin` | — | Path to PLANTS binary |
| `--ligand_name` | `UNL` | Residue name for ligands in PDB files |

\* `--smiles` and `--ligands` are mutually exclusive; one is required.

---

## Configuration File (`config.json`)

```json
{
    "receptor_pdb": "protein.pdb",
    "boltz_yaml": "protein.yaml",
    "output_dir": "path/to/results",
    "plants_config": "plants_input.conf",
    "plants_bin": "/path/to/PLANTS",
    "fname": "my_protein",
    "input_ligand_name": "UNL",
    "ligand_files": [
        "ligand1.pdb",
        "ligand2.pdb"
    ]
}
```

**`boltz_yaml`** (recommended): if set, BOLANTS automatically runs `boltz predict`
before docking if the output does not yet exist. No manual pre-run needed.

**`work_dir`** (optional): path to an existing `boltz_results_<fname>` directory if
you want to reuse a previous `boltz predict` run. Derived automatically from
`boltz_yaml` if omitted.

**`use_msa_server`** (optional, default `true`): set to `false` to run
`boltz predict` without `--use_msa_server` (useful in offline environments).

**`prepared_mols_file`** (optional): path to the `.pkl` file produced by `prepare.py`.
Used to look up molecular weights for the Boltz-2 affinity manifest.

**PLANTS config** (`plants_input.conf`):
```
scoring_function  chemplp
search_speed      speed1
bindingsite_center  X  Y  Z
bindingsite_radius  10.0
cluster_structures  5
cluster_rmsd        2.0
```

> `protein_file`, `ligand_file`, and `output_dir` are injected automatically —
> do not include them in the conf file.

---

## `run.py` Command Line Options

| Option | Description |
|---|---|
| `--num_workers N` | Number of parallel docking workers (default: 1) |
| `--float32_matmul_precision` | `highest` / `high` / `medium` — trade speed for VRAM |
| `--skip_docking` | Skip docking, run Boltz-2 scoring only |
| `--vina_override` | Re-run docking even if output already exists |
| `--boltz_override` | Re-run Boltz-2 scoring even if output already exists |
| `--keep_intermediate_files` | Keep intermediate files (mol2, cif, etc.) after run |

---

## Reference

Furui, K, & Ohue, M. Boltzina: Efficient and Accurate Virtual Screening via
Docking-Guided Binding Prediction with Boltz-2. AI for Accelerated Materials Design
— NeurIPS 2025. https://openreview.net/forum?id=OwtEQsd2hN
