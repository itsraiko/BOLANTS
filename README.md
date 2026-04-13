# BOLANTS

**BOLANTS** (BOLtz-2 + plANTS) is a molecular docking and affinity prediction pipeline that combines [PLANTS](https://github.com/purnawanpp/plants) docking with [Boltz-2](https://github.com/jwohlwend/boltz) scoring.

This project is a fork of [Boltzina](https://github.com/ohuelab/boltzina), replacing AutoDock Vina with PLANTS as the docking engine.

![Pipeline](https://arxiv.org/html/2508.17555v1/x1.png)

## What changed from Boltzina

| | Boltzina (original) | BOLANTS |
|---|---|---|
| **Docking engine** | AutoDock Vina | PLANTS |
| **Ligand format** | PDBQT | MOL2 |
| **Receptor format** | PDBQT (Meeko) | MOL2 (SPORES/obabel) |
| **Score source** | `REMARK VINA RESULT:` | `ranking.csv → TOTAL_SCORE` |
| **Config** | Vina grid txt | PLANTS `.conf` file |

## Installation

```bash
git clone https://github.com/itsraiko/BOLANTS.git
cd BOLANTS
pip install .
```

Install PLANTS from: https://github.com/purnawanpp/plants  
Install Maxit from: https://sw-tools.rcsb.org/apps/MAXIT/

Then make sure PLANTS is executable:
```bash
chmod +x /path/to/PLANTS
```

## Usage

```bash
python run.py config.json --num_workers 4 --float32_matmul_precision medium
```

For GPUs with limited VRAM (e.g. RTX 4060), use `--float32_matmul_precision medium` to avoid OOM errors.

## Configuration File

```json
{
    "work_dir": "path/to/boltz2_output",
    "output_dir": "path/to/results",
    "receptor_pdb": "protein.pdb",
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

**PLANTS config** (`plants_input.conf`):
```
scoring_function  chemplp
search_speed      speed1
bindingsite_center  X  Y  Z
bindingsite_radius  10.0
cluster_structures  5
cluster_rmsd        2.0
```

> `protein_file`, `ligand_file`, `output_dir` are added automatically — do not include them in the conf.

## Pipeline

```
1. receptor.pdb  →  receptor.mol2       (SPORES or obabel)
2. ligand.pdb    →  ligand.mol2         (SPORES or obabel)
3. PLANTS docking → docked_ligands.mol2 + ranking.csv
4. MOL2 → PDB → CIF                    (obabel + maxit)
5. Boltz-2 affinity scoring
6. boltzina_results.csv
```

## Command Line Options

| Option | Description |
|---|---|
| `--num_workers N` | Parallel docking workers |
| `--float32_matmul_precision` | `highest` / `high` / `medium` |
| `--skip_docking` | Skip docking, only run Boltz-2 scoring |
| `--vina_override` | Re-run docking even if output exists |
| `--boltz_override` | Re-run Boltz-2 scoring even if output exists |
| `--keep_intermediate_files` | Don't clean up intermediate files |

## Before Running

You need to run Boltz-2 first to generate `manifest.json` and constraints:

```bash
boltz predict protein.yaml --out_dir boltz_results_base --use_msa_server
```

Then set `work_dir` in config.json to that output directory.

## Reference

Furui, K, & Ohue, M. Boltzina: Efficient and Accurate Virtual Screening via Docking-Guided Binding Prediction with Boltz-2. AI for Accelerated Materials Design - NeurIPS 2025. https://openreview.net/forum?id=OwtEQsd2hN
