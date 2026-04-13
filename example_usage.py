#!/usr/bin/env python3
"""
Example usage of BOLANTS (BOLtz-2 + plANTS).

Minimal config.json example:
{
    "receptor_pdb": "protein.pdb",
    "boltz_yaml": "protein.yaml",       # auto-runs boltz predict if needed
    "plants_config": "plants.conf",
    "plants_bin": "/path/to/PLANTS",
    "output_dir": "results",
    "fname": "myprotein",
    "input_ligand_name": "UNL",
    "ligand_files": ["ligands/input_pdbs/MOL1.pdb"]
}

Run with:
    python run.py config.json --float32_matmul_precision medium
"""

from boltzina_main import Boltzina
from pathlib import Path


def example_usage():
    receptor_pdb = "sample/CDK2/1ckp_cdk2.pdb"
    ligand_files = ["sample/CDK2/ligands/input_pdbs/CDK2_active_0.pdb"]
    output_dir = "example_output"
    plants_config = "sample/CDK2/plants_input.conf"
    plants_bin = None  # auto-detected from PATH or common locations
    work_dir = "sample/CDK2/boltz_results_base"

    if not Path(receptor_pdb).exists():
        print(f"Receptor not found: {receptor_pdb}")
        return

    print("Initializing BOLANTS...")
    boltzina = Boltzina(
        receptor_pdb=receptor_pdb,
        output_dir=output_dir,
        config=plants_config,
        plants_bin=plants_bin,
        work_dir=work_dir,
        fname="1ckp_cdk2",
        input_ligand_name="UNL",
    )

    print("Running PLANTS docking + Boltz-2 affinity scoring...")
    boltzina.run(ligand_files)
    boltzina.save_results_csv()

    df = boltzina.get_results_dataframe()
    if not df.empty:
        print(f"\nTop result:")
        print(df[["ligand_name", "docking_score", "affinity_pred_value"]].head())
    else:
        print("No results generated.")


if __name__ == "__main__":
    example_usage()
