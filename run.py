#!/usr/bin/env python3
import json
import sys
import subprocess
import argparse
from boltzina_main import Boltzina
from pathlib import Path


def run_boltz_predict(boltz_yaml: str, work_dir: str, use_msa_server: bool) -> str:
    """Run 'boltz predict' and return the resulting work_dir path.

    boltz writes output to: <out_dir>/boltz_results_<yaml_stem>/
    """
    yaml_path = Path(boltz_yaml)
    if not yaml_path.exists():
        raise FileNotFoundError(f"boltz_yaml not found: {boltz_yaml}")

    fname = yaml_path.stem                          # e.g. "7slz"
    work_path = Path(work_dir)
    # work_dir points to the inner boltz output dir, so the --out_dir
    # passed to boltz is its parent.
    boltz_out_dir = work_path.parent
    boltz_out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [sys.executable, "-m", "boltz", "predict", str(yaml_path),
           "--out_dir", str(boltz_out_dir)]
    if use_msa_server:
        cmd.append("--use_msa_server")

    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    # boltz creates boltz_results_<fname> inside out_dir
    expected = boltz_out_dir / f"boltz_results_{fname}"
    if not expected.exists():
        raise RuntimeError(
            f"boltz predict finished but expected output dir not found: {expected}"
        )
    return str(expected)


def ensure_work_dir(config_dict: dict) -> str:
    """Return a valid work_dir, running boltz predict first if necessary."""
    work_dir = config_dict.get("work_dir", "")
    boltz_yaml = config_dict.get("boltz_yaml", "")
    use_msa_server = config_dict.get("use_msa_server", True)

    manifest_ok = (
        work_dir and
        (Path(work_dir) / "processed" / "manifest.json").exists()
    )

    if manifest_ok:
        return work_dir

    # work_dir is missing or manifest is absent — try to run boltz predict
    if not boltz_yaml:
        if not work_dir:
            raise ValueError(
                "config.json must have either 'work_dir' (pointing to an existing "
                "boltz predict output) or 'boltz_yaml' (to run boltz predict automatically)."
            )
        raise FileNotFoundError(
            f"manifest.json not found in work_dir '{work_dir}/processed/'. "
            "Either run 'boltz predict' manually, or add 'boltz_yaml' to config.json "
            "to let BOLANTS run it automatically."
        )

    # Derive work_dir from yaml stem if not provided
    if not work_dir:
        fname = Path(boltz_yaml).stem
        output_dir = config_dict.get("output_dir", "results")
        work_dir = str(Path(output_dir) / "boltz_results" / f"boltz_results_{fname}")

    print(f"manifest.json not found — running Boltz-2 prediction first...")
    work_dir = run_boltz_predict(boltz_yaml, work_dir, use_msa_server)
    return work_dir


def main():
    parser = argparse.ArgumentParser(description="BOLANTS: PLANTS docking + Boltz-2 affinity")
    parser.add_argument("config", type=str, help="Path to config.json")
    parser.add_argument("--batch_size", type=int, default=1)
    parser.add_argument("--num_workers", type=int, default=1)
    parser.add_argument("--vina_cpu", type=int, default=1, help="Unused (kept for compatibility)")
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--vina_override", action="store_true", help="Re-run PLANTS docking even if output exists")
    parser.add_argument("--boltz_override", action="store_true", help="Re-run Boltz-2 scoring even if output exists")
    parser.add_argument("--use_kernels", action="store_true")
    parser.add_argument("--skip_docking", action="store_true")
    parser.add_argument("--float32_matmul_precision", type=str, default="highest",
                        choices=["highest", "high", "medium"])
    parser.add_argument("--skip_trunk_and_structure", action="store_true")
    parser.add_argument("--keep_intermediate_files", action="store_true")
    parser.add_argument("--output_dir", type=str, default=None)
    args = parser.parse_args()

    with open(args.config, "r") as f:
        config_dict = json.load(f)

    # --- Auto boltz predict if needed ---
    work_dir = ensure_work_dir(config_dict)

    receptor_pdb          = config_dict["receptor_pdb"]
    ligand_files          = config_dict["ligand_files"]
    output_dir            = args.output_dir or config_dict["output_dir"]
    plants_config         = config_dict.get("plants_config") or config_dict.get("vina_config")
    input_ligand_name     = config_dict["input_ligand_name"]
    fname                 = config_dict["fname"]
    seed                  = args.seed if config_dict.get("seed") is None else config_dict["seed"]
    float32_matmul_prec   = config_dict.get("float32_matmul_precision", args.float32_matmul_precision)
    plants_bin            = config_dict.get("plants_bin", None)
    ligand_chain_id       = config_dict.get("ligand_chain_id", "B")
    scoring_only          = config_dict.get("scoring_only", False)
    prepared_mols_file    = config_dict.get("prepared_mols_file", None)
    predict_affinity_args = config_dict.get("predict_affinity_args", None)
    pairformer_args       = config_dict.get("pairformer_args", None)
    msa_args              = config_dict.get("msa_args", None)
    steering_args         = config_dict.get("steering_args", None)
    diffusion_process_args = config_dict.get("diffusion_process_args", None)
    run_trunk_and_structure = not args.skip_trunk_and_structure

    print("--------------------------------")
    print(f"Output directory: {output_dir}")
    print(f"Seed: {seed}")
    print(f"Mode: {'scoring only' if scoring_only else 'docking'}")
    print(f"Using float32 matmul precision: {float32_matmul_prec}")

    boltzina = Boltzina(
        receptor_pdb=receptor_pdb,
        output_dir=output_dir,
        config=plants_config,
        plants_bin=plants_bin,
        work_dir=work_dir,
        input_ligand_name=input_ligand_name,
        fname=fname,
        seed=seed,
        vina_override=args.vina_override,
        boltz_override=args.boltz_override,
        num_workers=args.num_workers,
        vina_cpu=args.vina_cpu,
        batch_size=args.batch_size,
        float32_matmul_precision=float32_matmul_prec,
        ligand_chain_id=ligand_chain_id,
        scoring_only=scoring_only,
        prepared_mols_file=prepared_mols_file,
        use_kernels=args.use_kernels,
        predict_affinity_args=predict_affinity_args,
        pairformer_args=pairformer_args,
        msa_args=msa_args,
        steering_args=steering_args,
        diffusion_process_args=diffusion_process_args,
        skip_docking=args.skip_docking,
        run_trunk_and_structure=run_trunk_and_structure,
        clean_intermediate_files=not args.keep_intermediate_files,
    )

    boltzina.run(ligand_files)

    print("Saving results to CSV...")
    boltzina.save_results_csv()
    print(boltzina.get_results_dataframe())


if __name__ == "__main__":
    main()
