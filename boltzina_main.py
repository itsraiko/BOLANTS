import os
import sys
import gc
import subprocess

# Ensure the conda env's bin is on PATH so shell=True subprocesses
# (pdb_chain, pdb_tidy, obabel, maxit …) are found.
_conda_bin = os.path.dirname(sys.executable)
if _conda_bin not in os.environ.get("PATH", ""):
    os.environ["PATH"] = _conda_bin + os.pathsep + os.environ.get("PATH", "")

# Add maxit to PATH and set RCSBROOT if not already set
_maxit_root = os.path.join(os.path.dirname(__file__), "maxit-v11.300-prod-src")
_maxit_bin  = os.path.join(_maxit_root, "bin")
if os.path.isdir(_maxit_bin) and _maxit_bin not in os.environ.get("PATH", ""):
    os.environ["PATH"]     = _maxit_bin + os.pathsep + os.environ.get("PATH", "")
    os.environ["RCSBROOT"] = _maxit_root
import pickle
import json
import csv
import copy
import torch
from tqdm import tqdm
from pathlib import Path
from typing import List, Dict, Any, Optional
from multiprocessing import Pool
import pandas as pd
from rdkit import Chem
import shutil
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

from boltzina.data.parse.mmcif import parse_mmcif
from boltzina.affinity.predict_affinity import load_boltz2_model, predict_affinity
from boltz.main import get_cache_path

class Boltzina:
    def __init__(
        self,
        receptor_pdb: str,
        output_dir: str,
        config: str,
        work_dir: Optional[str] = None,
        seed: Optional[int] = None,
        num_workers: int = 4,
        vina_cpu: int = 1,   # kept for API compatibility, unused by PLANTS
        plants_bin: Optional[str] = None,
        batch_size: int = 4,
        num_boltz_poses: int = 1,
        timeout: int = 300,
        vina_override: bool = False,
        boltz_override: bool = False,
        input_ligand_name: str = "MOL",
        base_ligand_name: str = "MOL",
        ligand_chain_id: str = "B",
        fname: Optional[str] = None,
        float32_matmul_precision: str = "highest",
        scoring_only: bool = False,
        skip_docking: bool = False,
        skip_run_structure: bool = True,
        use_kernels: bool = False,
        clean_intermediate_files: bool = True,
        prepared_mols_file: Optional[str] = None,
        predict_affinity_args: Optional[dict] = None,
        pairformer_args: Optional[dict] = None,
        msa_args: Optional[dict] = None,
        steering_args: Optional[dict] = None,
        diffusion_process_args: Optional[dict] = None,
        run_trunk_and_structure: bool = True, # Strongly recommended to be True
    ):
        self.receptor_pdb = Path(receptor_pdb)
        self.output_dir = Path(output_dir)
        self.config = Path(config)
        self.work_dir = Path(work_dir)
        self.seed = seed
        self.vina_override = vina_override
        self.boltz_override = boltz_override
        self.plants_bin = plants_bin or self._find_plants()
        self.num_workers = num_workers
        self.batch_size = batch_size
        self.results = []
        self.num_boltz_poses = num_boltz_poses
        self.pose_idxs = [str(pose_idx) for pose_idx in range(1, self.num_boltz_poses + 1)]
        self.input_ligand_name = input_ligand_name
        self.base_ligand_name = base_ligand_name
        self.ligand_chain_id = ligand_chain_id
        self.float32_matmul_precision = float32_matmul_precision
        self.scoring_only = scoring_only
        self.skip_docking = skip_docking
        self.skip_run_structure = skip_run_structure
        self.use_kernels = use_kernels
        self.timeout = timeout
        self.clean_intermediate_files = clean_intermediate_files
        self.predict_affinity_args = predict_affinity_args
        self.pairformer_args = pairformer_args
        self.msa_args = msa_args
        self.steering_args = steering_args
        self.diffusion_process_args = diffusion_process_args
        self.run_trunk_and_structure = run_trunk_and_structure
        # Create output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.prepared_mols_file = prepared_mols_file
        self.mol_dict = None
        self.vina_cpu = vina_cpu

        self.ligand_files = []
        # Prepare receptor MOL2 file for PLANTS
        if not self.scoring_only:
            self.receptor_mol2 = self._prepare_receptor_mol2()
        else:
            self.receptor_mol2 = self.receptor_pdb

        # Initialize cache directory and CCD
        self.cache_dir = Path(get_cache_path())
        self.ccd_path = self.cache_dir / 'ccd.pkl'
        self.ccd = None
        manifest_path = self.work_dir / "processed" / "manifest.json"
        with open(manifest_path, "r") as f:
            manifest = json.load(f)
        self.manifest = manifest

        self.fname = self._get_fname() if fname is None else fname
        torch.set_float32_matmul_precision(self.float32_matmul_precision)


    def _find_plants(self) -> Optional[str]:
        """Search for the PLANTS binary in PATH and common locations."""
        found = shutil.which("PLANTS") or shutil.which("plants")
        if found:
            return found
        common = [
            "/home/ege/Desktop/zeynep_plants/PLANTS",
            "/usr/local/bin/PLANTS",
            "/opt/PLANTS/PLANTS",
        ]
        for p in common:
            if os.path.isfile(p) and os.access(p, os.X_OK):
                return p
        return None

    def _prepare_receptor_mol2(self) -> Path:
        """Convert receptor PDB to MOL2 format for PLANTS using Open Babel."""
        receptor_mol2 = self.output_dir / "receptor.mol2"
        if receptor_mol2.exists() and not self.vina_override:
            print(f"Skipping receptor preparation for {receptor_mol2} because it already exists")
            return receptor_mol2

        spores = shutil.which("SPORES") or shutil.which("spores")
        if spores:
            # SPORES gives better atom-type assignments for PLANTS
            cmd = [spores, "--mode", "settypes", str(self.receptor_pdb), str(receptor_mol2)]
        else:
            cmd = ["obabel", str(self.receptor_pdb), "-O", str(receptor_mol2)]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Failed to convert receptor to MOL2: {e.stderr}")

        if not receptor_mol2.exists():
            raise RuntimeError(f"Failed to create receptor MOL2 file: {receptor_mol2}")

        return receptor_mol2

    def _load_ccd(self) -> Dict[str, Any]:
        if self.ccd_path.exists():
            with self.ccd_path.open('rb') as file:
                ccd = pickle.load(file)
                if self.base_ligand_name in ccd:
                    ccd.pop(self.base_ligand_name)
                return ccd
        else:
            return {}

    def _get_fname(self) -> str:
        return self.manifest["records"][0]["id"]

    def run(self, ligand_files: List[str]) -> None:
        if self.scoring_only:
            print("Running scoring only...")
            self.run_scoring_only(ligand_files)
            return
        self.ligand_files = ligand_files
        prep_tasks = []
        for idx, ligand_file in enumerate(ligand_files):
            ligand_path = Path(ligand_file)
            ligand_output_dir = self.output_dir / "out" / str(idx)
            ligand_output_dir.mkdir(parents=True, exist_ok=True)
            if (ligand_output_dir / "done").exists() and not self.vina_override:
                continue
            prep_tasks.append((idx, ligand_path, ligand_output_dir))
        print(f"Docking {len(ligand_files)} ligands with {self.num_workers} workers...")

        if not self.skip_docking:
            if self.num_workers == 1:
                for task in tqdm(prep_tasks, desc="Preparing ligands"):
                    self._prepare_ligand(task)
            else:
                ligand_num_workers = self.num_workers // self.vina_cpu
                with Pool(ligand_num_workers) as pool:
                    list(tqdm(pool.imap(self._prepare_ligand, prep_tasks), total=len(prep_tasks), desc="Preparing ligands"))
        else:
            print("Skipping docking...")

        self.ccd = self._load_ccd()
        if self.prepared_mols_file:
            with open(self.prepared_mols_file, "rb") as f:
                self.mol_dict = pickle.load(f)

        for idx, ligand_file in enumerate(ligand_files):
            all_exist = True
            for pose_idx in self.pose_idxs:
                fname = f"{self.fname}_{idx}_{pose_idx}"
                if not (self.output_dir / "boltz_out" / "predictions" / fname / f"affinity_{fname}.json").exists():
                    all_exist = False
                    break
            if not all_exist:
                ligand_output_dir = self.output_dir / "out" / str(idx)
                self._update_ccd_for_ligand(ligand_output_dir, Path(ligand_file))

        print("Preparing structures for scoring...")
        structure_tasks = []
        for idx, ligand_file in enumerate(ligand_files):
            ligand_path = Path(ligand_file)
            ligand_output_dir = self.output_dir / "out" / str(idx)
            docked_ligands_dir = ligand_output_dir / "docked_ligands"
            for pose_idx in self.pose_idxs:
                complex_file = docked_ligands_dir / f"docked_ligand_{pose_idx}_{self.ligand_chain_id}_complex_fix.cif"
                fname = f"{self.fname}_{ligand_output_dir.stem}_{pose_idx}"
                affinity_file = self.output_dir / "boltz_out" / "predictions" / fname / f"affinity_{fname}.json"
                pre_affinity_file = self.output_dir / "boltz_out" / "predictions" / fname / f"pre_affinity_{fname}.npz"
                if affinity_file.exists() and not self.boltz_override:
                    continue
                if pre_affinity_file.exists() and not self.boltz_override:
                    continue
                structure_tasks.append((complex_file, pose_idx, idx))

        print(f"Preparing {len(structure_tasks)} structures with {self.num_workers} workers...")
        (self.output_dir / "boltz_out" / "processed").mkdir(parents=True, exist_ok=True)

        for complex_file, pose_idx, ligand_idx in tqdm(structure_tasks, desc="Preparing structures"):
            self._prepare_structure(complex_file, pose_idx, ligand_idx)

        # clean CCD and mol_dict
        self.ccd = None
        self.mol_dict = None
        gc.collect()

        if self.clean_intermediate_files:
            for complex_file, pose_idx, ligand_idx in tqdm(structure_tasks, desc="Preparing structures"):
                    self._cleanup_preaffinity_intermediates(pose_idx, ligand_idx)

        record_ids = []
        for idx, ligand_file in enumerate(ligand_files):
            ligand_path = Path(ligand_file)
            ligand_output_dir = self.output_dir / "out" / str(idx)
            docked_ligands_dir = ligand_output_dir / "docked_ligands"
            for pose_idx in self.pose_idxs:
                fname = f"{self.fname}_{idx}_{pose_idx}"
                affinity_file = self.output_dir / "boltz_out" / "predictions" / fname / f"affinity_{fname}.json"
                pre_affinity_file = self.output_dir / "boltz_out" / "predictions" / fname / f"pre_affinity_{fname}.npz"
                if affinity_file.exists() and not self.boltz_override:
                    continue
                if not pre_affinity_file.exists(): # If preparation failed
                    continue
                record_ids.append(fname)

        print(f"Affinity prediction will be run for {len(record_ids)} records")
        self._update_manifest(record_ids)
        self._link_constraints(record_ids)

        print("Scoring poses with Boltzina...")
        # Execute scoring with torch multiprocessing
        self._score_poses()
        self._extract_results()

        # Clean up intermediate files after scoring
        if self.clean_intermediate_files:
            self._cleanup_scoring_intermediates()

    def _prepare_ligand(self, task_data):
        """Prepare ligand task for multiprocessing"""
        idx, ligand_path, ligand_output_dir = task_data

        try:
            ligand_mol2    = ligand_output_dir / "ligand.mol2"
            plants_out_dir = ligand_output_dir / "plants_out"
            if not (plants_out_dir / "docked_ligands.mol2").exists():
                # Convert ligand PDB to MOL2
                self._convert_to_mol2(ligand_path, ligand_mol2)

                # Run PLANTS docking
                self._run_plants(ligand_mol2, plants_out_dir)

            # Preprocess docked structures (split poses, build CIF complexes)
            self._preprocess_docked_structures(idx, plants_out_dir)

            # Clean up intermediate files
            if self.clean_intermediate_files:
                self._cleanup_plants_intermediates(ligand_output_dir)

            # Touch done file only on successful completion
            (ligand_output_dir / "done").touch()

        except Exception as e:
            print(f"Error processing ligand {ligand_path} (idx={idx}): {e}")
            # Remove done file if it exists (in case of partial failure)
            done_file = ligand_output_dir / "done"
            if done_file.exists():
                try:
                    done_file.unlink()
                except OSError:
                    pass

    def _convert_to_mol2(self, input_file: Path, output_file: Path) -> None:
        """Convert ligand PDB to MOL2 format for PLANTS using SPORES or Open Babel."""
        if output_file.exists() and not self.vina_override:
            return

        spores = shutil.which("SPORES") or shutil.which("spores")
        if spores:
            cmd = [spores, "--mode", "settypes", str(input_file), str(output_file)]
        else:
            cmd = ["obabel", str(input_file), "-O", str(output_file)]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            if not output_file.exists():
                raise RuntimeError(f"Failed to create MOL2 file: {output_file}")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Error converting {input_file} to MOL2: {e.stderr}")
        except Exception as e:
            raise RuntimeError(f"Unexpected error during MOL2 conversion: {e}")

    def _run_plants(self, ligand_mol2: Path, plants_out_dir: Path) -> None:
        """Run PLANTS docking. Output is written to plants_out_dir/docked_ligands.mol2."""
        if (plants_out_dir / "docked_ligands.mol2").exists() and not self.vina_override:
            print(f"Skipping PLANTS docking for {plants_out_dir} because it already exists")
            return

        if self.plants_bin is None:
            raise RuntimeError(
                "PLANTS binary not found. Set plants_bin= in Boltzina() or add PLANTS to PATH."
            )

        # Remove previous output dir if it exists (PLANTS refuses to overwrite)
        if plants_out_dir.exists():
            shutil.rmtree(plants_out_dir)

        # Write a per-ligand PLANTS config file next to the ligand
        config_path = ligand_mol2.parent / "plants.conf"
        with open(config_path, "w") as f:
            # Copy base settings from the user-supplied config (binding site, scoring, etc.)
            with open(self.config) as base:
                f.write(base.read().rstrip() + "\n\n")
            # Append paths and output dir
            f.write(f"protein_file  {self.receptor_mol2}\n")
            f.write(f"ligand_file   {ligand_mol2}\n")
            f.write(f"output_dir    {plants_out_dir}\n")

        cmd = [self.plants_bin, "--mode", "screen", str(config_path)]
        try:
            subprocess.run(cmd, check=True, timeout=self.timeout)
            if not (plants_out_dir / "docked_ligands.mol2").exists():
                raise RuntimeError(f"PLANTS failed to create docked_ligands.mol2 in {plants_out_dir}")
        except subprocess.TimeoutExpired:
            raise RuntimeError(f"PLANTS docking timed out for {ligand_mol2} after {self.timeout} seconds")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"PLANTS docking failed for {ligand_mol2}: {e.stderr}")
        except Exception as e:
            raise RuntimeError(f"Unexpected error during PLANTS docking: {e}")

    def _preprocess_docked_structures(self, ligand_idx: int, plants_out_dir: Path) -> None:
        ligand_output_dir = self.output_dir / "out" / str(ligand_idx)
        docked_ligands_dir = ligand_output_dir / "docked_ligands"
        docked_ligands_dir.mkdir(exist_ok=True)
        complex_fix_cifs = [
            docked_ligands_dir / f"docked_ligand_{pose_idx}_{self.ligand_chain_id}_complex_fix.cif"
            for pose_idx in self.pose_idxs
        ]

        if all(complex_fix_cif.exists() for complex_fix_cif in complex_fix_cifs) and not self.vina_override:
            return

        # Split PLANTS multi-pose MOL2 into individual PDB files via Open Babel
        docked_mol2 = plants_out_dir / "docked_ligands.mol2"
        cmd = [
            "obabel", str(docked_mol2), "-m", "-O",
            str(docked_ligands_dir / "docked_ligand_.pdb")
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Failed to split PLANTS docked MOL2 to PDB files: {e.stderr}")
        except Exception as e:
            raise RuntimeError(f"Unexpected error during PLANTS MOL2 to PDB conversion: {e}")

        # Process each docked pose
        for pose_idx in self.pose_idxs:
            pdb_file = docked_ligands_dir / f"docked_ligand_{pose_idx}.pdb"
            if not pdb_file.exists():
                continue

            ligand_output_dir = self.output_dir / "out" / str(ligand_idx)
            base_name = f"docked_ligand_{pose_idx}"
            self._process_pose(ligand_output_dir, base_name, pdb_file)

    def _process_pose(self, ligand_output_dir: Path, base_name: str, pdb_file: Path) -> None:
        docked_ligands_dir = ligand_output_dir / "docked_ligands"
        docked_ligands_dir.mkdir(exist_ok=True)
        prep_file = docked_ligands_dir / f"{base_name}_prep.pdb"
        complex_file = docked_ligands_dir / f"{base_name}_{self.ligand_chain_id}_complex.pdb"
        complex_cif = docked_ligands_dir / f"{base_name}_{self.ligand_chain_id}_complex.cif"
        complex_fix_cif = docked_ligands_dir / f"{base_name}_{self.ligand_chain_id}_complex_fix.cif"
        if complex_fix_cif.exists() and not self.vina_override:
            return
        # Process with pdb_chain and pdb_rplresname
        try:
            if self.input_ligand_name != self.base_ligand_name:
                cmd1 = f"pdb_chain -{self.ligand_chain_id} {pdb_file} | pdb_rplresname -\"{self.input_ligand_name}\":{self.base_ligand_name} | pdb_tidy > {prep_file}"
                subprocess.run(cmd1, shell=True, check=True)
            else:
                cmd1 = f"pdb_chain -{self.ligand_chain_id} {pdb_file} | pdb_tidy > {prep_file}"
                subprocess.run(cmd1, shell=True, check=True)

            if not prep_file.exists():
                raise RuntimeError(f"Failed to create prep file: {prep_file}")

            # Merge with receptor
            cmd2 = f"pdb_merge {self.receptor_pdb} {prep_file} | pdb_tidy > {complex_file}"
            subprocess.run(cmd2, shell=True, check=True)

            if not complex_file.exists():
                raise RuntimeError(f"Failed to create complex file: {complex_file}")

            # Convert to CIF
            cmd3 = [
                "maxit", "-input", str(complex_file), "-output", str(complex_cif), "-o", "1"
            ]
            subprocess.run(cmd3, check=True)

            if not complex_cif.exists():
                raise RuntimeError(f"Failed to create CIF file: {complex_cif}")

            # Fix CIF
            cmd4 = [
                "maxit", "-input", str(complex_cif), "-output", str(complex_fix_cif), "-o", "8"
            ]
            subprocess.run(cmd4, check=True)

            if not complex_fix_cif.exists():
                raise RuntimeError(f"Failed to create fixed CIF file: {complex_fix_cif}")

        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Error processing pose {pdb_file}: {e.stderr}")
        except Exception as e:
            raise RuntimeError(f"Unexpected error processing pose {pdb_file}: {e}")

    def _update_ccd_for_ligand(self, ligand_output_dir: Path, ligand_path: Optional[Path] = None):
        base_extra_mols_dir = self.output_dir / "boltz_out" / "processed" / "mols"
        base_extra_mols_dir.mkdir(exist_ok=True, parents=True)
        extra_mols_dir = ligand_output_dir / "boltz_out" / "mols"
        extra_mols_dir.mkdir(exist_ok=True, parents=True)

        if (extra_mols_dir / f"{self.base_ligand_name}.pkl").exists() and not self.vina_override:
            all_exist = True
            for pose_idx in self.pose_idxs:
                fname = f"{self.fname}_{ligand_output_dir.stem}_{pose_idx}"
                if not (base_extra_mols_dir / f"{fname}.pkl").exists():
                    all_exist = False
                    break
            if all_exist:
                return

        if self.mol_dict:
            mol = self.mol_dict[ligand_path.stem]
        else:
            mol = Chem.MolFromPDBFile(ligand_path)
            for atom in mol.GetAtoms():
                pdb_info = atom.GetPDBResidueInfo()
                if pdb_info:
                    atom_name = pdb_info.GetName().strip().upper()
                    atom.SetProp("name", atom_name)

        if mol is None:
            raise ValueError(f"Failed to read PDB file {ligand_path}")

        for pose_idx in self.pose_idxs:
            fname = f"{self.fname}_{ligand_output_dir.stem}_{pose_idx}"
            with open(base_extra_mols_dir / f"{fname}.pkl", "wb") as f:
                pickle.dump({self.base_ligand_name: mol}, f)
        with open(extra_mols_dir / f"{self.base_ligand_name}.pkl", "wb") as f:
            pickle.dump(mol, f)
        return

    def _link_constraints(self, record_ids: List[str]) -> None:
        source_constraints_file = self.work_dir / "processed" / "constraints" / f"{self.fname}.npz"
        target_constraints_dir = self.output_dir / "boltz_out" / "processed" / "constraints"
        target_constraints_dir.mkdir(exist_ok=True, parents=True)
        for record_id in record_ids:
            target_constraints_file = target_constraints_dir / f"{record_id}.npz"
            if not target_constraints_file.exists():
                shutil.copy(source_constraints_file, target_constraints_file)
        return

    def _update_manifest(self, record_ids: List[str]) -> None:
        manifest = copy.deepcopy(self.manifest)
        record = [record for record in manifest["records"] if record["id"] == self.fname][0]
        manifest["records"] = []
        for record_id in record_ids:
            new_record = copy.deepcopy(record)
            # for chain_id, _ in enumerate(new_record["chains"]):
            #     if new_record["chains"][chain_id]["msa_id"] != -1:
            #         new_record["chains"][chain_id]["msa_id"] = f"{record_id}_{chain_id}"
            new_record["id"] = record_id
            manifest["records"].append(new_record)
        with open(self.output_dir / "boltz_out" / "processed" / f"manifest.json", "w") as f:
            json.dump(manifest, f, indent=4)

    def _prepare_structure_parallel(self, task_data):
        complex_file, pose_idx, ligand_idx = task_data
        return self._prepare_structure(complex_file, pose_idx, ligand_idx)

    def _prepare_structure(self, complex_file: Path, pose_idx: str, ligand_idx: int) -> Optional[Path]:
        """Prepare structure by parsing MMCIF and saving structure data"""
        fname = f"{self.fname}_{ligand_idx}_{pose_idx}"
        pose_output_dir = self.output_dir / "boltz_out" / "predictions" / fname
        if not complex_file.exists():
            try:
                plants_out_dir = self.output_dir / "out" / str(ligand_idx) / "plants_out"
                self._preprocess_docked_structures(ligand_idx, plants_out_dir)
            except Exception as e:
                print(f"Error preparing structure for complex {complex_file} and pose {pose_idx}: {e}")
                return None
        output_path = pose_output_dir / f"pre_affinity_{fname}.npz"
        if output_path.exists() and not self.boltz_override:
            print(f"Skipping structure preparation for pose {pose_idx} because it already exists")
            return pose_output_dir
        pose_output_dir.mkdir(parents=True, exist_ok=True)
        extra_mols_dir = self.output_dir / "out" / str(ligand_idx) / "boltz_out" / "mols"

        # Build mols dict: ligand mol loaded directly into dict,
        # moldir set to the boltz cache so standard residues (MET, etc.) resolve correctly.
        mols_for_parsing = dict(self.ccd)
        ligand_mol_pkl = extra_mols_dir / f"{self.base_ligand_name}.pkl"
        if ligand_mol_pkl.exists():
            with open(ligand_mol_pkl, "rb") as _f:
                mols_for_parsing[self.base_ligand_name] = pickle.load(_f)

        try:
            # Parse MMCIF structure
            parsed_structure = parse_mmcif(
                path=str(complex_file),
                mols=mols_for_parsing,
                moldir=str(self.cache_dir / "mols"),  # boltz cache: has all standard CCD residues
                call_compute_interfaces=False
            )

            # Save structure data
            structure_v2 = parsed_structure.data

            structure_v2.dump(output_path)
            assert output_path.exists(), f"Failed to save structure data for pose {pose_idx} at {output_path}"
            return pose_output_dir

        except Exception as e:
            print(f"Error preparing structure for complex {complex_file} and pose {pose_idx}: {e}")
            return None

    def _score_poses(self):
        """Score a single pose"""
        work_dir = self.work_dir
        output_dir = self.output_dir / "boltz_out" / "predictions"
        extra_mols_dir = self.output_dir / "boltz_out" / "processed" / "mols"
        constraints_dir = self.output_dir / "boltz_out" / "processed" / "constraints"
        # Run Boltzina scoring directly with predict_affinity
        self.boltz_model = load_boltz2_model(skip_run_structure = self.skip_run_structure, use_kernels=self.use_kernels, run_trunk_and_structure=self.run_trunk_and_structure, predict_affinity_args=self.predict_affinity_args, pairformer_args=self.pairformer_args, msa_args=self.msa_args, steering_args=self.steering_args, diffusion_process_args=self.diffusion_process_args)
        predict_affinity(
            work_dir,
            model_module=self.boltz_model,
            output_dir=str(output_dir),  # boltz_out directory
            structures_dir=str(output_dir),
            constraints_dir=str(constraints_dir),
            extra_mols_dir=extra_mols_dir,
            manifest_path = self.output_dir / "boltz_out" / "processed" / "manifest.json",
            num_workers=4,
            batch_size=self.batch_size,
            seed = self.seed,
        )

    def _extract_docking_score(self, plants_out_dir: Path, pose_idx: int) -> Optional[float]:
        """Read TOTAL_SCORE for pose_idx (1-based) from PLANTS ranking.csv."""
        ranking_csv = plants_out_dir / "ranking.csv"
        try:
            with open(ranking_csv, newline="") as f:
                reader = csv.DictReader(f)
                for row_num, row in enumerate(reader, start=1):
                    if row_num == pose_idx:
                        return float(row["TOTAL_SCORE"])
        except Exception:
            pass
        return None

    def _extract_results(self):
        results = []
        for ligand_idx, ligand_file in enumerate(self.ligand_files):
            for pose_idx in self.pose_idxs:
                fname = f"{self.fname}_{ligand_idx}_{pose_idx}"
                pose_output_dir = self.output_dir / "boltz_out" / "predictions" / fname
                if not (pose_output_dir / f"affinity_{fname}.json").exists():
                    print(f"Skipping {fname} because it doesn't exist")
                    continue
                with open(pose_output_dir / f"affinity_{fname}.json", "r") as f:
                    affinity = json.load(f)
                affinity["ligand_name"] = ligand_file
                affinity["ligand_idx"] = ligand_idx
                affinity["docking_rank"] = pose_idx
                # Set docking score based on mode
                if self.scoring_only:
                    affinity["docking_score"] = None
                else:
                    plants_out_dir = self.output_dir / "out" / str(ligand_idx) / "plants_out"
                    affinity["docking_score"] = self._extract_docking_score(plants_out_dir, int(pose_idx))
                results.append(affinity)
        self.results = results

    def _cleanup_plants_intermediates(self, ligand_output_dir: Path) -> None:
        """Clean up intermediate files after PLANTS docking, keeping plants_out/."""
        # Remove ligand.mol2 (input to PLANTS, no longer needed)
        ligand_mol2 = ligand_output_dir / "ligand.mol2"
        if ligand_mol2.exists():
            try:
                ligand_mol2.unlink()
            except OSError:
                pass
        docked_ligands_dir = ligand_output_dir / "docked_ligands"
        if docked_ligands_dir.exists():
            for file_path in docked_ligands_dir.iterdir():
                if not file_path.name.endswith(f"_{self.ligand_chain_id}_complex_fix.cif"):
                    if file_path.is_file():
                        file_path.unlink()
                    elif file_path.is_dir():
                        shutil.rmtree(file_path)

    def _cleanup_preaffinity_intermediates(self, pose_idx: str, ligand_idx: int) -> None:
        """Clean up intermediate files after pre-affinity calculation"""
        # Check if pre_affinity file exists in predictions directory
        fname = f"{self.fname}_{ligand_idx}_{pose_idx}"
        predictions_dir = self.output_dir / "boltz_out" / "predictions" / fname
        pre_affinity_file = predictions_dir / f"pre_affinity_{fname}.npz"

        if pre_affinity_file.exists():
            # Get the ligand output directory
            ligand_output_dir = self.output_dir / "out" / str(ligand_idx)

            # Remove boltz_out directory under the ligand output directory
            ligand_boltz_out = ligand_output_dir / "boltz_out"
            if ligand_boltz_out.exists():
                try:
                    shutil.rmtree(ligand_boltz_out)
                except OSError:
                    pass

            # Remove docked_ligands directory
            docked_ligands_dir = ligand_output_dir / "docked_ligands"
            if docked_ligands_dir.exists():
                try:
                    shutil.rmtree(docked_ligands_dir)
                except OSError:
                    pass

    def _cleanup_scoring_intermediates(self) -> None:
        """Clean up intermediate files after scoring"""
        boltz_output_dir = self.output_dir / "boltz_out"

        # Remove processed directory and all its contents
        processed_dir = boltz_output_dir / "processed"
        if processed_dir.exists():
            try:
                import shutil
                shutil.rmtree(processed_dir)
            except OSError:
                pass

        # Remove pre_affinity files from predictions directory
        predictions_dir = boltz_output_dir / "predictions"
        if predictions_dir.exists():
            for pre_affinity_file in predictions_dir.glob("*/pre_affinity_*.npz"):
                # Extract the base name to check for corresponding affinity file
                affinity_file = pre_affinity_file.parent / pre_affinity_file.name.replace("pre_affinity_", "affinity_").replace(".npz", ".json")
                if affinity_file.exists():
                    try:
                        pre_affinity_file.unlink()
                    except OSError:
                        pass


    def save_results_csv(self, output_file: Optional[str] = None) -> None:
        if output_file is None:
            output_file = self.output_dir / "boltzina_results.csv"
        else:
            output_file = Path(output_file)

        if not self.results:
            print("No results to save")
            return

        fieldnames = [
            'ligand_name', 'ligand_idx', 'docking_rank', 'docking_score',
            'affinity_pred_value', 'affinity_probability_binary',
            'affinity_pred_value1', 'affinity_probability_binary1',
            'affinity_pred_value2', 'affinity_probability_binary2'
        ]

        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(self.results)

        print(f"Results saved to {output_file}")

    def get_results_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.results)

    def run_scoring_only(self, ligand_files: List[str]) -> None:
        """
        Run scoring-only mode for ligands with existing poses (no docking).
        Based on scoring_only.py logic.
        """
        self.ligand_files = ligand_files
        print(f"Running scoring-only mode for {len(self.ligand_files)} ligand poses...")

        if self.ccd is None:
            self.ccd = self._load_ccd()

        if self.prepared_mols_file:
            with open(self.prepared_mols_file, "rb") as f:
                self.mol_dict = pickle.load(f)

        # Process pose files
        for ligand_idx, pdb_file in enumerate(self.ligand_files):
            ligand_path = Path(pdb_file)
            ligand_output_dir = self.output_dir / "out" / str(ligand_idx)
            ligand_output_dir.mkdir(parents=True, exist_ok=True)
            base_name = ligand_path.stem
            self._process_pose(ligand_output_dir, base_name, ligand_path)

        # Update CCD for each ligand
        for ligand_idx, pdb_file in enumerate(self.ligand_files):
            ligand_path = Path(pdb_file)
            ligand_output_dir = self.output_dir / "out" / str(ligand_idx)
            all_exist = True
            for pose_idx in self.pose_idxs:
                fname = f"{self.fname}_{ligand_idx}_{pose_idx}"
                if not (self.output_dir / "boltz_out" / "predictions" / fname / f"affinity_{fname}.json").exists():
                    all_exist = False
                    break
            if not all_exist:
                ligand_output_dir = self.output_dir / "out" / str(ligand_idx)
                self._update_ccd_for_ligand(ligand_output_dir, ligand_path)

        # Process boltz input and prepare structures
        record_ids = []
        for ligand_idx, pdb_file in enumerate(self.ligand_files):
            ligand_path = Path(pdb_file)
            base_name = ligand_path.stem
            ligand_output_dir = self.output_dir / "out" / str(ligand_idx)
            complex_file = ligand_output_dir / "docked_ligands" / f"{base_name}_{self.ligand_chain_id}_complex_fix.cif"

            for pose_idx in self.pose_idxs:
                fname = f"{self.fname}_{ligand_output_dir.stem}_{pose_idx}"
                affinity_file = self.output_dir / "boltz_out" / "predictions" / fname / f"affinity_{fname}.json"
                pre_affinity_file = self.output_dir / "boltz_out" / "predictions" / fname / f"pre_affinity_{fname}.npz"
                if affinity_file.exists() and not self.boltz_override:
                    continue
                self._prepare_structure(complex_file, pose_idx, ligand_idx)
                if not pre_affinity_file.exists():
                    continue
                record_ids.append(fname)

        if self.clean_intermediate_files:
            for ligand_idx, pdb_file in enumerate(self.ligand_files):
                for pose_idx in self.pose_idxs:
                    self._cleanup_preaffinity_intermediates(pose_idx, ligand_idx)

        # Update manifest and link constraints
        self._update_manifest(record_ids)
        self._link_constraints(record_ids)

        # Score poses
        print("Scoring poses with Boltzina...")
        self._score_poses()

        # Extract results
        self._extract_results()

        # Clean up intermediate files after scoring
        if self.clean_intermediate_files:
            self._cleanup_scoring_intermediates()



def main():
    import argparse

    parser = argparse.ArgumentParser(description="Boltzina: Vina docking + Boltz scoring pipeline")
    parser.add_argument("--receptor", required=True, help="Receptor PDB file")
    parser.add_argument("--ligands", required=True, nargs="+", help="Ligand files (SDF/MOL2/SMI)")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--config", required=True, help="Vina config file")
    parser.add_argument("--work_dir", help="Working directory for Boltz results")
    parser.add_argument("--vina_override", action="store_true", help="Override existing Vina output directory")
    parser.add_argument("--boltz_override", action="store_true", help="Override existing Boltz output directory")
    parser.add_argument("--num_workers", type=int, default=4, help="Number of workers for parallel processing")
    parser.add_argument("--batch_size", type=int, default=4, help="Batch size for Boltz parallel processing")
    parser.add_argument("--num_boltz_poses", type=int, default=1, help="Number of Boltz poses to score")
    args = parser.parse_args()

    # Initialize Boltzina
    boltzina = Boltzina(
        receptor_pdb=args.receptor,
        output_dir=args.output_dir,
        config=args.config,
        work_dir=args.work_dir,
        vina_override=args.vina_override,
        boltz_override=args.boltz_override,
        num_workers=args.num_workers,
        batch_size=args.batch_size,
        num_boltz_poses=args.num_boltz_poses
    )

    # Run the pipeline
    boltzina.run(args.ligands)

    # Save results
    boltzina.save_results_csv()

    # Print summary
    df = boltzina.get_results_dataframe()
    print(f"\nProcessed {len(df)} poses from {df['ligand_idx'].nunique()} ligands")
    print(f"Best affinity score: {df['affinity_pred_value'].max():.4f}")


if __name__ == "__main__":
    main()
