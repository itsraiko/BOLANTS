#!/usr/bin/env python3
"""
BOLANTS prepare.py — Tek komutla tüm girdi dosyalarını hazırla

Kullanım:
    python prepare.py \
        --protein  protein.pdb \
        --smiles   ligands.smi \
        --center   -6.28 10.59 4.83 \
        --output   my_run

Opsiyonel:
    --radius      10.0          (bağlanma cebi yarıçapı, varsayılan 10.0 Å)
    --fname       myprotein     (proje adı, varsayılan: pdb dosyasının adı)
    --chain       A             (hangi protein zinciri, varsayılan: ilk zincir)
    --plants_bin  /path/PLANTS  (PLANTS binary yolu)
    --ligand_name UNL           (ligand residue adı, varsayılan: UNL)

Girdi SMILES dosyası formatı (her satır):
    SMILES  ISIM
    O=C(N...)  ZINC000343638897
    CC(=O)Nc1...  ZINC000012345678

Çıktı (--output klasörü):
    protein_clean.pdb
    protein.yaml
    plants_input.conf
    config.json
    ligands/
        input_pdbs/MOL001.pdb ...
        prepared_mols.pkl
"""

import argparse
import json
import pickle
import subprocess
import sys
from pathlib import Path

# 3-harf → 1-harf amino asit kodu
AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
    # Yaygın modifiye rezidüler
    "MSE": "M", "HSD": "H", "HSE": "H", "HSP": "H",
    "HIE": "H", "HID": "H", "HIP": "H",
    "CYX": "C", "CYM": "C",
}


def clean_protein(pdb_path: Path, out_path: Path, chain: str) -> None:
    """Sadece istenen zincirin ATOM kayıtlarını al."""
    kept = []
    seen_residues = {}  # (chain, resnum, icode) → resname, sıralı tutmak için

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            line_chain = line[21]
            if chain and line_chain != chain:
                continue
            res_name = line[17:20].strip()
            res_num  = line[22:26].strip()
            i_code   = line[26].strip()
            key = (line_chain, res_num, i_code)
            if res_name in AA3TO1:
                kept.append(line.rstrip())
                seen_residues[key] = res_name

    if not kept:
        raise ValueError(
            f"Zincir '{chain}' için ATOM kaydı bulunamadı: {pdb_path}\n"
            f"PDB dosyasındaki zincirler: {_list_chains(pdb_path)}"
        )

    with open(out_path, "w") as f:
        f.write("\n".join(kept) + "\n")

    print(f"  Protein temizlendi: {len(seen_residues)} rezidü, {len(kept)} atom → {out_path.name}")


def _list_chains(pdb_path: Path) -> list:
    chains = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                chains.add(line[21])
    return sorted(chains)


def extract_sequence(pdb_path: Path, chain: str) -> str:
    """ATOM kayıtlarından amino asit dizisini çıkar."""
    seen = {}  # (resnum, icode) → resname  (sıralı, tekrarsız)
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[21] != chain:
                continue
            res_name = line[17:20].strip()
            res_num  = line[22:26].strip()
            i_code   = line[26].strip()
            key = (int(res_num), i_code)
            if res_name in AA3TO1 and key not in seen:
                seen[key] = res_name

    if not seen:
        raise ValueError(f"Zincir '{chain}' için rezidü bulunamadı.")

    ordered = sorted(seen.keys())
    seq = "".join(AA3TO1[seen[k]] for k in ordered)
    print(f"  Dizi çıkarıldı: {len(seq)} amino asit")
    return seq


def read_smiles_file(smiles_path: Path) -> dict:
    """SMILES dosyasını oku: {name: smiles}"""
    mols = {}
    with open(smiles_path) as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            smiles = parts[0]
            name = parts[1] if len(parts) >= 2 else f"MOL{i:04d}"
            mols[name] = smiles
    if not mols:
        raise ValueError(f"SMILES dosyası boş: {smiles_path}")
    print(f"  {len(mols)} ligand okundu.")
    return mols


def write_yaml(out_path: Path, sequence: str, chain_id: str,
               ligand_smiles: str, ligand_chain_id: str = "B") -> None:
    content = f"""\
version: 1
sequences:
  - protein:
      id: {chain_id}
      sequence: {sequence}
  - ligand:
      id: {ligand_chain_id}
      smiles: "{ligand_smiles}"
"""
    out_path.write_text(content)
    print(f"  YAML yazıldı: {out_path.name}")


def write_plants_conf(out_path: Path, cx: float, cy: float, cz: float,
                      radius: float) -> None:
    content = f"""\
# BOLANTS binding site configuration
scoring_function  chemplp
search_speed      speed1

bindingsite_center  {cx:.3f}  {cy:.3f}  {cz:.3f}
bindingsite_radius  {radius:.1f}

cluster_structures  5
cluster_rmsd        2.0
"""
    out_path.write_text(content)
    print(f"  PLANTS config yazıldı: {out_path.name}")


def prepare_ligands(smiles_dict: dict, ligands_dir: Path) -> list:
    """ligand_preparation.py mantığını inline çalıştır."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from boltz.data.parse.schema import compute_3d_conformer
        Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
    except ImportError as e:
        raise ImportError(f"Eksik paket: {e}\n'pip install .' çalıştırdığından emin ol.") from e

    pdbs_dir = ligands_dir / "input_pdbs"
    pdbs_dir.mkdir(parents=True, exist_ok=True)

    mols_dict = {}
    failed = []

    for name, smiles in smiles_dict.items():
        out_pdb = pdbs_dir / f"{name}.pdb"
        try:
            mol2d = Chem.MolFromSmiles(smiles)
            if mol2d is None:
                raise ValueError("SMILES parse hatası")
            mol = Chem.AddHs(mol2d)
            ok = compute_3d_conformer(mol)
            if not ok:
                raise ValueError("3D konformer hesaplanamadı")
            mol = Chem.RemoveHs(mol)

            canonical_order = AllChem.CanonicalRankAtoms(mol)
            for atom, can_idx in zip(mol.GetAtoms(), canonical_order):
                atom_name = atom.GetSymbol().upper() + str(can_idx + 1)
                atom.SetProp("name", atom_name)
                info = Chem.AtomPDBResidueInfo()
                info.SetName(atom_name.rjust(4))
                info.SetResidueName("UNL")
                info.SetResidueNumber(1)
                info.SetChainId("A")
                info.SetIsHeteroAtom(True)
                atom.SetMonomerInfo(info)

            Chem.MolToPDBFile(mol, str(out_pdb))
            mols_dict[name] = mol

        except Exception as e:
            print(f"  [UYARI] {name} hazırlanamadı: {e}")
            failed.append(name)

    pkl_path = ligands_dir / "prepared_mols.pkl"
    with open(pkl_path, "wb") as f:
        pickle.dump(mols_dict, f)

    pdb_paths = [str(pdbs_dir / f"{n}.pdb") for n in mols_dict]
    print(f"  Ligandlar hazırlandı: {len(mols_dict)}/{len(smiles_dict)}"
          + (f"  ({len(failed)} başarısız)" if failed else ""))
    return pdb_paths


def write_config(out_path: Path, cfg: dict) -> None:
    with open(out_path, "w") as f:
        json.dump(cfg, f, indent=4)
    print(f"  config.json yazıldı: {out_path.name}")


def main():
    parser = argparse.ArgumentParser(
        description="BOLANTS — Tüm girdi dosyalarını otomatik hazırla"
    )
    parser.add_argument("--protein",   required=True,  help="Protein PDB dosyası")
    parser.add_argument("--smiles",    required=True,  help="SMILES dosyası (SMILES ISIM formatı)")
    parser.add_argument("--center",    required=True,  nargs=3, type=float,
                        metavar=("X", "Y", "Z"),        help="Bağlanma cebi merkezi (Å)")
    parser.add_argument("--radius",    default=10.0,   type=float, help="Bağlanma cebi yarıçapı (varsayılan: 10.0)")
    parser.add_argument("--output",    default="bolants_run",      help="Çıktı klasörü (varsayılan: bolants_run)")
    parser.add_argument("--fname",     default=None,               help="Proje adı (varsayılan: PDB dosyasının adı)")
    parser.add_argument("--chain",     default=None,               help="Protein zinciri (varsayılan: ilk zincir)")
    parser.add_argument("--plants_bin",default=None,               help="PLANTS binary yolu")
    parser.add_argument("--ligand_name", default="UNL",            help="Ligand residue adı (varsayılan: UNL)")
    args = parser.parse_args()

    protein_path = Path(args.protein).resolve()
    smiles_path  = Path(args.smiles).resolve()
    out_dir      = Path(args.output).resolve()
    cx, cy, cz   = args.center

    if not protein_path.exists():
        sys.exit(f"Hata: protein dosyası bulunamadı: {protein_path}")
    if not smiles_path.exists():
        sys.exit(f"Hata: SMILES dosyası bulunamadı: {smiles_path}")

    # Zinciri otomatik belirle
    chain = args.chain or _list_chains(protein_path)[0]
    fname = args.fname or protein_path.stem.lower()

    out_dir.mkdir(parents=True, exist_ok=True)
    ligands_dir = out_dir / "ligands"

    print(f"\n=== BOLANTS Hazırlık ===")
    print(f"Protein : {protein_path.name}  (zincir {chain})")
    print(f"SMILES  : {smiles_path.name}")
    print(f"Merkez  : ({cx}, {cy}, {cz})  r={args.radius} Å")
    print(f"Çıktı   : {out_dir}\n")

    # 1. Protein temizle
    print("[1/5] Protein temizleniyor...")
    clean_pdb = out_dir / f"{fname}_clean.pdb"
    clean_protein(protein_path, clean_pdb, chain)

    # 2. Dizi çıkar
    print("[2/5] Amino asit dizisi çıkarılıyor...")
    sequence = extract_sequence(clean_pdb, chain)

    # 3. SMILES dosyasını oku
    print("[3/5] Ligandlar okunuyor...")
    smiles_dict = read_smiles_file(smiles_path)
    first_smiles = next(iter(smiles_dict.values()))

    # 4. Dosyaları yaz
    print("[4/5] Yapılandırma dosyaları yazılıyor...")
    yaml_path   = out_dir / f"{fname}.yaml"
    conf_path   = out_dir / "plants_input.conf"
    config_path = out_dir / "config.json"

    write_yaml(yaml_path, sequence, chain, first_smiles)
    write_plants_conf(conf_path, cx, cy, cz, args.radius)

    # 5. Ligand PDB'leri hazırla
    print("[5/5] Ligand 3D konformerleri hesaplanıyor...")
    ligand_files = prepare_ligands(smiles_dict, ligands_dir)

    # config.json
    config = {
        "receptor_pdb":      str(clean_pdb),
        "boltz_yaml":        str(yaml_path),
        "output_dir":        str(out_dir / "results"),
        "plants_config":     str(conf_path),
        "plants_bin":        args.plants_bin,
        "fname":             fname,
        "input_ligand_name": args.ligand_name,
        "prepared_mols_file": str(ligands_dir / "prepared_mols.pkl"),
        "ligand_files":      ligand_files,
    }
    write_config(config_path, config)

    print(f"""
=== Hazırlık tamamlandı! ===

Çalıştırmak için:
    python run.py {config_path} --float32_matmul_precision medium

RTX 4060 gibi düşük VRAM'li GPU için:
    python run.py {config_path} --float32_matmul_precision medium --num_workers 2
""")


if __name__ == "__main__":
    main()
