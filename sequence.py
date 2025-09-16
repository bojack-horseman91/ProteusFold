import PeptideBuilder
from Bio.PDB import PDBParser, Structure, Model, PDBIO
from openmm.app import PDBFile
from pdbfixer import PDBFixer
from pathlib import Path
from tqdm.auto import tqdm


# Map three-letter → one-letter
_THREE2ONE = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V',

    # 21st and 22nd amino acids
    'SEC': 'U',  # selenocysteine
    'PYL': 'O',  # pyrrolysine
}

def extract_sequence_and_chainlist(pdb_path):
    """
    Returns:
        seq: full one-letter sequence
        chain_ids: list of chain IDs per residue, in same order
    """
    parser = PDBParser(QUIET=True, PERMISSIVE=True)
    structure = parser.get_structure('X', str(pdb_path))
    model = next(structure.get_models())

    three_letter = []
    chain_ids    = []
    for chain in model:
        for res in chain:
            if res.id[0]==' ' and res.get_resname() in _THREE2ONE:
                three_letter.append(res.get_resname())
                chain_ids.append(chain.id)

    seq = ''.join(_THREE2ONE[r] for r in three_letter)
    return seq, chain_ids

def assign_chain_ids(struct, chain_ids):
    """
    Given a PeptideBuilder structure with 1 chain,
    rename residues so each has its correct chain_id.
    """
    model = next(struct.get_models())
    chain = next(model.get_chains())
    residues = list(chain.get_residues())

    if len(residues) != len(chain_ids):
        raise ValueError(f"Chain length mismatch: built {len(residues)} vs chain_ids {len(chain_ids)}")

    # detach residues from old chain, create multiple chains grouped by chain_id
    from Bio.PDB import Chain
    new_chains = {}
    for res, cid in zip(residues, chain_ids):
        if cid not in new_chains:
            new_chains[cid] = Chain.Chain(cid)
        new_chains[cid].add(res)

    # build new model
    new_model = Model.Model(0)
    for c in new_chains.values():
        new_model.add(c)

    # new structure
    new_struct = Structure.Structure("full")
    new_struct.add(new_model)
    return new_struct
from Bio.PDB import Structure

def save_pdb_with_ter(structure: Structure, out_pdb: str):
    """
    Save a Biopython Structure to PDB format with:
      - Proper fixed‐width columns per PDB spec.
      - Atom serial numbers assigned sequentially.
      - TER record after each chain.
      - Final END record.
    """
    lines = []
    atom_serial = 1

    # Iterate over all models (or just the first if you prefer)
    for model in structure:
        for chain in model:
            last_residue = None
            for residue in chain:
                last_residue = residue
                resname = residue.get_resname()
                resseq  = residue.id[1]
                for atom in residue:
                    name      = atom.get_name()
                    x, y, z   = atom.coord
                    occupancy = atom.get_occupancy() or 1.00
                    bfactor   = atom.get_bfactor()   or 0.00
                    element   = atom.element.strip().upper().rjust(2)

                    # Build fixed-width ATOM line with an explicit blank altLoc (col 17)
                    line = (
                        f"ATOM  {atom_serial:5d} "   # cols 1–12
                        f"{name:^4}"                # cols 13–16
                        " "                         # col 17 = altLoc blank
                        f"{resname:>3}"             # cols 18–20
                        " "                         # col 21
                        f"{chain.id}"               # col 22
                        f"{resseq:4d}"              # cols 23–26
                        "    "                      # cols 27–30 blank
                        f"{x:8.3f}"
                        f"{y:8.3f}"
                        f"{z:8.3f}"
                        f"{occupancy:6.2f}"
                        f"{bfactor:6.2f}"
                        "          "               # cols 67–76 blank
                        f"{element:>2}"             # cols 77–78
                        "\n"
                    )
                    lines.append(line)
                    atom_serial += 1

            # Emit TER record after the last residue of this chain
            if last_residue is not None:
                resname = last_residue.get_resname()
                resseq  = last_residue.id[1]
                ter_line = (
                    f"TER   {atom_serial:5d}      "
                    f"{resname:>3}"   # cols 18–20
                    " "               # col 21
                    f"{chain.id}"     # col 22
                    f"{resseq:4d}"    # cols 23–26
                    + " " * (80 - 26) # pad to 80 cols
                    + "\n"
                )
                lines.append(ter_line)
                atom_serial += 1

    # Final END record
    lines.append("END".ljust(80) + "\n")

    # Write out the PDB file
    with open(out_pdb, 'w') as f:
        f.writelines(lines)

    print(f"💾 Saved structure with TER records to: {out_pdb}")




# def build_peptide_structure(seq, chain_id,
#                             phi=-135.0, psi=135.0, omega=180.0):
#     """
#     Build a single‐chain Structure for the sequence, then rename
#     its chain to chain_id.
#     """
#     if len(seq) < 1:
#         raise ValueError("Empty sequence!")
#     L = len(seq)
#     phi_vals   = [phi]   * (L-1)
#     psi_vals   = [psi]   * (L-1)
#     omega_vals = [omega] * (L-1)
#     struct = PeptideBuilder.make_structure(seq, phi_vals, psi_vals, omega_vals)

#     # Grab the one chain and rename
#     model = next(struct.get_models())
#     chain = next(model.get_chains())
#     chain.id = chain_id

#     # Wrap into fresh Structure object
#     new_struct = Structure.Structure("full")
#     new_model  = Model.Model(0)
#     new_model.add(chain)
#     new_struct.add(new_model)
#     return new_struct

def save_structure(structure, out_pdb):
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out_pdb))
    print(f"✅ Wrote: {out_pdb.name}")

def add_hydrogens(in_pdb, out_pdb):
    try:
        fixer = PDBFixer(filename=str(in_pdb))
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(True)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=7.0)
        with open(out_pdb, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
        print(f"✅ Hydrogens added: {out_pdb.name}")
    except Exception as e:
        print(f"⚠️ Skipping hydrogens for {in_pdb.name}: {e}")

import re
import re

import re
from pathlib import Path


def check_coordinate_range(structure):
    """
    Check the range of coordinates in the structure and warn if they exceed typical PDB limits.
    """
    min_coord = [float('inf')] * 3
    max_coord = [float('-inf')] * 3
    for atom in structure.get_atoms():
        coord = atom.coord
        for i in range(3):
            if coord[i] < min_coord[i]:
                min_coord[i] = coord[i]
            if coord[i] > max_coord[i]:
                max_coord[i] = coord[i]
    print(f"Coordinate range: x={min_coord[0]:.3f} to {max_coord[0]:.3f}, "
          f"y={min_coord[1]:.3f} to {max_coord[1]:.3f}, z={min_coord[2]:.3f} to {max_coord[2]:.3f}")
    for i in range(3):
        if min_coord[i] < -999.999 or max_coord[i] > 9999.999:
            print(f"⚠️ Warning: Coordinates out of standard PDB range for axis {i}: "
                  f"{min_coord[i]:.3f} to {max_coord[i]:.3f}")
if __name__ == "__main__":
    folder = Path('.')
    pdb_files = list(folder.glob('*.pdb'))  # or use '*.pdb' if needed

    for pdb_path in tqdm(pdb_files, desc="Processing PDB files"):
        protein_id = pdb_path.stem  # '1abc' from '1abc.pdb'
        ideal_pdb = pdb_path.with_name(f'ideal_{pdb_path.name}')
        ideal_fixed = pdb_path.with_name(f'ideal_fixed_{pdb_path.name}')

        # Skip if already processed
        if ideal_pdb.exists() and ideal_fixed.exists():
            print(f"⏭️  Skipping {pdb_path.name}, already processed.")
            continue

        # Skip ideal files
        if pdb_path.name.startswith(('ideal_', 'ideal_fixed_')):
            continue

        try:
            print(f"\n🔹 Processing {pdb_path.name}")
            seq, chain_ids = extract_sequence_and_chainlist(pdb_path)
            print(f"{pdb_path.name}: seq_length={len(seq)}, unique_chain_ids={set(chain_ids)}")

            # Build plain peptide
            L = len(seq)
            phi   = [-135.0] * (L - 1)
            psi   = [ 135.0] * (L - 1)
            omega = [180.0] * (L - 1)
            struct = PeptideBuilder.make_structure(seq, phi, psi, omega)

            # Assign real chain IDs
            struct_fixed = assign_chain_ids(struct, chain_ids)

            # Check coordinate range before saving
            check_coordinate_range(struct_fixed)

            # Save with TER
            save_pdb_with_ter(struct_fixed, ideal_pdb)

            # No need for fix_pdb_coordinates; save_pdb_with_ter already handles it
            add_hydrogens(ideal_pdb, ideal_fixed)


        except Exception as e:
            print(f"❌ Error processing {pdb_path.name}: {e}")
            continue


