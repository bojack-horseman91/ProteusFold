from pathlib import Path
from tqdm import tqdm
import pandas as pd
from biopandas.pdb import PandasPdb

# --- configuration: minimal columns + rename mapping ---
KEEP_COLS = [
   'record_name', 'atom_number', 'atom_name', 'alt_loc',
   'residue_name', 'chain_id', 'residue_number', 'insertion',
   'x_coord', 'y_coord', 'z_coord', 'occupancy', 'b_factor', 'element_symbol'
]

RENAME_COLS = {
    'residue_name': 'Residue_Name',
    'residue_number': 'Residue_ID',
    'chain_id': 'Chain',
    'atom_name': 'Atom_Name',
    'element_symbol': 'Element',
    'x_coord': 'X',
    'y_coord': 'Y',
    'z_coord': 'Z',
    'b_factor': 'B_Factor'
}

def pdb_to_csv(pdb_path: Path, csv_path: Path, keep_cols=KEEP_COLS, rename_cols=RENAME_COLS):
    """
    Read a PDB with PandasPdb and write a minimal CSV with renamed columns.
    Missing keep_cols are created as pd.NA so output CSV always has a consistent set.
    """
    ppdb = PandasPdb()
    ppdb.read_pdb(str(pdb_path))
    if 'ATOM' not in ppdb.df or ppdb.df['ATOM'] is None:
        raise ValueError(f"No ATOM dataframe found in {pdb_path}")

    df = ppdb.df['ATOM'].copy()

    # ensure all keep columns exist (fill missing with NA)
    for col in keep_cols:
        if col not in df.columns:
            df[col] = pd.NA

    # select and copy the minimal columns
    out = df[keep_cols].copy()

    # --- replacement for deprecated errors='ignore' ---
    # Try converting whole column to numeric; if conversion raises, leave original values unchanged.
    try:
        out['atom_number'] = pd.to_numeric(out['atom_number'])
    except Exception:
        # keep original column as-is (same behavior as errors='ignore')
        pass

    # rename columns as requested (only columns present will be renamed)
    out = out.rename(columns=rename_cols)

    # write CSV
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(csv_path, index=False)
    return csv_path

# ==== Main script ====
if __name__ == "__main__":
    input_folder = Path(".")
    output_folder = Path("..") / "protein_data_frames"
    output_folder.mkdir(parents=True, exist_ok=True)  # Create folder if it doesn't exist

    # File patterns to match
    patterns = ["ideal_fixed_*.pdb", "first_model_*.pdb"]

    # Loop through files
    all_files = []
    for pattern in patterns:
        all_files.extend(sorted(list(input_folder.glob(pattern))))

    if not all_files:
        print("No PDB files found for patterns:", patterns)

    for pdb_path in tqdm(all_files, desc="PDB -> CSV"):
        csv_path = output_folder / (pdb_path.stem + ".csv")
        try:
            out_file = pdb_to_csv(pdb_path, csv_path)
            # tqdm.write(f"Saved: {out_file}")
        except Exception as e:
            tqdm.write(f"⚠️ Skipped {pdb_path} due to error: {e}")
