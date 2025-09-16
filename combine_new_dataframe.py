from pathlib import Path
from tqdm.auto import tqdm
import pandas as pd

# --- folders ---
node_folder = Path("..") / "protein_data_frames"
combined_folder = Path("..") / "combined_dataframe"
combined_folder.mkdir(parents=True, exist_ok=True)

# file lists
ideal_files = list(node_folder.glob("ideal_fixed_*.csv"))
first_files = list(node_folder.glob("first_model_*.csv"))

def extract_protein_id(filename, prefix):
    stem = filename.stem  # e.g. ideal_fixed_2lci_nodes
    return stem.replace(prefix, "")

ideal_ids = {extract_protein_id(f, "ideal_fixed_") for f in ideal_files}
first_ids = {extract_protein_id(f, "first_model_") for f in first_files}
common_ids = ideal_ids.intersection(first_ids)

# --- align function (your version, kept behavior but robustified slightly) ---
def align_chains_and_renumber_per_chain(df_first, df_ideal,
                                       chain_col='Chain',
                                       resid_col='Residue_ID',
                                       resname_col='Residue_Name',
                                       keep_old=True,
                                       verbose=False):
    df1 = df_first.copy()
    df2 = df_ideal.copy()

    if keep_old:
        df1[f'{chain_col}_old'] = df1[chain_col]
        df1[f'{resid_col}_old'] = df1[resid_col]
        df2[f'{chain_col}_old'] = df2[chain_col]
        df2[f'{resid_col}_old'] = df2[resid_col]

    # Unique residues (preserve order encountered in df1 / df2)
    unique_first = df1[[chain_col, resid_col, resname_col]].drop_duplicates(keep='first').reset_index(drop=True)
    unique_ideal = df2[[chain_col, resid_col, resname_col]].drop_duplicates(keep='first').reset_index(drop=True)

    # assign temporary IDs in each unique set (1..N)
    unique_first = unique_first.reset_index().rename(columns={'index': 'tmp_idx_first'})
    unique_first['temp_resid'] = range(1, len(unique_first) + 1)
    unique_ideal = unique_ideal.reset_index().rename(columns={'index': 'tmp_idx_ideal'})
    unique_ideal['temp_resid'] = range(1, len(unique_ideal) + 1)

    if verbose:
        print(f"Unique residues in df_first: {len(unique_first)}; df_ideal: {len(unique_ideal)}")

    # check counts
    if len(unique_first) != len(unique_ideal):
        raise ValueError(f"Unique residue counts differ: df_first has {len(unique_first)}, df_ideal has {len(unique_ideal)}")

    # merge temp_resid back into full frames using the three-key match
    df1 = df1.merge(unique_first[[chain_col, resid_col, resname_col, 'temp_resid']], on=[chain_col, resid_col, resname_col], how='left')
    df2 = df2.merge(unique_ideal[[chain_col, resid_col, resname_col, 'temp_resid']], on=[chain_col, resid_col, resname_col], how='left')

    if df1['temp_resid'].isna().any() or df2['temp_resid'].isna().any():
        raise ValueError("Some rows did not receive temp_resid during merge. Check keys.")

    # map temp_resid -> chain from df1 (reference)
    temp_to_chain = unique_first.set_index('temp_resid')[chain_col].to_dict()
    df2[chain_col] = df2['temp_resid'].map(temp_to_chain)

    # renumber residue IDs sequential within each chain based on temp_resid order
    for df in (df1, df2):
        df[resid_col] = df.groupby(chain_col)['temp_resid'].rank(method='dense').astype(int)

    # drop helper column
    df1 = df1.drop(columns=['temp_resid'])
    df2 = df2.drop(columns=['temp_resid'])

    return df1, df2

# --- main loop ---
for protein_id in tqdm(sorted(common_ids), desc="Proteins"):
    ideal_csv = node_folder / f"ideal_fixed_{protein_id}.csv"
    first_csv = node_folder / f"first_model_{protein_id}.csv"
    combined_csv = combined_folder / f"combined_{protein_id}.csv"

    if combined_csv.exists():
        continue

    try:
        df_ideal = pd.read_csv(ideal_csv)
        df_first = pd.read_csv(first_csv)

        # Apply alignment and per-chain renumbering
        df_first, df_ideal = align_chains_and_renumber_per_chain(df_first, df_ideal,
                                                                 chain_col='Chain',
                                                                 resid_col='Residue_ID',
                                                                 resname_col='Residue_Name',
                                                                 keep_old=False,
                                                                 verbose=False)

        # Normalize Atom_Name (H1 -> H) in first model (as you requested)
        df_first.loc[df_first['Atom_Name'] == 'H1', 'Atom_Name'] = 'H'

        # Ensure Chain dtype consistent
        df_first['Chain'] = df_first['Chain'].astype(str)
        df_ideal['Chain'] = df_ideal['Chain'].astype(str)

        # Prepare df_ideal subset with only join keys + coords (so we don't accidentally pull other cols)
        ideal_subset = df_ideal[['Residue_Name', 'Residue_ID', 'Chain', 'Atom_Name', 'Element', 'X', 'Y', 'Z']].copy()

        # LEFT JOIN: keep all rows from df_first, bring ideal coords when match exists
        df_merged = df_first.merge(
            ideal_subset,
            on=['Residue_Name', 'Residue_ID', 'Chain', 'Atom_Name', 'Element'],
            how='left',
            suffixes=('', '_ideal')  # this will keep X,Y,Z from left as X,Y,Z and right as X_ideal if names clash
        )

        # After merge, if right-hand coords are present they will be named 'X_ideal' etc.
        # If pandas created 'X' and 'X_ideal', we're good. But if it kept 'X' only (no match), create empty ideal cols.
        for c in ['X_ideal', 'Y_ideal', 'Z_ideal']:
            # if the right-hand coord exists under plain 'X' (unlikely because left already has X),
            # check and rename; otherwise ensure column exists
            if c not in df_merged.columns:
                # the right-side columns would be named 'X_ideal' only if merge produced them;
                # If not present, try to find 'X' from right by seeing duplicated columns (pandas keeps left)
                df_merged[c] = pd.NA

        # If merge accidentally created duplicate columns like 'X_x'/'X_y', handle common variants:
        # prefer the _ideal variant if present; if not, try to detect suffixes created by pandas
        for base in ['X', 'Y', 'Z']:
            possible_ideal_names = [f"{base}_ideal", f"{base}_y", f"{base}_right"]
            for pname in possible_ideal_names:
                if pname in df_merged.columns and f"{base}_ideal" not in df_merged.columns:
                    df_merged = df_merged.rename(columns={pname: f"{base}_ideal"})

        # Final column ordering: all original df_first columns, then X_ideal,Y_ideal,Z_ideal
        first_cols = list(df_first.columns)
        # remove any accidental duplicates
        final_cols = first_cols + [c for c in ['X_ideal', 'Y_ideal', 'Z_ideal'] if c in df_merged.columns]

        df_out = df_merged.loc[:, final_cols].copy()

        # Diagnostics
        n_first = len(df_first)
        n_merged = len(df_out)
        n_missing_coords = df_out[['X_ideal','Y_ideal','Z_ideal']].isna().any(axis=1).sum() if {'X_ideal','Y_ideal','Z_ideal'}.issubset(df_out.columns) else None
        # print(f"{protein_id}: first rows={n_first}, output rows={n_merged}, missing ideal coords (rows)={n_missing_coords}")

        # save
        df_out.to_csv(combined_csv, index=False)

    except Exception as e:
        print(f"⚠️ Skipped protein {protein_id} due to error: {e}")

    # break
