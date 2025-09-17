import pandas as pd
from pathlib import Path
from tqdm.auto import tqdm

combined_folder = Path(".")
all_combined_files = list(combined_folder.glob("combined_*.csv"))
print(all_combined_files)
all_dfs = []
for file in tqdm(all_combined_files):
    # Extract protein id: combined_XXXX_nodes.csv -> XXXX
    protein_id = file.stem.replace("combined_", "")
    
    df = pd.read_csv(file)
    df["protein_id"] = protein_id  # add protein id column
    
    all_dfs.append(df)

# Concatenate into one big dataframe
big_df = pd.concat(all_dfs, ignore_index=True)

# Save final combined dataframe
big_df.to_csv(combined_folder / "all_proteins_combined.csv", index=False)
print(f"✅ Combined dataframe saved to: {combined_folder / 'all_proteins_combined.csv'}")
