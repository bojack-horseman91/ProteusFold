import mdtraj as md
import networkx as nx
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm.auto import tqdm


from mendeleev import element
import pandas as pd

from mendeleev import element

elements = ['C', 'H', 'N', 'O', 'S', 'Se', 'D']  # Add 'D'

data = []
for symbol in elements:
    if symbol == 'D':
        # Get hydrogen (H) properties first
        e = element('H')
        # Copy hydrogen properties and override deuterium-specific ones
        data.append({
            'symbol': 'D',
            'atomic_number': 1,
            'atomic_weight': 2.014,  # Deuterium atomic mass
            'electronegativity_pauling': e.en_pauling,
            'covalent_radius (pm)': e.covalent_radius,
            'vdw_radius (pm)': e.vdw_radius,
            'electron_affinity': e.electron_affinity,
            'valency': e.nvalence('H'),
            'density (g/cm3)': 0.180,  # Density of deuterium gas at STP
            'dipole_polarizability': e.dipole_polarizability,
            'thermal_conductivity': 0.174,  # Approximate for D2 gas
            'neutrons': 1,  # Deuterium has 1 neutron
            'specific heat capacity': 5.21,  # J/g·K, D2 gas approx
            # 'electrophilicity': None,  # No defined electrophilicity for atom
            # 'lattice_constant': None,    # Not applicable for atomic D
            'atomic_radius': e.atomic_radius,
        })
    else:
        e = element(symbol)
        data.append({
            'symbol': e.symbol,
            'atomic_number': e.atomic_number,
            'atomic_weight': e.atomic_weight,
            'electronegativity_pauling': e.en_pauling,
            'covalent_radius (pm)': e.covalent_radius,
            'vdw_radius (pm)': e.vdw_radius,
            'electron_affinity': e.electron_affinity,
            'valency': e.nvalence('H'),
            'density (g/cm3)': e.density,
            'dipole_polarizability': e.dipole_polarizability,
            'thermal_conductivity': e.thermal_conductivity,
            'neutrons': e.neutrons,
            'specific heat capacity': e.specific_heat_capacity,
            # 'electrophilicity': e.electrophilicity(),
            # 'lattice_constant': e.lattice_constant,
            'atomic_radius': e.atomic_radius,
        })

atomic_info = pd.DataFrame(data).set_index('symbol')


print(atomic_info)
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd

aas = 'ACDEFGHIKLMNOPQRSTUVWY'
data = []

for aa in aas:
    analysed = ProteinAnalysis(aa)
    try:
        # special handling for U and O
        if aa == 'U':
            gravy_value = 2.5  # Selenocysteine, approximated
        elif aa == 'O':
            gravy_value = -1.5
        else:
            gravy_value = analysed.gravy()

        data.append({
            'amino_acid': aa,
            'molecular_weight': analysed.molecular_weight(),
            'aromaticity': analysed.aromaticity(),
            'gravy (hydrophobicity)': gravy_value,
            'isoelectric_point': analysed.isoelectric_point(),
            'Charge': analysed.charge_at_pH(7.0),
            'Charge-ph4': analysed.charge_at_pH(4.0),
            'Charge-ph10': analysed.charge_at_pH(10.0),
        })
    except Exception as e:
        print(f"⚠️ Error processing amino acid '{aa}': {e}")
        raise ValueError(f"{aa} caused an error: {e}")


residue_info = pd.DataFrame(data).set_index('amino_acid')
print(residue_info)


three_to_one = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V','SEC':
                'U','PYL':'O'}

aromatic_residues = {'PHE','TYR','TRP','HIS'}
aromatic_atoms = {
  'PHE': {'CG','CD1','CD2','CE1','CE2','CZ'},
  'TYR': {'CG','CD1','CD2','CE1','CE2','CZ'},
  'TRP': {'CG','CD1','CD2','CE2','CE3','CZ2','CZ3','CH2'},
  'HIS': {'CG','ND1','CD2','CE1','NE2'}
}
from sklearn.preprocessing import MinMaxScaler

def normalize_df(df):
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    scaler = MinMaxScaler()
    df_scaled = df.copy()
    df_scaled[numeric_cols] = scaler.fit_transform(df[numeric_cols])
    return df_scaled

# Normalize atomic_info and residue_info
atomic_info_normalized = normalize_df(atomic_info)
residue_info_normalized = normalize_df(residue_info)

# (Optional) check
print("\n✅ Normalized atomic_info:")
print(atomic_info_normalized.head())

print("\n✅ Normalized residue_info:")
print(residue_info_normalized.head())
atomic_info_normalized.to_csv('atomic_properties.csv')
residue_info_normalized.to_csv('residue_properties.csv')
