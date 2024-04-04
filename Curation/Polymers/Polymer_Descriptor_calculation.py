import rdkit 
from psmiles import PolymerSmiles as PS
import pandas as pd
import numpy as np

df =pd.read_csv('/home/matt/Proj/MatSciCuration_TBV/Data/Polymers/PolymerData_Properties_Inspected.csv')

    

def smiles_to_morgan_fingerprints_properties(df: pd.DataFrame, radius: int = 2, n_bits: int = 2048) -> pd.DataFrame:
    # Extract SMILES strings and convert to RDKit molecule objects
    psmiles_list = df['PSMILES'].tolist()
    mol_list = [PS(psmiles) for psmiles in psmiles_list]

    # Remove None values (invalid SMILES)
    valid_mol_indices = [i for i, mol in enumerate(mol_list) if mol is not None]
    mol_list = [mol for mol in mol_list if mol is not None]

    n_mols = len(mol_list)

    fingerprint_matrix = np.zeros((n_mols, n_bits))

    for i, mol in enumerate(mol_list):
        fingerprint = mol.fingerprint('rdkit')
        fingerprint_matrix[i] = np.array(fingerprint)

    # Create a DataFrame with fingerprints and filter original DataFrame to only valid molecules
    fingerprint_df = pd.DataFrame(fingerprint_matrix, columns=[f'fp_{i}' for i in range(n_bits)], index=valid_mol_indices)
    filtered_df = df.loc[valid_mol_indices].reset_index(drop=True)

    # Concatenate fingerprints and properties
    result_df = pd.concat([filtered_df, fingerprint_df], axis=1)

    return result_df


result_df = smiles_to_morgan_fingerprints_properties(df)

result_df.to_csv('/home/matt/Proj/MatSciCuration_TBV/Data/Polymers/PolymerDescriptors.csv',index=False)