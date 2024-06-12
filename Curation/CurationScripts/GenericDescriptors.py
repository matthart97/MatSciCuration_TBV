"""
This script should be used when one needs to colculcte desciprtors for data curation purposes
"""

import rdkit 
from psmiles import PolymerSmiles as PS
import pandas as pd
import numpy as np
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import List

class GenericFeaturization():


#TODO add functons for rdkit descriptors as well 


    def MorganFingerprintsFromSmiles(self, df: pd.DataFrame, radius: int = 2, n_bits: int = 1024) -> pd.DataFrame:
        # Extract SMILES strings and convert to RDKit molecule objects
        smiles_list = df['NormalizedSMILES'].tolist()
        mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

        # Remove None values (invalid SMILES)
        valid_mol_indices = [i for i, mol in enumerate(mol_list) if mol is not None]
        mol_list = [mol for mol in mol_list if mol is not None]

        n_mols = len(mol_list)

        fingerprint_matrix = np.zeros((n_mols, n_bits))

        for i, mol in enumerate(mol_list):
            fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
            fingerprint_matrix[i] = np.array(fingerprint)

        # Create a DataFrame with fingerprints and filter original DataFrame to only valid molecules
        fingerprint_df = pd.DataFrame(fingerprint_matrix, columns=[f'fp_{i}' for i in range(n_bits)], index=valid_mol_indices)
        filtered_df = df.loc[valid_mol_indices].reset_index(drop=True)

        # Concatenate fingerprints and properties
        result_df = pd.concat([filtered_df, fingerprint_df], axis=1)

        return result_df





    def MorganganFingerprintFromPSmiles(df: pd.DataFrame, radius: int = 2, n_bits: int = 2048) -> pd.DataFrame:
        # Extract SMILES strings and convert to RDKit molecule objects
        psmiles_list = df['NormalizedPSMILES'].tolist()
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
