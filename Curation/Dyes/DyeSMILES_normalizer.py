from rdkit import Chem
import pandas as pd 
import numpy as np
def normalize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Chem.MolToSmiles(mol)
    else:
        return None

# Load your dataset
df = pd.read_csv('/home/matt/Proj/MatSciCuration_TBV/Data/Dyes/DyeData_Uncurated.csv')  # assuming you have a pandas DataFrame

# Apply the function to normalize SMILES strings
df['SMILES'] = df['SMILES'].apply(normalize_smiles)

# Remove rows where SMILES strings couldn't be normalized
df = df.dropna(subset=['SMILES'])

df.to_csv('/home/matt/Proj/MatSciCuration_TBV/Data/Dyes/Dyes_Normallized.csv',index=False)