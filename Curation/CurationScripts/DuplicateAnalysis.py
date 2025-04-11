
"""
This calss will allow one to analyze duplicates given a set of criteria



"""

from scipy.spatial.distance import squareform, pdist
from tqdm import tqdm
from rdkit.Chem import AllChem
import pandas as pd 
from rdkit import Chem
from typing import List
import numpy as np
import re
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolToSmiles
from scipy.spatial.distance import pdist, squareform
from tqdm import tqdm








from scipy.spatial.distance import squareform, pdist
from tqdm import tqdm
from rdkit.Chem import AllChem
import pandas as pd 
from rdkit import Chem
from typing import List
import numpy as np
import re

class DuplicateAnalyzer:
    def FindDuplicateStrings(self, StrList, ReturnSuspiciousStrings=True):
        """
        For seeing duplicates in a data set based on names (requires list of names)
        """
        def standardize_name(name):
            return name.strip().lower()

        name_dict = {}

        if ReturnSuspiciousStrings:
            for name in StrList:
                standardized_name = standardize_name(name)
                if standardized_name in name_dict:
                    name_dict[standardized_name].append(name)
                else:
                    name_dict[standardized_name] = [name]
            
            # Find duplicates
            duplicates = {k: v for k, v in name_dict.items() if len(v) > 1}
            return duplicates
        
        else:
            unique_names = []
            duplicates = []

            for name in StrList:
                standardized_name = standardize_name(name)
                if standardized_name in name_dict:
                    name_dict[standardized_name].append(name)
                    # Add the duplicate names to the duplicates list
                    duplicates.extend(name_dict[standardized_name])
                else:
                    name_dict[standardized_name] = [name]
                    unique_names.append(name)

            # Remove duplicates from the unique list
            unique_names = [name for name in unique_names if name not in duplicates]
            return unique_names

    def FindDuplcateSmilesStrings(self, StrList, ReturnSuspiciousStrings=True):
        def canonicalize(smiles):
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return Chem.MolToSmiles(mol)
            return None

        if ReturnSuspiciousStrings:
            canonical_smiles_dict = {}
            suspiciously_similar = []

            for smiles in StrList:
                canonical_smiles = canonicalize(smiles)
                if canonical_smiles:
                    if canonical_smiles in canonical_smiles_dict:
                        canonical_smiles_dict[canonical_smiles].append(smiles)
                        # Add all suspiciously similar SMILES to the list
                        suspiciously_similar.extend(canonical_smiles_dict[canonical_smiles])
                    else:
                        canonical_smiles_dict[canonical_smiles] = [smiles]

            return list(set(suspiciously_similar))
        
        else:
            canonical_smiles_dict = {}
            unique_smiles = []

            for smiles in StrList:
                canonical_smiles = canonicalize(smiles)
                if canonical_smiles:
                    if canonical_smiles not in canonical_smiles_dict:
                        canonical_smiles_dict[canonical_smiles] = smiles
                        unique_smiles.append(smiles)

            return unique_smiles

    def find_duplicate_fingerprints(self, df, return_suspicious_datapoints=True):
        Fingerprint_matrix = df.filter(regex='fp_')
        dist = pd.DataFrame(squareform(pdist(Fingerprint_matrix)))

        # Change matrix to numpy for faster runtime
        Dist = dist.to_numpy()
        flaglist = []
        q = len(Dist)

        # Iterate over top triangle of distance matrix 
        for i in tqdm(range(1, q-1)):
            for j in range(i+1, q):
                z = Dist[i][j]
                if z < 0.0:
                    flaglist.append((i, j))
        
        if return_suspicious_datapoints:
            return flaglist
        else:
            # Create a set of indices to be removed
            indices_to_remove = {j for _, j in flaglist}
            # Return DataFrame with duplicate entries removed
            return df.drop(df.index[list(indices_to_remove)])





















"""
class DuplicateAnalyzer():

    def FindDuplicateStrings(StrList,ReturnSuspiciousStrings=True):
        


        def standardize_name(name,name_dict):
            return name.strip().lower()

        if ReturnSuspiciousStrings:
        
            # Use a dictionary to store standardized names
            



            for name in StrList:
                standardized_name = standardize_name(name)
                if standardized_name in name_dict:
                    name_dict[standardized_name].append(name)
                else:
                    name_dict[standardized_name] = [name]
            
            # Find duplicates
            duplicates = {k: v for k, v in name_dict.items() if len(v) > 1}
            
            return duplicates
        
        else:

            name_dict = {}
            unique_names = []
            duplicates = []

            for name in StrList:
                standardized_name = standardize_name(name)
                if standardized_name in name_dict:
                    name_dict[standardized_name].append(name)
                    # Add the duplicate names to the duplicates list
                    duplicates.extend(name_dict[standardized_name])
                else:
                    name_dict[standardized_name] = [name]
                    unique_names.append(name)

            # Remove duplicates from the unique list
            unique_names = [name for name in unique_names if name not in duplicates]

            return unique_names
        


    def FindDuplcateSmilesStrings(self,StrList, ReturnSuspiciousStrings=True):


            def canonicalize(smiles):
                mol = MolFromSmiles(smiles)
                if mol:
                    return MolToSmiles(mol)
                return None
            if ReturnSuspiciousStrings:
                canonical_smiles_dict = {}
                suspiciously_similar = []

                for smiles in StrList:
                    canonical_smiles = canonicalize(smiles)
                    if canonical_smiles:
                        if canonical_smiles in canonical_smiles_dict:
                            canonical_smiles_dict[canonical_smiles].append(smiles)
                            # Add all suspiciously similar SMILES to the list
                            suspiciously_similar.extend(canonical_smiles_dict[canonical_smiles])
                        else:
                            canonical_smiles_dict[canonical_smiles] = [smiles]

                return list(set(suspiciously_similar))
            

            else:
                canonical_smiles_dict = {}
                unique_smiles = []

                for smiles in StrList:
                    canonical_smiles = canonicalize(smiles)
                    if canonical_smiles:
                        if canonical_smiles not in canonical_smiles_dict:
                            canonical_smiles_dict[canonical_smiles] = smiles
                            unique_smiles.append(smiles)

                return unique_smiles
            

    def find_duplicate_fingerprints(df, return_suspicious_datapoints=True):

        
        Fingerprint_matrix = df.filter(regex='fp_')
        dist = pd.DataFrame(squareform(pdist(Fingerprint_matrix)))

        # Change matrix to numpy for faster runtime
        Dist = dist.to_numpy()
        flaglist = []
        q = len(Dist)

        # Iterate over top triangle of distance matrix 
        for i in tqdm(range(1, q-1)):
            for j in range(i+1, q):
                z = Dist[i][j]
                if z < 0.0:
                    flaglist.append((i, j))
        
        if return_suspicious_datapoints:
            return flaglist
        else:
            # Create a set of indices to be removed
            indices_to_remove = {j for _, j in flaglist}
            # Return DataFrame with duplicate entries removed
            return df.drop(df.index[list(indices_to_remove)])
"""