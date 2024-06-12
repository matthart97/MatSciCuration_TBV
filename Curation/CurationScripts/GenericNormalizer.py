import os
from rdkit import Chem
import pandas as pd 
import numpy as np
from psmiles import PolymerSmiles as PS





class GenericNormalizer:
 



    def NormalizeSmiles(self, smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
        except:
            return None
        if mol is not None:
            try:
                return Chem.MolToSmiles(mol)
            except:
                return None
        else:
            return None
        

    def NormalizeSmilesinDataframe(self,df):
        PossibleSMI= ['SMILES','Smiles','smiles']
        columns = df.columns

        if set(PossibleSMI).intersection(columns):
            smicol= next((col for col in PossibleSMI if col in columns), None)
        else:
            raise Exception("The SMILES column should be clearly labaled in the dataframe")
        
        df['NormalizedSmiles'] = df[smicol].apply(self.NormalizeSmiles)
        # Remove rows where SMILES strings couldn't be normalized
        df = df.dropna(subset=['NormalizedSmiles'])
        # eliminate those with unsepcified elements
        df = df[~df[smicol].str.contains(r'\*', na=False)]

        return df


    def NormalizePSmiles(self,psmiles):
        assert type(psmiles) == str, "PSMILES should be a string"
        if "*" not in psmiles:
            raise ValueError( "The enetered string appears not to be in PSMILES format")
        return PS(psmiles).canonicalize 


    def NormalizePSmilesinDataframe(self,df):
        PossibleSMI= ['PSMILES','PSmiles','Psmiles','psmiles','pSMILES','pSmiles']
        columns = df.columns

        if set(PossibleSMI).intersection(columns):
            smicol = set(PossibleSMI).intersection(columns)
        else:
            raise Exception("The PSMILES column should be clearly labaled in the dataframe")
        
        df['NormalizedSmiles'] = df[smicol].apply(self.NormalizePSmiles)
        # Remove rows where PSMILES strings couldn't be normalized
        df = df.dropna(subset=['NormalizedPSmiles'])

        return df




        
    
        


