import os
from rdkit import Chem
import pandas as pd 
import numpy as np
from psmiles import PolymerSmiles as PS





class GenericNormalizer:


    def NormalizeSmiles(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return Chem.MolToSmiles(mol)
        else:
            return None
        

    def NormalizeSmilesinDataframe(self,df):
        PossibleSMI= ['SMILES','Smiles','smiles']
        columns = df.columns()

        if set(PossibleSMI).intersection(columns):
            smicol = set(PossibleSMI).intersection(columns)
        else:
            raise Exception("The SMILES column should be clearly labaled in the dataframe")
        
        df['NormalizedSmiles'] = df[smicol].apply(self.NormalizeSmiles)
        # Remove rows where SMILES strings couldn't be normalized
        df = df.dropna(subset=['NormalizedSmiles'])

        return df


    def NormalizePSmiles(self,psmiles):
        assert type(psmiles) == str, "PSMILES should be a string"
        if "*" not in psmiles:
            raise ValueError( "The enetered string appears not to be in PSMILES format")
        return PS(psmiles).canonicalize 


    def NormalizeSmilesinDataframe(self,df):
        PossibleSMI= ['PSMILES','PSmiles','Psmiles','psmiles','pSMILES','pSmiles']
        columns = df.columns()

        if set(PossibleSMI).intersection(columns):
            smicol = set(PossibleSMI).intersection(columns)
        else:
            raise Exception("The PSMILES column should be clearly labaled in the dataframe")
        
        df['NormalizedSmiles'] = df[smicol].apply(self.NormalizePSmiles)
        # Remove rows where PSMILES strings couldn't be normalized
        df = df.dropna(subset=['NormalizedPSmiles'])

        return df




        
    
        


