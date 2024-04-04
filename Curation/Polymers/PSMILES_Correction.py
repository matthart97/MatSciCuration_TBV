# this code ensures that the sterosk contained in the mines PSMILES is correct 

# This will canonicalize all the PSMILES data I've collected
from psmiles import PolymerSmiles as PS
import pandas as pd 
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import List
import numpy as np

df = pd.read_csv("/home/matt/Proj/MatSciCuration_TBV/Data/Polymers/PolymerData_uncurated.csv")
print(df.columns)
def move_second_asterisk(s):
    count = 0
    for i, char in enumerate(s):
        if char == '*':
            count += 1
            if count == 2:
                new_s = s[:i] + s[i+1:] + '*'
                return new_s
    return s



#run the smiles on a corrected code 
psmi = df['PSMILES']
valid = []
for p in psmi:
    try:
        try: 
            PS(p)
            valid.append(p)
        except:
            p1 = move_second_asterisk(p)
            try:
                PS(p1)
                valid.append(p1)
            except:
                valid.append(None)
    except:
        valid.append(None)

df['corrected'] = valid
df = df.dropna(subset=['corrected'])


df['PSMILES'] = df['corrected'].values

df = df.drop(['corrected'],axis=1)


df.to_csv("/home/matt/Proj/MatSciCuration_TBV/Data/Polymers/PolymerData_PSMILES_corrected.csv",index=False)