# this code will unpack the rdata found in the high T conductgin polymer paper 


import pyreadr
import csv
import pandas as pd

import os 

PATH = os.getcwd()+"Polymer_lake/"

data = pyreadr.read_r("/../../../DataCollection/PolymerDataCollection/HighTCond_Polymer_iqspr/initial_SMILES.RData")

df = pd.DataFrame(data['smis.uni'])
df.to_csv(PATH+"highTdata.csv")
	