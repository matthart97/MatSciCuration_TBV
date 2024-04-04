
# Read the first dataset
import pandas as pd
import os


PATH = os.getcwd()
data_path = PATH +"Dye_lake/"

df1 = pd.read_csv(PATH+"deep4chem.csv") #changed chromophre to SMILES
df2 = pd.read_csv(PATH+"ocelot_data.csv")
df3 = pd.read_csv(PATH+"WDVDataComplete_uncurated.csv")

df1['source'] = 'D4C' #deep4chem
df2['source'] = 'OCT'
df3['source'] = 'WDV' 
 
# Add column names for the second dataset
df2.columns = [
    'identifier', 'VerticalIonizationE', 'AdiabaticIonizationE', 'VerticalElectronAffinity', 'AdiabaticElectronAffinity', 'HLGap', 'lowest-lyingSinglet', 'lowest-lyingTriplet', 'HoleReorganizationE', 'CatIonRelax2', 'CationRelax1', 'ElectronReorganizationE', 'AnIonRelax1', 'AnionRelax2', 'lumo', 'homo', 'SMILES'
]




# Merge the two datasets on the 'SMILES' column
merged_df = df1.merge(df2, on='SMILES', how='outer')

merged_df = merged_df.merge(df3,on='SMILES', how='outer')
# Write the merged dataset to a CSV file
merged_df.to_csv('DyeData_Uncurated.csv', index=False)