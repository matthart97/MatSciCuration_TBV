
# Read the first dataset
import pandas as pd
import os


PATH = os.getcwd()
data_path = PATH +"/Dye_lake/"

df1 = pd.read_csv(data_path+"deep4chem.csv") #changed chromophre to SMILES

exclude = ['SMILES','Tag']
df1.rename(columns=lambda x: x + ('SolventDependent') if x not in exclude else x, inplace=True)
df1 = df1.drop(['Tag'], axis =1)


df2 = pd.read_csv(data_path+"ocelot_data.csv")


df3 = pd.read_csv(data_path+"WDVDataComplete_uncurated.csv")

df3.replace({1.0: 'yes', 0.0: 'unknown', -1.0: 'no'}, inplace=True)

df1['source'] = 'D4C' #deep4chem
df2['source'] = 'OCT' #ocelot 
df3['source'] = 'WDV' #world dye variety
 

print(df2.columns)
# Add column names for the second dataset
df2.columns = [
    'identifier', 'VerticalIonizationE', 'AdiabaticIonizationE', 'VerticalElectronAffinity', 'AdiabaticElectronAffinity', 'HLGap', 'lowest-lyingSinglet', 'lowest-lyingTriplet', 'HoleReorganizationE', 'CatIonRelax2', 'CationRelax1', 'ElectronReorganizationE', 'AnIonRelax1', 'AnionRelax2', 'lumo', 'homo', 'SMILES', 'source'
]


# Merge the two datasets on the 'SMILES' column
merged_df = df1.merge(df2, on='SMILES', how='outer')

merged_df = merged_df.merge(df3,on='SMILES', how='outer')
# Write the merged dataset to a CSV file
merged_df.to_csv('DyeDataUncurated.csv', index=False)