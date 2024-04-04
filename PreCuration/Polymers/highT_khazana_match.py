import pandas as pd



import os 

PATH = os.getcwd()+"Polymer_lake/"

# Read the input dat

# Read the CSV file into a DataFrame
df = pd.read_csv(PATH+'Khazana.csv')

# Rename the 'smiles' column to 'PSMILES'
df.rename(columns={'smiles': 'PSMILES'}, inplace=True)

# Pivot the DataFrame to create a sparse matrix
sparse_matrix = df.pivot(index='PSMILES', columns='property', values='value')

# Fill NaN values with 0


# Write the sparse matrix to a new CSV file
sparse_matrix.to_csv('Polymer_lake/KhazanaMatrix.csv')


# match the polymers in khazan with thhose in highT 


df1 = pd.read_csv(PATH+"highTdata.csv")
df1.rename(columns = {"smis.uni":"PSMILES"})
df0 = pd.read_csv('Polymer_lake/KhazanaMatrix.csv') 



DF = df0.merge(df1, on = 'PSMILES')
DF.to_csv(PATH+"TandKhaz.csv",index = False)