import pandas as pd
import numpy as np 



df = pd.read_csv('/home/matt/Proj/MatSciCuration_TBV/Data/Dyes/Dyes_Normallized.csv')


# Group by the 'PSMILES' column and combine the rows using the 'first' function to fill NaNs with available values
combined_df = df.groupby('SMILES').agg(lambda x: x.first_valid_index()).reset_index()

# Replace the index values in the resulting DataFrame with the actual values from the original DataFrame
for column in combined_df.columns:
    if column != 'SMILES':
        combined_df[column] = combined_df[column].apply(lambda x: df.loc[x, column] if not np.isnan(x) else np.nan)

combined_df['non_nan_count'] = combined_df.apply(lambda x: x.count(), axis=1)

# Sort the DataFrame based on the non_nan_count column in descending order
sorted_df = combined_df.sort_values(by='non_nan_count', ascending=False)

# Reset the index and drop the non_nan_count column
sorted_df = sorted_df.reset_index(drop=True).drop(columns=['non_nan_count'])


sorted_df.to_csv('/home/matt/Proj/MatSciCuration_TBV/Data/Dyes/Dyes_Matched.csv',index=False)