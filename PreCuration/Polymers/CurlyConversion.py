import json

import json


import pandas as pd
import os


PATH = os.getcwd() + 'Polymer_lake/'


def convert_curly_smiles_to_psmiles(curly_smiles):
    psmiles = curly_smiles.replace("{n+}", "*")
    index = psmiles.find("{-}")
    if index != -1:
        psmiles = psmiles[:index - 1] + '*' + psmiles[index - 1] + psmiles[index + 3:]
    psmiles = psmiles.replace("(n+}", "*")
    return psmiles
# Replace this with the path to your JSON file
json_file = 'Crow_scraped.json'

with open(json_file, 'r') as file:
    data = json.load(file)

for item in data:
    curly_smiles = item["NAMES AND IDENTIFIERS OF POLYMER"][5]["CurlySMILES"]
    psmiles = convert_curly_smiles_to_psmiles(curly_smiles)
    item["NAMES AND IDENTIFIERS OF POLYMER"][5]["PSMILES"] = psmiles
    del item["NAMES AND IDENTIFIERS OF POLYMER"][5]["CurlySMILES"]

# Write the modified data back to the JSON file
with open('CROW_PSMILES.json', 'w') as file:
    json.dump(data, file, indent=4)
