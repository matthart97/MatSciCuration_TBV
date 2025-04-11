# Report on the formation of the Dye Datalake and the polymer datalake


The lakes were made via combining the csv files from each data source into a single folder. 

Files in each data lake were then pooled into a single flat csv file and exported for curation.

The File Dyepooling.py does this operation in addition to assigning the source of each data point.


ON a more complicated note, the polymer datalake was constructed using the sources in the data collection report, but all structural representations of the SMILES strings had to first be converted to SMILES format.

Contained in the polyemrs recuration section are several python files for converting everything to PSMILES.