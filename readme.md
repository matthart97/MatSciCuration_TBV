# Trust, but Verify.

"Healthy distrust makes a good basis for cooperation"- Joseph Stalin


This repository should act as ageneral guide to how one should go about collecting and curating data in the materials and chemical sciences.
We hope that this guide will inspire people from across the sceintific informatics disciplines to adopt a standard, reproducible pipeline for all of thier studies. Furthermore, we hope that this repository should encourage people to upload their data curation pipline to a github repository when a sceintific manuscript is submitted to a major journal.

We use the case studies of dyes and polymers as an example of our datacuration pipline

# Overall Pipeline

The Process begins with designing and formally defining your data curation process to ensure reproducibility.

Data Collection
Clearly define the research problem
Decide which data sources to collect from
Decide how data will be collected (direct download, SQL query, etc.) and take note.
Construction of a Data Lake
Consolidate all collected data into one location
Create schema for keeping track of the origin of each datapoint
Workflow Design
Design the Workflow to answer the following questions:

How will the data from different sources be integrated into one database?
How will chemical structures be stored? Common choices include SMILES, InChI, sdf formats.
How will information on experimental conditions be stored?
How should chemical characterization information be handled?
How will the properties of interest be processed and stored?
Validation Rules:
Vaildation rules should be defined after the curation workflow has been constructed. It should answer the following

What elements will require attention during curation (processing data, structure, etc.)?
What defines a “bad” datapoint? How does one know when a datapoint requires attention?
Are there any restrictions that should be placed on the chemical formulations included? For instance, one may only want to include inorganic formulas in their data.
Are there any physical equations associated with the property you are trying to predict? Do these place any logical or physical restrictions on your datapoints?
What software will be used to curate the data? Do custom scripts need to be written?
Will the data be manually inspected?
If the data is too large to be manually inspected, can the inspection process be automated?
Curation Reports
Creating a GitHub repository for data curation is an excellent way to maintain curation version control. Initially, the data lake and a pre-curation report should be uploaded to GitHub. The report should describe all the decisions made by the curator and should at minimum include:

The number of raw data points acquired.
The data sources and which datapoints were taken from them.
All meta-data collected from experiments or data processing steps.
An explanation of the validation rules.
The intended final database schema.
Details of the pipeline for data curation.
Reccommended Data Curation Workflow:
Normalize all chemical structures to the same representation format.
Assign a unique ID to each chemical concept in the data. (Name to structure-property correspondence)
Correct and misrepresented compounds.
Check all quantitative and qualitative endpoints for reasonableness and phyicality
Simple featurization
Duplicate analysis (via names, structures, and features)
Calculate modelability indices
Novelty detaction
Curation Report
Case studies are detailed in this repository. Resulting databases will be placed here in addision to their own repositories. 


# **A note from the authors**


While one can strive for perfection, it's impossible to meet this standard for data curation. A accurate database is one that is contunually updated and monitored.
Mistakes can happen at any stage of the curation process. Therefore, we ask that community members point out and make efforts to correct any errors in the databasescreated by this work, should they be found :)
