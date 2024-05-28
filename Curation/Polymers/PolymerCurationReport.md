# Report on Data Collection: 



In order to display our data curation pipeline, data was collected on a few different topics relevent to materials science and chemical engineering.

The topics were of personal interest to the researchers, and are also databases that are generally lacking in the sceintific literature.

- Polymers 
- Organic Dyes/ Organic Optical Materials 

In addition to this, we chose to evaluate a few different case studies published under the relevent topics.

The Case studies selected for investigation were:

- **Random forest machine learning models for interpretable X-ray absorption near-edge structure spectrum-property relationships**, https://www.nature.com/articles/s41524-020-00376-6, https://doi.org/10.1038/s41524-020-00376-6
- **Predicting the Band Gaps of Inorganic Solids by Machine Learning**,https://pubs.acs.org/doi/10.1021/acs.jpclett.8b00124,https://doi.org/10.1021/acs.jpclett.8b00124
- **Predicting Synthesizability using Machine Learning on Databases of Existing Inorganic Materials**, https://pubs.acs.org/doi/10.1021/acsomega.2c04856, https://doi.org/10.1021/acsomega.2c04856
  
However, the data collection process here just relies on downloading thier git repositories

All data sources are highlighted below

# Data collection on polymers 

## CROW Database:

This is a knowledgbase on polymersciecne and polymer engineering overall. It contains properties and polymer categories for (**number of polymers that I find**). Data was connected though the use of a web scraper, which can be found in the /PolymerDataCollection/CROW/ folder.

Unfortunately, as of 2024 This database is now defunct. We were lucky to capture it during our data curation process.

## Data from the publication "Machine-learning-assisted discovery of polymers with highthermal conductivity using a molecular design algorithm"

Data can be found on the associated github https://github.com/stewu5/HighTCond_Polymer_iqspr.git

Paper can be found at https://www.nature.com/articles/s41524-019-0203-2, https://doi.org/10.1038/s41524-019-0203-2

This paper contains 38,310 polymeric structure-property associations, as well as polymeric structures. The properties included are:
- Thermal Conductivity 
- glass transition tempertaure
- Melting temperature 
- Density

Much of this information comes from the PolyInfo Database.

## Khazana
Contains 6233 datapoints on computed perperties of polymers such as atomization energy (eV/atom),crystallization tendency (%), band gap chain (eV), band gap bulk (eV), electron affinity (eV), and others.


The data in this project was collected from **Polymer informatics with multi-task learning**, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8085610/, 10.1016/j.patter.2021.100238







# Data Collected on Dyes 

## Deep4chem 
 The publically available version of the Deep4Chem database, used in "Deep Learning Optical Spectroscopy Based on Experimental Database: Potential Applications to Molecular Design" (https://pubs.acs.org/doi/10.1021/jacsau.1c00035)



## Ocelot database

OCELOT is an online archive for
Organic Crystals in Electronic and Light-Oriented Technologies

https://oscar.as.uky.edu/database/

25253 data points on optically active chromophores, including SMILES and other DFT calculated proeprties

## World Dye variety

https://www.worlddyevariety.com/

This is a virtual catalog of dyes that has not been maintained since 2014. In an effort to keep the information contained, we used NLP an dwebscraping to collect a large amount of unstructured data on dyes 