# Curation Report for the Creation of the "Dyebrary"


## Definition of the Problem

We seek to create a database of dyes and dye-like molecules. The database should serve as the foundation for any data-driven tool for the contriustion of technologies centered around the optical properties of organic molecules. Therefore, we want to integrate several databasee that describe abroad range of optical properties of compounds. 



## Data Sources and Data Collection Process

The sources for the dyebrary consist of 3 main sources:

- [OCELOT](https://oscar.as.uky.edu/)
  - Contains 25,000 compounds relating to the engineering of organic crystals for optical applications
- [Deep4Chem](http://deep4chem.korea.ac.kr/)
  - Contains over 20,00 dye-solvent combinations and the properties revelvent to spectroscopy 
- [World dye variety](https://www.worlddyevariety.com/) 
  - Contains a sereis of product sheets on about 3000 dyes
- Product sheets frrom various companies
  

### OCELOT

The OCELOT database was downlaoded directly from the website 

### Deep4Chem

Downloaded directly from the source 

### World Dye Variety 

Used a python web scraper to create structured csv files about each dye givien what was on the website. Language model was used to extract information from text.


## Data Curation Process


### Dye Data Pooling
Using the script Dyepooling.py, the data was consolidataed in a sinlge csv file. 


### Dye Curation


For more details on running the code, there is a ipynb notebook displaying the entire curation process under Curation/Dyes/CurationDisplayDyes.ipynb


The final database and result of the curation display is the Dyebrary. Found at /Dyebrary/Dyebrary.csv

