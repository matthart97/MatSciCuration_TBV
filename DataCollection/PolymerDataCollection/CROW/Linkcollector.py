import requests
from bs4 import BeautifulSoup
import json
import requests
import csv
#url = 'https://www.polymerdatabase.com/polymer%20index/polyacrylamides.html'
#page = requests.get(url)
#html = page.text


# this code will collect all the links needed to scape polymer data from the CROW database


# first, collect all the links that lead to types of polymers. This collection is small enough to do manually

Upper_links= ['https://www.polymerdatabase.com/polymer%20index/home.html','https://www.polymerdatabase.com/polymer%20index/C-D%20index.html',
'https://www.polymerdatabase.com/polymer%20index/E-F%20index.html',
'https://www.polymerdatabase.com/polymer%20index/G-L%20index.html',
'https://www.polymerdatabase.com/polymer%20index/G-L%20index2.html',
'https://www.polymerdatabase.com/polymer%20index/M-P%20index.html',
'https://www.polymerdatabase.com/polymer%20index/M-P%20index2.html',
'https://www.polymerdatabase.com/polymer%20index/S-V%20index.html',
'https://www.polymerdatabase.com/polymer%20index/S-V%20index2.html']


# the next step is to collect all of the links, corresponding to polymer families, directed in each of the upper links



#empty list to hold families
family_links = []


for i in Upper_links:
    #establish bs4 html 
    page = requests.get(i)

    html = page.text
    soup = BeautifulSoup(html, "html.parser")
    
    links = soup.find_all("a")
    for link in links:
        href = link.get("href")
        if href.startswith("poly") or href.startswith("../poly") or href.startswith("Poly") or href.startswith("epoxy"):
            family_links.append(href)
            

# store all the polymer family links in a different file 
with open("family_links.csv", "w", newline="", encoding="utf-8") as csvfile:
    csv_writer = csv.writer(csvfile)
    

    for poly_link in family_links:
        csv_writer.writerow([poly_link])



# read in the file containting the family links and make them amenable to web scaping

with open("/home/matt/Proj/MatSciCuration_TBV/DataCollection/PolymerDataCollection/CROW/family_links.csv") as f:
    lines = f.readlines()

newlinks = []
for i in range(0,len(lines)):
    if "../polymer index/" in lines[i]:
        g = lines[i].replace("../polymer index/","")
        newlinks.append(g)
    else:
        newlinks.append(lines[i])
        

with open("/home/matt/Proj/MatSciCuration_TBV/DataCollection/PolymerDataCollection/CROW/family_links.csv",'w') as f:
    f.writelines(newlinks)



# use all the polymer types as a method of getting links to each polymer 


with open("/home/matt/Proj/MatSciCuration_TBV/DataCollection/PolymerDataCollection/CROW/family_links.csv") as f:
    lines = f.readlines()


# list to hold links to specific polymers
polymers = []


for line in lines:
    line = str(line).strip()
    #print (line)
    
    link ="https://www.polymerdatabase.com/polymer%20index/"+line



    
    page = requests.get(link)


    

    if page.ok == True:
        

        html = page.text

        soup = BeautifulSoup(html, "html.parser")

        

        links = soup.find_all("a")

        for link in links:
            href = link.get("href")
            if href.startswith("../polymers"):
                polymers.append(href)
            

    with open("polymer_links.csv", "w", newline="", encoding="utf-8") as csvfile:
        csv_writer = csv.writer(csvfile)
        

        for poly_link in polymers:
            csv_writer.writerow([poly_link])



