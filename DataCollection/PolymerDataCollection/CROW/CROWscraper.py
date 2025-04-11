from bs4 import BeautifulSoup
import json
import requests



# this code will scrap polymer data from the CROW database based off links 

url = 'https://www.polymerdatabase.com/polymers/polyacrylamide.html'
page = requests.get(url)
soup = BeautifulSoup(page.text, 'html.parser')


"""

# Extract the required information

def polyextract(soup):
    polymer_name = soup.find("b").text.strip()

    sections = {}
    for div in soup.find_all("div", class_="datagrid"):
        section_title = div.find("h3").text.strip()
        table = div.find("table")
        
        section_data = []
        for row in table.find_all("tr"):
            cells = [cell.text.strip() for cell in row.find_all("td")]
            row_data = {}
            for i in range(0, len(cells), 2):
                row_data[cells[i]] = cells[i + 1] if i + 1 < len(cells) else None
            section_data.append(row_data)
        
        sections[section_title] = section_data

    data = {
        'polymer name': polymer_name,
        **sections
    }
    return data




# Save the data to a JSON file
with open("scraped_data.json", "w", encoding="utf-8") as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

print("Data saved to scraped_data.json")

"""





def polyextract(soup):
    polymer_name = soup.find("b").text.strip()

    sections = {}
    for div in soup.find_all("div", class_="datagrid"):
        section_title = div.find("h3").text.strip()
        table = div.find("table")
        
        section_data = []
        for row in table.find_all("tr"):
            cells = [cell.text.strip() for cell in row.find_all("td")]
            row_data = {}
            for i in range(0, len(cells), 2):
                row_data[cells[i]] = cells[i + 1] if i + 1 < len(cells) else None
            section_data.append(row_data)
        
        sections[section_title] = section_data

    data = {
        'polymer name': polymer_name,
        **sections
    }
    return data



# Scaping on all from link collection 

with open("/home/matt/Proj/MatSciCuration_TBV/DataCollection/PolymerDataCollection/CROW/polymer_links.csv") as f:
    lines = f.readlines()



Data = []

for line in lines:
    line = line.replace("../","")
    line = str(line).strip()
    #print (line)
    
    link ="https://www.polymerdatabase.com/"+line


    print(link)
    

    page = requests.get(link)


    if page.ok == True:
        

        html = page.text

        soup = BeautifulSoup(html, "html.parser")

        data = polyextract(soup)

        

        Data.append(data)

with open("scraped_data.json", "w", encoding="utf-8") as f:
    json.dump(Data, f, ensure_ascii=False, indent=4)

print("Data saved to scraped_data.json") 





