import requests
import re
from bs4 import BeautifulSoup
def GetPicCAS(link):
    page= requests.get(link)
    soup= soup.find_all("img", {"class": "size-full"})
    img_link =  soup[0]['src']
    CAS = re.search('\d+-\d+-\d+',soup[0]['alt'])

    return CAS, img_link
    
