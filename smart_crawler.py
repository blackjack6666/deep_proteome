from urllib.request import Request, urlopen
import time
import random
from datetime import datetime
from bs4 import BeautifulSoup
import requests
import json
from collections import defaultdict


def get_smart_info(protein_list:list):
    """
    -----
    crawl smart db to get domain positions
    -----
    :param protein_list: a list of protein uniprot IDs
    :return:
    """
    info_dict = {}
    for prot in protein_list:
        domain_dict = defaultdict(set)
        time.sleep(1)
        req = Request('https://smart.embl.de/smart/show_motifs.pl?ID='+prot,
                      headers={'User-Agent': 'Mozilla/5.0'})
        webpage = urlopen(req).read().decode('utf-8')
        web_split = webpage.split('domNfo=')[1].split('};')[0]+'}'  # scrap domain info.
        split_dict = json.loads(web_split) # convert dict string into dict structure
        # print (split_dict)
        for each in split_dict:
            print (each, split_dict[each]['n'], split_dict[each]['st'], split_dict[each]['en'])
            domain_dict[split_dict[each]['n']].add((int(split_dict[each]['st']),int(split_dict[each]['en'])))
        info_dict[prot] = domain_dict
    return info_dict


def download_page(url):
    response = requests.get(url)
    response.raise_for_status()
    return response.text


def main(url):
    content = download_page(url)
    soup = BeautifulSoup(content, 'html.parser')
    result = {}
    table = soup.find('table', {"class":"k-focusable k-selectable"})
    print (table)
    for row in table.find_all('tr'):
        row_header = row.th.get_text()
        row_cell = row.td.get_text()
        result[row_header] = row_cell
    # with open('book_table.json', 'w') as storage_file:
    #     storage_file.write(json.dumps(result))

if __name__=='__main__':
    prot_list = ['Q8TER0']
    print (get_smart_info(prot_list))

    # main('https://smart.embl.de/smart/show_motifs.pl?ID=Q8TER0')