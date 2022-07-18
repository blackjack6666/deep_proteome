"""
some function/features for MatrisomeDB2.0
"""


from urllib.request import Request, urlopen
import time
import random
import pandas as pd
from datetime import datetime
from bs4 import BeautifulSoup
import requests
import json
from collections import defaultdict
from multiprocessing_naive_algorithym import *
from aho_corasick import automaton_trie,automaton_matching
from protein_coverage import fasta_reader
from tsv_reader import psm_reader


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
            # print (each, split_dict[each]['n'], split_dict[each]['st'], split_dict[each]['en'])
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


def peptide_map(psm_dict,protein_dict):
    """
    map peptides to proteome and return freq array dict
    :param psm_dict: {'peptide': frequency(int)}
    :return:
    """
    id_list,seq_list = extract_UNID_and_seq(protein_dict)
    seqline = creat_total_seq_line(seq_list,sep='|')
    zeroline = zero_line_for_seq(seqline)
    sep_pos_array = separator_pos(seqline)

    prot_freq_dict = {}
    aho_result = automaton_matching(automaton_trie([pep for pep in psm_dict.keys()]),seqline)
    for pos in aho_result:

        zeroline[pos[0]:pos[1] + 1] += psm_dict[pos[2]]  # map PSMs instead of peptides

    for i in range(len(sep_pos_array) - 1):  # iterate from the left of zeroline all the way to the right
        prot_freq_dict[id_list[i]] = zeroline[sep_pos_array[i] + 1:sep_pos_array[i + 1]]

    return prot_freq_dict


def plot_domain_coverage(prot_freq_dict,domain_pos_dict, protein_entry:str):
    """
    plot the domain coverage, sum of spec counts for all domain and ranking
    :param prot_freq_dict: return by peptide_map
    :param domain_pos_dict: return by get_smart_info
    :param protein_entry: target protein to plot
    :return:
    """
    from bokeh.models import HoverTool, ColumnDataSource
    from bokeh.palettes import Spectral7
    from bokeh.transform import factor_cmap
    from bokeh.plotting import figure
    from bokeh.io import save, output_file, show

    freq_array = prot_freq_dict[protein_entry]
    domain_dict = domain_pos_dict[protein_entry]
    df = pd.DataFrame(columns=['pos_start_end','sum_spec_count','domain_name'])
    idx = 0
    for each_domain in domain_dict:
        for each_tp in domain_dict[each_domain]:
            start, end = each_tp[0], each_tp[1]
            sum_spec_count = np.sum(freq_array[start-1:end])/(end-start)  # spec count sum for each domain entry, average by length
            start_end = str(start)+'_'+str(end)
            df.at[idx,'pos_start_end'] = start_end
            df.at[idx,'sum_spec_count'] = sum_spec_count
            df.at[idx,'domain_name'] = each_domain
            idx += 1

    domain_name = sorted(df.domain_name.unique())
    print (domain_name)
    domain_sum = df.groupby('domain_name')['sum_spec_count'].sum()
    print (domain_sum)
    source = ColumnDataSource(data=dict(domain=domain_name, counts=domain_sum, color=Spectral7))
    #
    # cmap = factor_cmap('domain_name', palette=Spectral7, factors=sorted(df.domain_name.unique()))
    p = figure(x_range=domain_name, plot_height=350, toolbar_location='right',
               y_axis_label="normalized total spec count",
               title=f'{protein_entry} domain coverage')

    p.vbar(x='domain', top='counts', width=0.9, line_color="white",
           color='color', legend="domain", source=source)
    p.xgrid.grid_line_color = None
    p.legend.orientation = "horizontal"
    p.legend.location = "top_center"
    # p.vbar(x=df.index,
    #        top=df.sum_spec_count,
    #        fill_alpha = .5,
    #        fill_color = 'salmon',
    #        line_alpha = .5,
    #        line_color='green',
    #        line_dash='dashed')
    # # save(p, "test_bokeh.html", title="Test")
    show(p)
    return df


def plot_domain_coverage2(prot_freq_dict,domain_pos_dict, protein_entry:str):
    """
    same as plot_domain_coverage, but group by domain name
    :param prot_freq_dict:
    :param domain_pos_dict:
    :param protein_entry:
    :return:
    """
    from bokeh.models import HoverTool, ColumnDataSource, FactorRange
    from bokeh.palettes import Spectral7
    from bokeh.transform import factor_cmap
    from bokeh.plotting import figure
    from bokeh.io import save, output_file, show
    from math import pi


    freq_array = prot_freq_dict[protein_entry]
    domain_dict = domain_pos_dict[protein_entry]

    # a dataframe to store data
    df = pd.DataFrame(columns=['pos_start_end', 'sum_spec_count', 'domain_name'])
    idx = 0
    for each_domain in domain_dict:
        for each_tp in domain_dict[each_domain]:
            start, end = each_tp[0], each_tp[1]
            sum_spec_count = np.sum(freq_array[start - 1:end]) / (
                        end - start)  # spec count sum for each domain entry, average by length
            start_end = str(start) + '_' + str(end)
            df.at[idx, 'pos_start_end'] = start_end
            df.at[idx, 'sum_spec_count'] = sum_spec_count
            df.at[idx, 'domain_name'] = each_domain
            idx += 1

    # group = df.groupby(('domain_name'))
    # source = ColumnDataSource(group)
    # print (','.join(source.column_names))

    # plotting with bokeh, embed into html later
    # color map for each domain, hex color randomly generated
    color_map = {each:color_generator() for each in df.domain_name.unique()}

    factors = [(domain, pos)for domain, pos in zip(df.domain_name, df.pos_start_end)]  # x axis
    y = df.sum_spec_count  # y axis

    p = figure(x_range=FactorRange(*factors), plot_height=400, plot_width=len(factors)*40,
               y_axis_label="normalized total spec count",
               title=f'{protein_entry} domain coverage')
    p.vbar(x=factors, top=y, width=0.5, alpha=0.5,color=[color_map[tp[0]] for tp in factors])
    p.y_range.start = 0
    p.x_range.range_padding = 0.05
    p.xgrid.grid_line_color = None
    # p.xaxis.axis_label = 'whatever'
    p.xaxis.major_label_text_font_size = "10pt"
    p.xaxis.group_text_font_size = "15pt"
    # p.xaxis.major_label_text_font = 'Arial'
    p.xaxis.major_label_orientation = pi/2
    show(p)


def color_generator():
    import random
    r = lambda: random.randint(0,255)
    return '#%02X%02X%02X' % (r(),r(),r())

if __name__=='__main__':
    prot_list = ['Q8TER0','E9PWQ3','P11276']
    info_dict = get_smart_info(prot_list)
    psm_tsv = 'D:/data/Naba_deep_matrisome/07232021_secondsearch/SNED1_1080D/psm.tsv'
    psm_dict = psm_reader(psm_tsv)[0]
    protein_dict = fasta_reader('D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta')
    protein_freq_dict = peptide_map(psm_dict,protein_dict)

    plot_domain_coverage2(protein_freq_dict,info_dict,'P11276')

    # main('https://smart.embl.de/smart/show_motifs.pl?ID=Q8TER0')