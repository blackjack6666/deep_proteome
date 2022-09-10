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
from tsv_reader import psm_reader, protein_info_from_fasta
from bokeh.models import HoverTool, ColumnDataSource, FactorRange, LinearColorMapper,ColorBar,BasicTicker,PrintfTickFormatter, Plot, Rect, Legend, LegendItem,SingleIntervalTicker, Label, LabelSet, TableColumn, DataTable, HTMLTemplateFormatter
from bokeh.palettes import Spectral7, Viridis, Plasma, Blues9, Turbo256
from bokeh.transform import factor_cmap
from bokeh.plotting import figure
from bokeh.io import save, output_file, show
from bokeh.embed import components
from math import pi
from subprocess import call

# set pythonhashseed to static to make sure hashfunciton produce same hash value

ptm_map_dict = {'Q\\[129\\]':'Gln deamidation','N\\[115\\]':'ASN deamidation',
                'Q\\[111\\]':'Gln to pyroGln','C\\[143\\]':'selenocysteine',
                'M\\[147\\]':'Met oxidation','P\\[113\\]':'Pro hydroxylation',
                'K\\[144\\]':'Lys hydroxylation'}

def get_smart_info(protein_list:list):
    """
    -----
    crawl smart db to get domain positions
    -----
    :param protein_list: a list of protein uniprot IDs
    :return:
    """
    time_start = time.time()
    info_dict = {}
    for prot in protein_list:
        print (f'crawling {prot} from SMART...')
        domain_dict = defaultdict(list)
        time.sleep(0.1)
        req = Request('https://smart.embl.de/smart/show_motifs.pl?ID='+prot,
                      headers={'User-Agent': 'Mozilla/5.0'})
        webpage = urlopen(req).read().decode('utf-8')
        try:
            web_split = webpage.split('domNfo=')[1].split('};')[0]+'}'  # scrap domain info.
            split_dict = json.loads(web_split) # convert dict string into dict structure
            # print (split_dict)
            for each in split_dict:
                # print (each, split_dict[each]['n'], split_dict[each]['st'], split_dict[each]['en'])
                domain_dict[split_dict[each]['n']].append((int(split_dict[each]['st']),int(split_dict[each]['en'])))
            info_dict[prot] = domain_dict
        except IndexError:
            print (f'{prot} has no domain info from SMART')
            continue
    print (f'crawler takes {time.time()-time_start}')
    return info_dict


def smart_scrap(protein_list:list):
    """
    separate confidently predicted domains and overlapped domains
    :param protein_list:
    :return:
    """

    time_start = time.time()
    conf_predict_info_dict, overlapped_info_dict = {}, {}
    for prot in protein_list:
        time.sleep(1)
        conf_domain_dict, overlap_domain_dict = defaultdict(set),defaultdict(set)
        req = Request('https://smart.embl.de/smart/show_motifs.pl?ID=' + prot,
                      headers={'User-Agent': 'Mozilla/5.0'})
        webpage = urlopen(req).read().decode('utf-8')
        conf_web_split = webpage.split('var visible_DS = new kendo.data.DataSource(')[1].split(']')[0]

        for start, end, domain in zip(re.findall('\"st\"\:\"\d+', conf_web_split),
                                      re.findall('\"en\"\:\"\d+', conf_web_split),
                                      re.findall('\>\w+\s*\w+\<', conf_web_split)):
            conf_domain_dict[domain[1:-1]].add((int(start.split('"')[-1]), int(end.split('"')[-1])))

        overlapped_web_split = webpage.split('var hidden_DS = new kendo.data.DataSource(')[1].split(']')[0]
        for start, end, domain in zip(re.findall('\"st\"\:\"\d+', overlapped_web_split),
                                      re.findall('\"en\"\:\"\d+', overlapped_web_split),
                                      re.findall('\>\w+\s*\w+\<', overlapped_web_split)):
            overlap_domain_dict[domain[1:-1]].add((int(start.split('"')[-1]),int(end.split('"')[-1])))
        conf_predict_info_dict[prot] = conf_domain_dict
        overlapped_info_dict[prot] = overlap_domain_dict
        print (f'{prot} crawl takes {time.time()-time_start}')
    return conf_predict_info_dict,overlapped_info_dict


def my_replace(match_obj):
    match_obj = match_obj.group()
    matched_aa = match_obj[0]
    if matched_aa != 'n':
        return matched_aa  # gives back the first element of matched object as string
    else:
        # if first match is n, then n acetylation, get rid of n
        print (match_obj)
        return ''


def peptide_map(psm_dict,protein_dict):
    """
    map peptides to proteome and return freq array dict
    :param psm_dict: {'peptide': frequency(int)}
    :return:
    """
    start = time.time()
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
    print (f'peptide mapping took {time.time()-start}s')

    return prot_freq_dict


def ptm_map(psm_list,protein_dict):
    """
    map psms with ptm to protein sequence
    :param psm_list: a list of PSMs with PTM, P[100]EPTIDE
    :param protein_dict:
    :return:
    """
    time_start = time.time()
    regex_pat = '\w{1}\[\d+\.?\d+\]'  # universal ptm pattern

    regex_set = set()
    peptide_psm_dict = defaultdict(list)  # append all psm into a dictionary, {peptide:[psm1,psm2,...]}
    id_ptm_idx_dict = {}  # {protein_id:{ptm1:nonzero_index_array,ptm2:nonzero_index_array,...}}
    id_ptm_freq_dict = {}  # {protein_id:{ptm1:freq_array,ptm2:freq_array,...}}

    for each in psm_list:
        each_reg_sub = re.sub(regex_pat, my_replace, each)
        peptide_psm_dict[each_reg_sub].append(each)
        match = re.findall(regex_pat, each)
        if match:
            for ptm in match:
                regex_set.add(ptm.replace('[','\[').replace(']','\]').replace('.','\.'))
    print (regex_set)
    print (peptide_psm_dict['MSEAAPAAPAAAPPAEK'])
    # aho mapping
    id_list, seq_list = extract_UNID_and_seq(protein_dict)
    seq_line = creat_total_seq_line(seq_list, sep="|")
    ptm_index_line_dict = {each:zero_line_for_seq(seq_line) for each in regex_set}
    separtor_pos_array = separator_pos(seq_line)
    aho_result = automaton_matching(automaton_trie([pep for pep in peptide_psm_dict]), seq_line)

    # ptm mapping
    for tp in aho_result:
        matched_pep = tp[2]  # without ptm site
        for psm in peptide_psm_dict[matched_pep]:
            for ptm in regex_set:  # check each ptm, mask other ptms
                new_psm = re.sub('n?\[\d+\.?\d+\]', '', psm.replace(ptm.replace('\\', ''), '*')).\
                    replace('*', ptm.replace('\\', ''))
                ptm_mod = set(re.findall(ptm, new_psm))
                if ptm_mod:
                    for ele in ptm_mod:
                        ### count multiple ptms in a peptide seq
                        num_of_mod = len(
                            re.findall(ele.replace('[', '\[').replace(']', '\]').replace('.', '\.'), new_psm))

                        PTM_index = [m.start() for m in
                                     re.finditer(ele.replace('[', '\[').replace(']', '\]').replace('.', '\.'), new_psm)]
                        PTM_index_clean = [ind - num * (len(ele) - 1) for ind, num in zip(PTM_index, range(num_of_mod))]
                        for indx in PTM_index_clean:
                            ptm_index_line_dict[ptm][tp[0] + indx] += 1

    # get ptm freq array and index
    for i in range(len(separtor_pos_array) - 1):
        id_ptm_idx_dict[id_list[i]] = {ptm:
            np.nonzero(ptm_index_line_dict[ptm][separtor_pos_array[i] + 1:separtor_pos_array[i + 1]])[0]
            for ptm in regex_set}
        id_ptm_freq_dict[id_list[i]] = {ptm:
            ptm_index_line_dict[ptm][separtor_pos_array[i] + 1:separtor_pos_array[i + 1]]
            for ptm in regex_set}
    print (f'PTM mapping took {time.time()-time_start}s')
    return id_ptm_idx_dict, id_ptm_freq_dict


def plot_domain_coverage(prot_freq_dict,domain_pos_dict, protein_entry:str):
    """
    plot the domain coverage, sum of spec counts for all domain and ranking
    :param prot_freq_dict: return by peptide_map
    :param domain_pos_dict: return by get_smart_info
    :param protein_entry: target protein to plot
    :return:
    """

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


def plot_domain_coverage2(prot_freq_dict,domain_pos_dict, protein_entry:str, mode='sum'):
    """
    same as plot_domain_coverage, but group by domain name
    :param prot_freq_dict:
    :param domain_pos_dict:
    :param protein_info_dict:
    :param protein_entry:
    :param mode: plot mode, default is sum, could be 'coverage'
    :return:
    """
    time_start = time.time()
    freq_array = prot_freq_dict[protein_entry]
    domain_dict = domain_pos_dict[protein_entry]

    # a dataframe to store data, might not be needed
    df = pd.DataFrame(columns=['pos_start_end', 'sum_spec_count', 'domain_name','coverage'])
    idx = 0
    for each_domain in domain_dict:
        for each_tp in domain_dict[each_domain]:
            start, end = each_tp[0], each_tp[1]
            # spec count sum for each domain entry, average by length, or coverage in 100
            sum_spec_count = np.sum(freq_array[start - 1:end]) / (
                        end - start)
            coverage = np.count_nonzero(freq_array[start - 1:end])/(end-start+1)
            start_end = str(start) + '_' + str(end)
            df.at[idx, 'pos_start_end'] = start_end
            df.at[idx, 'sum_spec_count'] = sum_spec_count
            df.at[idx, 'domain_name'] = each_domain
            df.at[idx, 'coverage'] = coverage
            idx += 1

    # group = df.groupby(('domain_name'))
    # source = ColumnDataSource(group)
    # print (','.join(source.column_names))

    # plotting with bokeh, embed into html later
    # color map for each domain, hex color randomly generated
    color_map = {each:color_generator() for each in df.domain_name.unique()}

    factors = [(domain, pos)for domain, pos in zip(df.domain_name, df.pos_start_end)]  # x axis

    y = df.sum_spec_count  # y axis
    colors = [color_map[tp[0]] for tp in factors] # color legend
    source = ColumnDataSource(data=dict(x=factors,
                                        y=y,
                                        color=colors,
                                        pos=df.pos_start_end,
                                        domain=df.domain_name,
                                        cov=df.coverage))

    # hover tools
    TOOLTIPS = [
        ("domain", "@domain"),
        ("start end position", "@pos"),
        ("normalized spec count", "@y"),
        ("sequence coverage", "@cov")
    ]
    y_axis_label = "domain coverage"

    p = figure(x_range=FactorRange(*factors), plot_height=400, plot_width=len(factors)*35, tooltips=TOOLTIPS,
               y_axis_label=y_axis_label,
               title='Domain coverage')
    # bar chart and line plot on the same figure
    p.vbar(x='x',top='y',color='color',source=source,legend_label='normalized spec count')
    p.line(x='x',y='cov',source=source, line_color='red',line_dash='dotdash',line_dash_offset=1,legend_label='sequence coverage',
           line_width=2)
    p.y_range.start = 0
    p.x_range.range_padding = 0.05
    p.xgrid.grid_line_color = None
    # p.xaxis.axis_label = 'whatever'
    p.xaxis.major_label_text_font_size = "8pt"
    p.xaxis.group_text_font_size = "10pt"
    # p.xaxis.major_label_text_font = 'Arial'
    p.xaxis.major_label_orientation = pi/2
    # show(p)
    print (f'time used {time.time()-time_start}s')
    # return the js script and div for html embedding
    return components(p)


def domain_cov_ptm(prot_freq_dict, ptm_map_result, domain_pos_dict,protein_entry:str):
    """
    -----
    draw rectangular box as protein domains and alpha as coverage on bokeh,
    and label PTMs.
    -----
    :param protein_freq_dict:
    :param domain_pos_dict:
    :param protein_entry:
    :return:
    """
    time_start = time.time()
    freq_array = prot_freq_dict[protein_entry]
    domain_dict = domain_pos_dict[protein_entry]
    protein_len = len(freq_array)

    ## coverage every 10 aa
    pos_cov_dict = {}
    bin_width = 10
    bar_shrink_raio = 5 # make bar shorter
    bar_bottom = 0.8
    interval=np.arange(0,protein_len,bin_width)
    for i in interval[:-1]:
        coverage = np.count_nonzero(freq_array[i:i+bin_width])/bin_width
        pos_cov_dict[i+bin_width/2] = coverage/bar_shrink_raio + bar_bottom  # move bar up
    pos_cov_dict[interval[-1]+bin_width/2] = np.count_nonzero(freq_array[interval[-1]:protein_len])/(protein_len-interval[-1])/bar_shrink_raio+bar_bottom
    source_cov = ColumnDataSource(dict(x=[k for k in pos_cov_dict.keys()],y=[v for v in pos_cov_dict.values()],
                                       label=['{:.1f}%'.format((v-0.8)*100) for v in pos_cov_dict.values()]))

    ## extract domain position and calculate domain coverage
    info_list = []
    for each_domain in domain_dict:
        for each_tp in domain_dict[each_domain]:
            start, end = each_tp[0], each_tp[1]
            coverage = np.count_nonzero(freq_array[start - 1:end]) / (end - start + 1)
            info_list.append((start,end,coverage,each_domain))

    start, end, coverage, domain_list = zip(*info_list)

    # x coordinates and widths of rectangular
    x,width = zip(*[((end_-start_)/2+start_,end_-start_) for end_,start_ in zip(end,start)])
    # hash color to show each domain
    color_map_dict = {domain:hashcolor(domain) for domain in set(domain_list)}
    color_list = [color_map_dict[each] for each in domain_list]

    source = ColumnDataSource(dict(x=x,w=width,color=color_list,domain=domain_list,
                                   start=start,end=end,
                                   coverage=coverage))

    ## prepare data for PTM labeling
    id_ptm_indx_dict = ptm_map_result[0]
    ptm_index_dict = id_ptm_indx_dict[protein_entry]
    ptm_index_sort = sorted([(idx,each) for each in ptm_index_dict for idx in ptm_index_dict[each]])

    # bokeh plot, hovertool
    hover = HoverTool(names=['rec'],tooltips=[('domain', '@domain'), ('start position', '@start'),('end position','@end'),('domain coverage','@coverage{:.1%}'),])
    # initiate bokeh figure
    p = figure(x_range=(-10,protein_len),
               y_range=(0,2),
               tools=['pan', 'box_zoom', 'wheel_zoom', 'save',
                      'reset', hover],
               plot_height=500, plot_width=1600,
               toolbar_location='right',
               title='',
               x_axis_label='amino acid position')

    # plot domains as rectangular and alpha shows coverage
    p.rect(x="x", y=0.6, width='w', height=50,
           source=source,
           fill_color='color',
    #       fill_alpha='coverage',
           line_width=2,
           line_color='black',
           height_units="screen",
           name='rec'
           )
    # reference domains, alpha=1
    # p.rect(x="x", y=1, width='w', height=10,
    #        source=source,
    #        fill_color='color',
    #        line_width=2,
    #        line_color='black',
    #        height_units="screen",
    #        name='rec'
    #        )

    # sequence coverage bar charts
    p.vbar(x='x',width=bin_width,top='y',bottom=bar_bottom,source=source_cov,color='#D3D3D3', name='seq cov')
    # cov_label = LabelSet(x='x',y='y',text='label',text_font_size='8px',
    #                      x_offset=-13.5, y_offset=0, source=source_cov)
    # p.add_layout(cov_label)
    cov_bar_legend_top, cov_bar_legend_bottom = bar_bottom+1/bar_shrink_raio,bar_bottom
    cov_bar_x_coor = -8
    # p.vbar(x=[cov_bar_x_coor],width=bin_width,top=[cov_bar_legend_top],bottom=[cov_bar_legend_bottom],
    #        color='#D3D3D3',name='seq_cov_legend')
    # label annotations
    for y_coor, text in zip([cov_bar_legend_bottom,cov_bar_legend_top],['0%','100%']):
        label_cov = Label(x=cov_bar_x_coor,y=y_coor,x_offset=0, y_offset=-5,text=text,text_font_size='10px',text_align='left')
        p.add_layout(label_cov)
    seq_cov_title = Label(x=-5, y=1.1, text='Sequence coverage binned by every 10 aa',text_font_size='12px',text_align='left',text_color='#A9A9A9')
    p.add_layout(seq_cov_title)

    # line shows whole protein length
    p.line(x=[0,protein_len],y=[0.6,0.6],line_width=10,color='#000000',alpha=0.8,name='line')

    # adjusted PTM text coordinates calculation
    numpy_zero_array = np.zeros((200, protein_len)) # mask numpy array for text plotting
    ptm_x,ptm_y = [], []
    new_ptm_x, new_ptm_y = [], [] # adjusted text coordinates to prevent overlap
    ptms = []
    for tp in ptm_index_sort:
        each_idx, ptm = tp
        # ptms.append(ptm.replace('\\', ''))
        ptms.append(ptm)
        ptm_x.append(each_idx)
        ptm_y.append(0.5)
        x_offset,y_offset = 0,0
        while True: # keep moving down if text are too close
            nonzero_count = np.count_nonzero(numpy_zero_array[190+y_offset:200+y_offset,each_idx+x_offset:each_idx+130+x_offset])
            if nonzero_count == 0:
                # print (ptm,each_idx,x_offset,y_offset)
                new_ptm_x.append(each_idx+x_offset)
                new_ptm_y.append((25+y_offset)/200*2)
                numpy_zero_array[190+y_offset:200+y_offset,each_idx+x_offset:each_idx+130+x_offset] += 1
                break
            else:
                # print ('moving down')
                # x_offset += 50 # value to move right
                y_offset -= 12 # value to move down

    # label ptm and connect to protein domains
    for x,y,x_,y_,ptm in zip(new_ptm_x,new_ptm_y,ptm_x,ptm_y,ptms):
        p.line(x=[x_+1,x+1],y=[y_,y+0.1],line_width=1,color='black',alpha=0.3) # connect domain with text
        label = Label(x=x,y=y,text=ptm_map_dict[ptm]+'\n'+str(x_+1),text_font_size='10px', text_align='center', text_font='Tahoma')
        p.add_layout(label)
    # dummy glyphs to help draw legend
    legend_gly = [p.line(x=[1, 1], y=[1, 1], line_width=15, color=c, name='dummy_for_legend')
                    for c in [v for v in color_map_dict.values()]]

    legend = Legend(title='Domains', background_fill_color='white',
                    border_line_color='black',border_line_width=3,
                    border_line_alpha=0.7,
                    items=[LegendItem(label=lab, renderers=[gly])
                           for lab, gly in zip([d for d in color_map_dict.keys()],legend_gly)])
    # alpha color bar, domain coverage
    # color_mapper = LinearColorMapper(palette=Blues9[::-1], low=0, high=1)
    # color_bar = ColorBar(color_mapper=color_mapper,location=(0, 0),ticker=SingleIntervalTicker(interval=0.1))
    # p.add_layout(color_bar,'right')

    p.add_layout(legend)
    p.xgrid.visible = False
    p.ygrid.visible = False
    p.yaxis.visible = False
    print (f'bokeh graph took {time.time()-time_start}s')
    # show(p)
    return components(p)


def ptm_domain_htmap(ptm_map_result, domain_pos_dict, protein_entry:str):
    """
    use bokeh to generate a heatmap showing ptm coverage by domains
    :param ptm_map_result: return from ptm_map, two returned values
    :param domain_pos_dict: return from get_smart_info, SMART web crawler
    :param protein_info_dict: {uniprotid:(gene,description,cls)}
    :param protein_entry: target Uniprot ID
    :return:
    """
    time_start = time.time()
    id_ptm_indx_dict, id_ptm_freq_array_dict = ptm_map_result
    ptm_indx_dict, ptm_freq_array_dict = id_ptm_indx_dict[protein_entry],id_ptm_freq_array_dict[protein_entry]
    domain_dict = domain_pos_dict[protein_entry]

    # create pd dataframe for plotting
    x = [domain+'_'+str(tp[0])+'_'+str(tp[1]) for domain in domain_dict for tp in domain_dict[domain]]
    y = [k.replace('\\','') for k in ptm_indx_dict.keys()]
    # print (y)
    df = pd.DataFrame(index=y,columns=x)
    df.columns.name = 'domain_position'
    df.index.name = 'PTMs'
    for domain in domain_dict:
        for tp in domain_dict[domain]:
            start,end = tp[0], tp[1]
            for each_ptm in ptm_freq_array_dict:
                ptm_occur = np.count_nonzero(ptm_freq_array_dict[each_ptm][start-1:end])
                df.at[each_ptm.replace('\\',''),domain+'_'+str(start)+'_'+str(end)] = ptm_occur

    df_stack = pd.DataFrame(df.stack(),columns=['value']).reset_index() # domain and PTMs for each row

    # plotting
    # colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce", "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]
    mapper = LinearColorMapper(palette=Plasma[5], low=df_stack.value.min(), high=df_stack.value.max())

    p = figure(x_range=x, y_range=y, plot_height=400, plot_width=len(x)*30,
               tooltips=[('domain position and PTM', '@domain_position @PTMs'), ('PTM occurrence', '@value')],
               y_axis_label="PTMs",
               title='PTM heatmap')

    p.rect(x="domain_position", y="PTMs", width=1, height=1,
           source=df_stack,
           fill_color={'field': 'value', 'transform': mapper},
           line_color=None)

    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "10px"
    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = pi / 3

    color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="7px",
                         ticker=BasicTicker(desired_num_ticks=len(Plasma[5])),
                         formatter=PrintfTickFormatter(format="%d"),
                         label_standoff=6, border_line_color=None)
    p.add_layout(color_bar, 'right')

    # show(p)
    print (f'time used {time.time()-time_start}s')
    return components(p)


def ptm_table_bokeh(id_ptm_idx_dict, protein_dict, protein_info_dict):
    """
    plot a bokeh table to show ptms frequency
    :param id_ptm_idx_dict:
    :param protein_dict:
    :param uniprot_gene_dict:
    :return:
    """

    info_dict = defaultdict(list)
    for prot in id_ptm_idx_dict:
        info_dict['Uniprot ID'].append(prot)
        info_dict['Gene'].append(protein_info_dict[prot][0])
        info_dict['Length'].append(len(protein_dict[prot]))
        info_dict['url'].append('%s.html' % prot)  # clickable protein id to other page
        for ptm in id_ptm_idx_dict[prot]:
            info_dict[ptm_map_dict[ptm]].append(len(id_ptm_idx_dict[prot][ptm]))
    df = pd.DataFrame(info_dict)
    # print (df)
    source = ColumnDataSource(df)

    columns = []
    columns.append(TableColumn(field='Uniprot ID',title='Uniprot ID',
                               formatter=HTMLTemplateFormatter(template='<a href="<%= url %>"><%= value %></a>')))
    columns += [TableColumn(field=each,title=each)
                for each in df.columns[1:] if each!='url']
    table = DataTable(source=source,columns=columns, width=1000, height=600, editable=True)
    # show(table)

    return components(table)


def color_generator():
    import random
    r = lambda: random.randint(0,255)
    return '#%02X%02X%02X' % (r(),r(),r())


def combine_bokeh(domain_bokeh_return, ptm_bokeh_return, protein_info_dict, html_out='test.html', UniportID=''):
    """
    generate the domain coverage bar graph and ptm domain heatmap in the same html
    :param domain_bokeh_return: script, div from plot_domain_coverage2
    :param ptm_bokeh_return: script, div from ptm_domain_htmap
    :param protein_info_dict: {uniprotid:(gene,description,cls)}
    :return:
    """
    # load bokeh js scripts and divs
    coverage_script, coverage_div = domain_bokeh_return
    ptm_script, ptm_div = ptm_bokeh_return
    smart_url = 'https://smart.embl.de/smart/show_motifs.pl?ID='+UniportID

    # print (coverage_script,coverage_div)
    # read html template
    html_template_f = open('F:/matrisomedb2.0/bokeh_html_template.html')
    html_template = html_template_f.read()
    html_template_f.close()

    # write new html with js script and divs
    new_html = html_template.replace('<!-- COPY/PASTE domain coverage SCRIPT HERE -->',coverage_script).\
        replace('<!-- COPY/PASTE ptm SCRIPT HERE -->',ptm_script).replace('<!-- INSERT domain DIVS HERE -->',coverage_div)\
        .replace('<!-- INSERT ptm DIVS HERE -->',ptm_div)

    new_html = new_html.replace('<!-- UniprotID -->',protein_info_dict[UniportID][0]+'  ('+protein_info_dict[UniportID][1]+')').replace('<!-- SMART URL -->',smart_url)
    with open(html_out, 'w',newline='\n') as f_w:
        f_w.write(new_html)

    return new_html


def bokeh_to_html(domain_cov_ptm_bokeh, ptm_table_bokeh, protein_info_dict, html_out='test.html',UniprotID=''):
    # load bokeh js scripts and divs
    time_start = time.time()
    js_script, div = domain_cov_ptm_bokeh
    ptm_js_script, ptm_div = ptm_table_bokeh
    smart_url = 'https://smart.embl.de/smart/show_motifs.pl?ID=' + UniprotID

    # read html template
    html_template_f = open('F:/matrisomedb2.0/bokeh_html_template.html')
    html_template = html_template_f.read()
    html_template_f.close()

    # write new html with js script and divs
    # domain graph div
    new_html = html_template.replace('<!-- COPY/PASTE domain coverage SCRIPT HERE -->', js_script).\
        replace('<!-- INSERT domain DIVS HERE -->',div)
    # ptm table div
    new_html = new_html.replace('<!-- COPY/PASTE ptm SCRIPT HERE -->',ptm_js_script).\
        replace('<!-- INSERT ptm DIVS HERE -->',ptm_div)

    new_html = new_html.replace('<!-- UniprotID -->',
                                protein_info_dict[UniprotID][0] + '  (' + protein_info_dict[UniprotID][
                                    1].rstrip(' ') + ')').replace('<!-- SMART URL -->', smart_url)
    with open(html_out, 'w', newline='\n') as f_w:
        f_w.write(new_html)
    print (f'bokeh to html {time.time()-time_start}s')
    return new_html


def hashcolor(s):

    return Turbo256[hash(s) % 256]


if __name__ == '__main__':
    from tsv_reader import modified_peptide_from_psm
    import json

    # SMART web crawler to extract domain info
    # prot_list = ['P11276']
    # info_dict = get_smart_info(prot_list)

    # read domain from json file
    with open('F:/matrisomedb2.0/smart_domain.json') as f_o:
        info_dict = json.load(f_o)

    # protein_info
    protein_info_dict = protein_info_from_fasta('D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta')

    # use subset of protein dict
    protein_dict = fasta_reader('D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta')
    ecm_df = pd.read_excel('F:/matrisomedb2.0/AnnotatedDeDuplicated.xlsx', index_col=0)
    ecm_list = ecm_df[ecm_df['Division'] != 'Non-matrisome'].index.tolist()
    ecm_protein_dict = {prot: protein_dict[prot] for prot in ecm_list}

    # peptide mapping
    psm_tsv = 'D:/data/Naba_deep_matrisome/07232021_secondsearch/SNED1_1080D/psm.tsv'
    psm_dict = psm_reader(psm_tsv)[0]
    protein_freq_dict = peptide_map(psm_dict,ecm_protein_dict)

    # ptm mapping
    # psm_list = modified_peptide_from_psm(psm_tsv)
    psm_list = ['n[42.01]M[100]LQTLFLTM[142]LTLALVK[0.984]']
    ptm_map_result = ptm_map(psm_list,ecm_protein_dict)[0]
    print (ptm_map_result['Q8BPB5'])
    # updated 8/3/22
    # domain_bokeh_return = domain_cov_ptm(protein_freq_dict,ptm_map_result, info_dict,protein_entry='E9PZ16')

    # update 8/22/22
    # ptm_bokeh_return = ptm_table_bokeh(ptm_map_result[0],ecm_protein_dict,protein_info_dict)

    # updated 8/3/22
    # bokeh_to_html(domain_bokeh_return,
    #               ptm_bokeh_return,
    #               protein_info_dict,
    #               html_out='F:/matrisomedb2.0/newbokeh_test_E9PZ16.html',
    #               UniprotID='E9PZ16')

    # domain coverage
    # domain_coverage_bokeh = plot_domain_coverage2(protein_freq_dict,info_dict,'E9PWQ3')

    # main('https://smart.embl.de/smart/show_motifs.pl?ID=Q8TER0')
    # ptm mapping
    # psm_list = modified_peptide_from_psm(psm_tsv)
    # ptm_map_result = ptm_map(psm_list,protein_dict)
    # ptm_domain_bokeh = ptm_domain_htmap(ptm_map_result,info_dict,'E9PWQ3')

    # combine domain coverage and ptm heatmaps in one html
    # combine_bokeh(domain_coverage_bokeh,
    #               ptm_domain_bokeh,
    #               protein_info_dict,
    #               html_out='F:/matrisomedb2.0/bokeh_test_E9PWQ3.html',
    #               UniportID='E9PWQ3')


    # info_dict = smart_scrap(['P20152'])
    # print (info_dict[0],'\n',info_dict[1])