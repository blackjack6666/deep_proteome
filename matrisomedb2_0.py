from multiprocessing_naive_algorithym import *
from aho_corasick import automaton_trie, automaton_matching
from protein_coverage import fasta_reader
from tsv_reader import protein_info_from_fasta
from bokeh.models import HoverTool, ColumnDataSource, FactorRange, LinearColorMapper,ColorBar,BasicTicker,PrintfTickFormatter, Plot, Rect, Legend, LegendItem,SingleIntervalTicker, Label, LabelSet, TableColumn, DataTable, HTMLTemplateFormatter
from bokeh.palettes import Spectral7, Viridis, Plasma, Blues9, Turbo256
from bokeh.transform import factor_cmap
from bokeh.plotting import figure
from bokeh.io import save, output_file, show
from bokeh.embed import components

ptm_map_dict = {'Q\\[129\\]':'Gln deamidation','N\\[115\\]':'ASN deamidation',
                'Q\\[111\\]':'Gln to pyroGln','C\\[143\\]':'selenocysteine',
                'M\\[15\\.9949\\]':'Met oxidation','P\\[15\\.9949\\]':'Pro hydroxylation',
                'K\\[15\\.9949\\]':'Lys hydroxylation','n\\[42\\.0106\\]':'N-term acetylation',
                'C\\[57\\.0215\\]':'Cys redu-alky'}


def my_replace(match_obj):
    match_obj = match_obj.group()
    matched_aa = match_obj[0]
    if matched_aa != 'n':
        return matched_aa  # gives back the first element of matched object as string
    else:
        # if first match is n, then n acetylation, get rid of n
        return ''


def peptide_map(psm_dict,protein_dict):
    """
    map peptides to proteome and return freq array dict
    :param psm_dict: {'peptide': [psm1,psm2,...]}
    :return:
    """
    prot_freq_dict = {}
    coverage_dict = {}
    prot_psm_dict = defaultdict(int)
    prot_psm_list_dict = defaultdict(list)

    start = time.time()
    id_list,seq_list = extract_UNID_and_seq(protein_dict)
    seqline = creat_total_seq_line(seq_list,sep='|')
    zeroline = zero_line_for_seq(seqline)
    sep_pos_array = separator_pos(seqline)
    pos_id_dict = read_position_ID_into_dict(id_list,seq_list,zeroline)

    map_start = time.time()
    aho_result = automaton_matching(automaton_trie([pep for pep in psm_dict.keys()]),seqline)
    print (f'aho tree building and mapping took {time.time()-map_start}')
    for pos in aho_result:
        map_pep = pos[2]
        # zeroline[pos[0]:pos[1] + 1] += psm_dict[pos[2]]  # map PSMs instead of peptides
        zeroline[pos[0]:pos[1] + 1] += len(psm_dict[map_pep])
        prot_psm_dict[pos_id_dict[pos[0]]] += len(psm_dict[map_pep])
        prot_psm_list_dict[pos_id_dict[pos[0]]] += psm_dict[map_pep]
    for i in range(len(sep_pos_array) - 1):  # iterate from the left of zeroline all the way to the right
        freq_array = zeroline[sep_pos_array[i] + 1:sep_pos_array[i + 1]]

        prot_freq_dict[id_list[i]] = freq_array
        coverage_dict[id_list[i]] = np.count_nonzero(freq_array)/len(freq_array)*100

    print (f'script took {time.time()-start}s')

    return prot_freq_dict, prot_psm_dict, prot_psm_list_dict, coverage_dict


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
    if 'C\\[57\\.0215\\]' in regex_set:
        regex_set.remove('C\\[57\\.0215\\]')

    # print (regex_set)
    # aho mapping
    id_list, seq_list = extract_UNID_and_seq(protein_dict)
    seq_line = creat_total_seq_line(seq_list, sep="|")
    ptm_index_line_dict = {each:zero_line_for_seq(seq_line) for each in regex_set}
    separtor_pos_array = separator_pos(seq_line)
    aho_result = automaton_matching(automaton_trie([pep for pep in peptide_psm_dict]), seq_line)

    # ptm mapping
    ptm_start = time.time()
    for tp in aho_result:
        matched_pep = tp[2]  # without ptm site
        for psm in peptide_psm_dict[matched_pep]:
            for ptm in regex_set:  # check each ptm, mask other ptms
                new_psm = re.sub('\[\d+\.?\d+\]', '', psm.replace(ptm.replace('\\', ''), '*')).\
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
    print (f'ptm mapping took {time.time()-ptm_start}')

    # get ptm freq array and index
    for i in range(len(separtor_pos_array) - 1):
        id_ptm_idx_dict[id_list[i]] = {ptm:
            np.nonzero(ptm_index_line_dict[ptm][separtor_pos_array[i] + 1:separtor_pos_array[i + 1]])[0]
            for ptm in regex_set}
        id_ptm_freq_dict[id_list[i]] = {ptm:
            ptm_index_line_dict[ptm][separtor_pos_array[i] + 1:separtor_pos_array[i + 1]]
            for ptm in regex_set}
    print (f'ptm script took {time.time()-time_start}s')
    return id_ptm_idx_dict, id_ptm_freq_dict


def seq_cov_gen(freq_array,ptm_dict,protein_seq):
    """
    generate seq coverage map in html format (inside body)
    :param freq_array: mapped 1d np zero array for one protein
    :param protein_seq:
    :param protein_id:
    :param ptm_dict: {ptm1:[index list],ptm2:[index list]}
    :param gene_name:
    :return:
    """
    # marks for html
    """mark {
	background-color: #2C3E50;
	color: white;
}
mark1 {
	background-color: #2471A3;
	color: white;
}
mark2 {
	background-color: #AED6F1;
	color: white;
}
mark3 {
	background-color: #95A5A6;
	color: white;
}
mark4 {
	background-color: #D7DBDD;
	color: white;
}
mark5 {
	background-color: white;
	color: red;
}
    """

    seq_cov = np.count_nonzero(freq_array)/len(freq_array)*100
    split_seq = np.arange(0, len(protein_seq), 165)
    split_seq = np.append(split_seq,len(protein_seq))
    max_freq = np.max(freq_array)
    ptm_freq_array = set([idx for ptm in ptm_dict for idx in ptm_dict[ptm]])
    # print (ptm_freq_array)

    seq_cov_str = ''
    for i in range(len(split_seq) - 1):
        for j in range(split_seq[i], split_seq[i + 1]):
            if j not in ptm_freq_array:
                if freq_array[j] == 0:
                    seq_cov_str += protein_seq[j]
                elif 1 <= freq_array[j] < 0.2 * max_freq:
                    seq_cov_str += '<mark4>' + protein_seq[j] + '</mark4>'
                elif 0.2 * max_freq <= freq_array[j] < 0.4 * max_freq:
                    seq_cov_str += '<mark3>' + protein_seq[j] + '</mark3>'
                elif 0.4 * max_freq <= freq_array[j] < 0.6 * max_freq:  # color legend changeable
                    seq_cov_str += '<mark2>' + protein_seq[j] + '</mark2>'
                elif 0.6 * max_freq <= freq_array[j] < 0.8 * max_freq:  # color legend changeable
                    seq_cov_str += '<mark1>' + protein_seq[j] + '</mark1>'
                else:
                     seq_cov_str += '<mark>' + protein_seq[j] + '</mark>'
            else:
                seq_cov_str += '<mark5>' + protein_seq[j] + '</mark5>'
        seq_cov_str += '\n'
    return seq_cov_str,seq_cov


def hashcolor(s):

    return Turbo256[hash(s) % 256]


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

    ptm_index_dict = ptm_map_result[protein_entry]
    ptm_index_sort = sorted([(idx,each) for each in ptm_index_dict for idx in ptm_index_dict[each]])

    # bokeh plot, hovertool
    hover = HoverTool(names=['rec'],tooltips=[('domain', '@domain'), ('start position', '@start'),('end position','@end'),('domain coverage','@coverage{:.1%}'),])
    # initiate bokeh figure
    p = figure(x_range=(-10,protein_len),
               y_range=(0,2),
               tools=['pan', 'box_zoom', 'wheel_zoom', 'save',
                      'reset', hover],
               plot_height=500, plot_width=1200,
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
    numpy_zero_array = np.zeros((800, protein_len)) # mask numpy array for text plotting
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
        x_move = int(protein_len/12)
        while True: # keep moving down if text are too close
            nonzero_count = np.count_nonzero(numpy_zero_array[790+y_offset:800+y_offset,each_idx+x_offset:each_idx+x_move+x_offset])
            if nonzero_count == 0:
                # print (ptm,each_idx,x_offset,y_offset)
                new_ptm_x.append(each_idx+x_offset)
                new_ptm_y.append((25+y_offset)/200*2)
                numpy_zero_array[790+y_offset:800+y_offset,each_idx+x_offset:each_idx+x_move+x_offset] += 1
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


def html_compile(html_template,
                 domain_return,
                 seq_cov_gen_return,
                 uniprot_id,
                 protein_info_dict,
                 html_out='test.html'):

    domain_js_script, domain_div = domain_return
    seq_cov_str, coverage = seq_cov_gen_return
    smart_url = 'https://smart.embl.de/smart/show_motifs.pl?ID=' + uniprot_id

    # read html template
    html_template_f = open(html_template)
    html_template = html_template_f.read()
    html_template_f.close()

    # write new html with js script and divs
    # domain graph div
    new_html = html_template.replace('<!-- COPY/PASTE domain coverage SCRIPT HERE -->', domain_js_script). \
        replace('<!-- INSERT domain DIVS HERE -->', domain_div)

    # seq coverage
    new_html = new_html.replace('<!-- COPY/PASTE seq coverage value HERE-->',str(coverage))\
        .replace('<!-- COPY/PASTE seq coverage str HERE-->',seq_cov_str)

    # id and gene
    new_html = new_html.replace('<!-- UniprotID -->',
                                protein_info_dict[uniprot_id][0] + '  (' + protein_info_dict[uniprot_id][
                                    1].rstrip(' ') + ')').replace('<!-- SMART URL -->', smart_url)

    with open(html_out, 'w', newline='\n') as f_w:
        f_w.write(new_html)


def peptide_reader(input_file):

    """
    read psm into a dictionary, {peptide:[psm1,psm2...]}
    :param input_file:
    :return:
    """
    regex_pat = '\w{1}\[\d+\.?\d+\]'
    psm_dict = defaultdict(list)

    with open(input_file,'r') as f_o:
        f_split = f_o.read().replace('\r','').split('\n')

        for each in f_split:
            if '.' in each:
                pep_mod = each[2:-2]
                reg_sub = re.sub(regex_pat, my_replace, pep_mod)

                psm_dict[reg_sub].append(pep_mod)
    return psm_dict


if __name__ == '__main__':
    import json
    import pickle

    # pep_file = 'F:/matrisomedb2.0/peptide_list.txt'
    # psm_dict = peptide_reader(pep_file)
    # pickle.dump(psm_dict,open('F:/matrisomedb2.0/psm_dict_allpeptides_0906.p','wb'))

    # psm_dict = pickle.load(open('F:/matrisomedb2.0/psm_dict_allpeptides_0906.p','rb'))

    # print (sum([len(psm_dict[each]) for each in psm_dict]))

    protein_dict = fasta_reader('F:/matrisomedb2.0/mat.fasta')

    protein_info_dict = protein_info_from_fasta('F:/matrisomedb2.0/mat.fasta')
    # protein_freq_dict= peptide_map(psm_dict,protein_dict)[0]
    # pickle.dump(protein_freq_dict,open('F:/matrisomedb2.0/prot_freq_dict_0907.p','wb'),protocol=5)
    # prot_psm_list_dict = peptide_map(psm_dict,protein_dict)[2]
    # pickle.dump(prot_psm_list_dict, open('F:/matrisomedb2.0/prot_psm_list_dict_0908.p', 'wb'), protocol=5)
    protein_psm_list_dict = pickle.load(open('F:/matrisomedb2.0/prot_psm_list_dict_0908.p', 'rb'))


    prot_freq_dict = pickle.load(open('F:/matrisomedb2.0/prot_freq_dict_0907.p','rb'))
    # print ('prot_freq_dict reading done.')

    # ptm mappig
    psm_list = protein_psm_list_dict['Q61001']
    print (psm_list)
    ptm_map_result = ptm_map(psm_list,protein_dict)[0]

    # pickle.dump(ptm_map_result,open('F:/matrisomedb2.0/ptm_map_result_0907.p','wb'),protocol=5)

    seq_cov = seq_cov_gen(prot_freq_dict['Q61001'],ptm_map_result['Q61001'], protein_dict['Q61001'])
    #
    # read domain from json file
    with open('F:/matrisomedb2.0/smart_domain.json') as f_o:
        info_dict = json.load(f_o)
    prot_len_dict = {each:len(protein_dict[each]) for each in info_dict if each in protein_dict}
    print (sorted(prot_len_dict.items(), key=lambda x:x[1]))

    domain_cov = domain_cov_ptm(prot_freq_dict,
                                ptm_map_result,
                                info_dict,
                                'Q61001')

    html_template = 'F:/matrisomedb2.0/domain_seq_cov_html_template.html'

    html_compile(html_template,domain_cov,seq_cov,'Q61001',protein_info_dict,
                 html_out='F:/matrisomedb2.0/Q61001_test.html')

    # crawl SMART db