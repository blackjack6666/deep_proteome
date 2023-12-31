"""
source code of some new functions in MatrisomeDB2.0
"""

from multiprocessing_naive_algorithym import *
from aho_corasick import automaton_trie, automaton_matching
from bokeh.models import HoverTool, ColumnDataSource, FactorRange, LinearColorMapper,ColorBar,BasicTicker,PrintfTickFormatter, Plot, Rect,CustomJS, Legend, LegendItem,SingleIntervalTicker, Label, LabelSet, TableColumn, DataTable, HTMLTemplateFormatter, Button
from bokeh.palettes import Spectral7, Viridis, Plasma, Blues9, Turbo256
from bokeh.transform import factor_cmap
from bokeh.plotting import figure,output_file,save
from bokeh.io import save, output_file, show
from bokeh.embed import components
from bokeh.layouts import column
import pandas as pd

ptm_map_dict = {'Q\\[129\\]':'Gln deamidation','N\\[115\\]':'ASN deamidation',
                'Q\\[111\\]':'Gln to pyroGln','C\\[143\\]':'selenocysteine',
                'M\\[15\\.9949\\]':'Met oxidation','P\\[15\\.9949\\]':'Pro hydroxylation',
                'K\\[15\\.9949\\]':'Lys hydroxylation','n\\[42\\.0106\\]':'N-term acetylation',
                'C\\[57\\.0215\\]':'Cys redu-alky','R\\[0\\.9840\\]':'Arg deamidation','Y\\[79\\.9663\\]':'Phospho-Tyr',
                'T\\[79\\.9663\\]':'Phospho-Thr', 'S\\[79\\.9663\\]':'Phospho-Ser'}
regex_list = ['M\\[15\\.9949\\]','P\\[15\\.9949\\]','K\\[15\\.9949\\]','n\\[42\\.0106\\]','C\\[57\\.0215\\]',
              'R\\[0\\.9840\\]','T\\[79\\.9663\\]','S\\[79\\.9663\\]','Y\\[79\\.9663\\]']


button_ptm_table_js_code = """
function table_to_csv(source) {
    const columns = Object.keys(source.data)
    const nrows = source.get_length()
    const lines = [columns.join(',')]

    for (let i = 0; i < nrows; i++) {
        let row = [];
        for (let j = 0; j < columns.length-1; j++) {
            const column = columns[j]
            row.push(source.data[column][i].toString())
        }
        lines.push(row.join(','))
    }
    return lines.join('\n').concat('\n')
}


const filename = prot+'_'+'PTMs_table.csv'
const filetext = table_to_csv(source)
const blob = new Blob([filetext], { type: 'text/csv;charset=utf-8;' })

//addresses IE
if (navigator.msSaveBlob) {
    navigator.msSaveBlob(blob, filename)
} else {
    const link = document.createElement('a')
    link.href = URL.createObjectURL(blob)
    link.download = filename
    link.target = '_blank'
    link.style.visibility = 'hidden'
    link.dispatchEvent(new MouseEvent('click'))
}
"""


def fasta_reader(fasta_file_path):

    with open(fasta_file_path, 'r') as file_open:
        file_split = file_open.read().split('\n>')

    return {each.split('\n')[0].split('|')[1]: ''.join(each.split('\n')[1:]) for each in file_split}


def protein_info_from_fasta(fasta_path):
    """
    get protein name, gene name, entry name, and description
    :param fasta_path:
    :return:
    """
    info_dict = {}
    with open(fasta_path,'r') as f:
        for line in f:
             if line.startswith('>'):
                protein_id = line.split('|')[1]
                cls = line.split('|')[0].split('>')[1]
                # print (protein_id)
                description = ' '.join(line.split('OS=')[0].split(' ')[1:])

                gene_name = line.split('GN=')[1].split(' ')[0].rstrip('\n') if 'GN=' in line else 'N/A'
                info_dict[protein_id] = (gene_name,description,cls)
    return info_dict

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

    # ptm mapping, n-term mod would give 2 index, need to fix. ---> might fix by add n? at the begining of regex str
    ptm_start = time.time()
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


def ptm_map2(peptide_psm_dict, protein_dict, regex_set=regex_list):
    time_start = time.time()

    id_ptm_idx_dict = {}  # {protein_id:{ptm1:nonzero_index_array,ptm2:nonzero_index_array,...}}
    regex_set = set(regex_set)
    if 'C\\[57\\.0215\\]' in regex_set:
        regex_set.remove('C\\[57\\.0215\\]')

    print(regex_set)
    # aho mapping
    id_list, seq_list = extract_UNID_and_seq(protein_dict)
    seq_line = creat_total_seq_line(seq_list, sep="|")
    ptm_index_line_dict = {each: zero_line_for_seq(seq_line) for each in regex_set}
    separtor_pos_array = separator_pos(seq_line)
    aho_result = automaton_matching(automaton_trie([pep for pep in peptide_psm_dict]), seq_line)

    # ptm mapping, n-term mod would give 2 index, need to fix ---> might fix by add n? in regex pattern to find n-term
    ptm_start = time.time()
    for tp in aho_result:
        matched_pep = tp[2]  # without ptm site
        for psm in peptide_psm_dict[matched_pep]:
            for ptm in regex_set:  # check each ptm, mask other ptms
                new_psm = re.sub('n?\[\d+\.?\d+\]', '', psm.replace(ptm.replace('\\', ''), '*')). \
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
    print(f'ptm mapping took {time.time() - ptm_start}')

    # get ptm freq array and index
    for i in range(len(separtor_pos_array) - 1):
        id_ptm_idx_dict[id_list[i]] = {ptm:
                                           np.nonzero(ptm_index_line_dict[ptm][
                                                      separtor_pos_array[i] + 1:separtor_pos_array[i + 1]])[0]
                                       for ptm in regex_set}

    print(f'ptm script took {time.time() - time_start}s')
    return id_ptm_idx_dict


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
    print (len(freq_array),len(protein_seq))
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


def domain_cov_ptm_csv(prot_freq_dict, ptm_map_result, domain_pos_dict,protein_entry:str,output_name):
    """
    -----
    output CSVs for domain coverage and domain ptms
    :param prot_freq_dict:
    :param ptm_map_result:
    :param domain_pos_dict:
    :param protein_entry:
    :return:
    """

    if protein_entry not in domain_pos_dict:
        return f'{protein_entry} no SMART domain available'
    else:
        freq_array = prot_freq_dict[protein_entry]
        domain_dict = domain_pos_dict[protein_entry]
        ptm_index_dict = ptm_map_result[protein_entry]
        info_list = []
        ptm_info_list = []
        for each_domain in domain_dict:
            for each_tp in domain_dict[each_domain]:
                start, end = each_tp[0], each_tp[1]
                # calculate domain coverage
                coverage = np.count_nonzero(freq_array[start - 1:end]) / (end - start + 1)*100
                info_list.append([start, end, coverage, each_domain])
                # domain ptms
                for ptm in ptm_index_dict:
                    ptm_index_ary = ptm_index_dict[ptm]
                    domain_ptm_idx = np.where((ptm_index_ary>=start-1)&(ptm_index_ary<=end-1),ptm_index_ary,0)
                    ptm_info_list.append([start,end,ptm_map_dict[ptm],np.array2string(domain_ptm_idx[domain_ptm_idx!=0]+1,separator=', ')[1:-1],each_domain])
        if info_list == []:
            return f'{protein_entry} no SMART domain available'
        else:
            df_cov = pd.DataFrame(info_list,columns=['start','end','seq cov%','domain name'])
            df_ptm = pd.DataFrame(ptm_info_list,columns=['start','end','PTM','position','domain name'])
            df_cov.to_csv(output_name+'_DomainCov.csv')
            df_ptm.to_csv(output_name+'_DomainPTM.csv')


def domain_cov_ptm(prot_freq_dict, ptm_map_result, domain_pos_dict,protein_entry:str, data_source='sample'):
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
    if protein_entry not in domain_pos_dict:
        return '','No SMART domain available'
    else:
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
        if info_list ==[]:
            return '', 'No SMART domain available'
        else:
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
            ptm_count = sum([len(ptm_index_dict[ptm]) for ptm in ptm_index_dict])
            # bokeh plot, hovertool
            hover = HoverTool(names=['rec'],tooltips=[('domain', '@domain'), ('start position', '@start'),('end position','@end'),('domain coverage','@coverage{:.1%}'),])
            # initiate bokeh figure
            p = figure(x_range=(-10,protein_len),
                       y_range=(0,2),
                       tools=['pan', 'box_zoom', 'wheel_zoom', 'save',
                              'reset', hover],
                       plot_height=500, plot_width=700,
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

            if data_source=='sample': # if plotting from sample-specific data, label PTMs below domains
                # adjusted PTM text coordinates calculation
                numpy_zero_array = np.zeros((1500, protein_len)) # mask numpy array for text plotting
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
                    x_move = int(protein_len/8)  # 12 to 8 when plot width is 1200 to 700
                    while True: # keep moving down if text are too close
                        nonzero_count = np.count_nonzero(numpy_zero_array[1490+y_offset:1500+y_offset,each_idx+x_offset:each_idx+x_move+x_offset])
                        if nonzero_count == 0:
                            # print (ptm,each_idx,x_offset,y_offset)
                            new_ptm_x.append(each_idx+x_offset)
                            new_ptm_y.append((25+y_offset)/200*2)
                            numpy_zero_array[1490+y_offset:1500+y_offset,each_idx+x_offset:each_idx+x_move+x_offset] += 1
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
            else: # if from global data, do not label PTMs
                label = Label(x=1,y=0.10,text=str(ptm_count)+' PTMs occurrence in total\nSee details from PTMs table', text_font_size = '15px', text_font = 'Tahoma')
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


def domain_cov_ptm2(prot_freq_dict, ptm_map_result, domain_pos_dict,protein_entry:str):
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
    if protein_entry not in domain_pos_dict:
        return '','No SMART domain available'
    else:
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
        ptm_clean_index_dict = {ptm:ptm_clean(ptm_index_dict[ptm]) for ptm in ptm_index_dict}


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
        ptm_line_y = np.linspace(0.1,0.5,len(ptm_index_dict))
        ptm_y_map_dict = {ptm:y for ptm, y in zip([k for k in ptm_index_dict.keys()],ptm_line_y)}

        for ptm in ptm_index_dict:
            for x in ptm_index_dict[ptm]:
                p.line(x=[x,x],y=[0.5, 0.5-ptm_y_map_dict[ptm]])

        for ptm in ptm_clean_index_dict:
            for label_x in ptm_clean_index_dict[ptm]:
                label = Label(x=label_x, y=ptm_y_map_dict[ptm], text=ptm_map_dict[ptm] + '\n' + str(label_x + 1), text_font_size='10px',
                              text_align='center', text_font='Tahoma')
                p.add_layout(label)
        #     # ptms.append(ptm.replace('\\', ''))
        #     ptms.append(ptm)
        #     ptm_x.append(each_idx)
        #     ptm_y.append(0.5)
        #     x_offset,y_offset = 0,0
        #     x_move = int(protein_len/12)
        #     while True: # keep moving down if text are too close
        #         nonzero_count = np.count_nonzero(numpy_zero_array[790+y_offset:800+y_offset,each_idx+x_offset:each_idx+x_move+x_offset])
        #         if nonzero_count == 0:
        #             # print (ptm,each_idx,x_offset,y_offset)
        #             new_ptm_x.append(each_idx+x_offset)
        #             new_ptm_y.append((25+y_offset)/200*2)
        #             numpy_zero_array[790+y_offset:800+y_offset,each_idx+x_offset:each_idx+x_move+x_offset] += 1
        #             break
        #         else:
        #             # print ('moving down')
        #             # x_offset += 50 # value to move right
        #             y_offset -= 12 # value to move down
        #
        # # label ptm and connect to protein domains
        # for x,y,x_,y_,ptm in zip(new_ptm_x,new_ptm_y,ptm_x,ptm_y,ptms):
        #     p.line(x=[x_+1,x+1],y=[y_,y+0.1],line_width=1,color='black',alpha=0.3) # connect domain with text
        #     label = Label(x=x,y=y,text=ptm_map_dict[ptm]+'\n'+str(x_+1),text_font_size='10px', text_align='center', text_font='Tahoma')
        #     p.add_layout(label)
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


def psmlist_todict(psm_list):
    regex_pat = '\w{1}\[\d+\.?\d+\]'
    psm_dict = defaultdict(list)

    for psm in psm_list:
        reg_sub = re.sub(regex_pat, my_replace, psm)
        psm_dict[reg_sub].append(psm)
    return psm_dict


def ptm_clean(ptm_indx_list):
    threshold = 10
    ptm_indx_list = sorted(ptm_indx_list)
    diff_list = np.insert(np.diff(ptm_indx_list), 0, 100)
    diff_list = [int(i / 10) for i in diff_list]
    order_list = list(zip(ptm_indx_list, diff_list))
    print(order_list)
    new_ptm_list = []
    for each in order_list:
        if each[1] > 0:
            new_ptm_list.append(each[0])
    return new_ptm_list


def ptm_table_bokeh(id_ptm_idx_dict, protein_dict, protein_info_dict, output_html):
    """
    plot a bokeh table to show ptms frequency
    :param id_ptm_idx_dict:
    :param protein_dict:
    :param uniprot_gene_dict:
    :return:
    """
    import pandas as pd

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
    columns.append(TableColumn(field='Uniprot ID', title='Uniprot ID'))
    columns += [TableColumn(field=each,title=each)
                for each in df.columns[1:] if each!='url']
    table = DataTable(source=source,columns=columns, width=1000, height=600, editable=True)
    output_file(output_html)
    save(table)


def ptm_table_bokeh2(id_ptm_idx_dict, output_base, js_code):
    """
    plot a bokeh table to show ptms frequency
    :param id_ptm_idx_dict:
    :param protein_dict:
    :param uniprot_gene_dict:
    :return:
    """
    import pandas as pd

    for prot in id_ptm_idx_dict:
        print (prot)
        ptm_index_dict = id_ptm_idx_dict[prot]
        ptms, positions = [ptm_map_dict[k] for k in ptm_index_dict.keys()], [[each+1 for each in ptm_index_dict[ptm]] for ptm in ptm_index_dict]
        df = pd.DataFrame(dict(PTMs=ptms, Positions=positions))
        # print (df)
        source = ColumnDataSource(df)

        columns = [TableColumn(field=each,title=each)
                    for each in df.columns]
        table = DataTable(source=source,columns=columns, width=2000, height=250, editable=True)
        # button to download table as tsv/excel
        button = Button(label="Download PTM table as CSV", button_type="success", width=300)
        button.js_on_event("button_click", CustomJS(args=dict(source=source, prot=prot),
                                                    code=js_code))
        controls = column(table,button)

        output_file(output_base+prot+'_ptmtable.html')
        save(controls)


def ptm_table_bokeh3(sample_data, output_base, js_code):
    """
    show ptm positions for each protein
    :param id_ptm_idx_dict:
    :param protein_dict:
    :param uniprot_gene_dict:
    :return:
    """
    import pandas as pd

    for sample in sample_data:


        ptm_id_index_dict = sample_data[sample]['ptm']

        for prot in ptm_id_index_dict:
            print (prot)
            # seq = protein_dict[prot]
            info_dict = {}
            # info_dict['Amino acid sequence'] = [aa for aa in seq]
            # info_dict['Amino acid position'] = range(1,len(seq)+1)

            for ptm in ptm_id_index_dict[prot]:

                info_dict[ptm_map_dict[ptm]] = [each+1 for each in ptm_id_index_dict[prot][ptm]]


            df = pd.DataFrame(dict(PTMs=[k for k in info_dict.keys()],Positions=[v for v in info_dict.values()]))
            # print (df)
            source = ColumnDataSource(df)

            columns = [TableColumn(field=each,title=each)
                        for each in df.columns]
            table = DataTable(source=source,columns=columns, width=1000, height=400, editable=True)

            # add button to download table as tsv/excel
            button = Button(label="Download PTM table as CSV", button_type="success", width=300)
            button.js_on_event("button_click", CustomJS(args=dict(source=source,prot=sample+'_'+prot),
                                                        code=js_code))
            controls = column(table, button)

            output_file(output_base+sample.replace('/','-')+'_'+prot+'_ptmtable.html')
            save(controls)


if __name__ == '__main__':

    import json
    import pickle
    import os

    # single protein test
    """
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
    """

    ### matrisomeDB batch
    with open('F:/matrisomedb2.0/smart_domain_0908.json') as f_o:
        info_dict = json.load(f_o)

    protein_dict = fasta_reader('F:/matrisomedb2.0/mat.fasta')

    protein_info_dict = protein_info_from_fasta('F:/matrisomedb2.0/mat.fasta')

    # global_protein_psm_dict = json.load(open('F:/matrisomedb2.0/global_protein_psm.dict_fromsample.json','r'))
    # sample_prot_psm_dict = json.load(open('F:/matrisomedb2.0/sample_protein_psm_dict_3.json','r')) # sample name has illgal charcters
    # pickle.dump(sample_prot_psm_dict['Pancreatic Ductal Adenocarcinoma Xenograft (BxPC3)']['P21980'],open(r'F:\matrisomedb2.0/data_for_test/sample_psm_list.p','wb'),protocol=5)

    ## map glob psm into freq_array_dict and ptm index dict
    # all_psm = pickle.load(open('F:/matrisomedb2.0/all_psm.p','rb'))
    # all_psm = pickle.load(open('F:/matrisomedb2_0_revise/all_psm_result9.p','rb'))
    # print (len(all_psm))
    # all_psm_dict = psmlist_todict(all_psm)
    # print ('convet psm list to pep_psm dict')
    # glob_prot_freq_dict = peptide_map(all_psm_dict,protein_dict)[0]
    # pickle.dump(glob_prot_freq_dict, open('F:/matrisomedb2_0_revise/glob_prot_freq_dict_result9.p', 'wb'),protocol=5)
    # glob_prot_ptm_ind_dict = ptm_map2(all_psm_dict,protein_dict,regex_list)
    # pickle.dump(glob_prot_ptm_ind_dict,open('F:/matrisomedb2_0_revise/glob_prot_ptm_ind_dict.p','wb'),protocol=5)

    ## map sample psm into freq array dict and ptm index dict
    # sample_prot_psm_dict = json.load(open('F:/matrisomedb2_0_revise/sample_protein_psm_dict_result9.json','r'))
    # sample_data = {}
    # for sample in sample_prot_psm_dict:
    #     print (sample)
    #     sample_dict = {}
    #     psm_list = []
    #     for prot in sample_prot_psm_dict[sample]:
    #         psm_list+=sample_prot_psm_dict[sample][prot]
    #     if len(psm_list)!=0:
    #         psm_dict = psmlist_todict(psm_list)
    #         prot_freq_dict = peptide_map(psm_dict,protein_dict)[0]
    #         sample_dict['freq'] = prot_freq_dict
    #         prot_ptm_ind_dict = ptm_map2(psm_dict,protein_dict,regex_list)
    #         sample_dict['ptm'] = prot_ptm_ind_dict
    #         sample_data[sample] = sample_dict
    # pickle.dump(sample_data,open('F:/matrisomedb2_0_revise/sample_result9.data','wb'),protocol=5)

    # glob_prot_freq_dict = pickle.load(open('F:/matrisomedb2.0/glob_prot_freq_dict.p','rb'))
    # glob_ptm_map = pickle.load(open('F:/matrisomedb2.0/glob_prot_ptm_ind_dict.p','rb'))
    # html_tempalte = open(r'F:\matrisomedb2.0\test/domain_seq_cov_html_template_0909.html', 'r')
    # html_tempalte_read = html_tempalte.read()
    # html_tempalte.close()

    ## domain cov/PTM CSVs
    folder = 'F:/matrisomedb2_0_revise/domain_cov_ptm_csv/'
    # for prot in glob_prot_freq_dict:
    #     print (prot)
    #     domain_cov_ptm_csv(glob_prot_freq_dict,glob_ptm_map,info_dict,prot,folder+prot)

    sample_data = pickle.load(open('F:/matrisomedb2_0_revise/sample_result9.data','rb'))
    for sample in sample_data:
        print (sample)
        prot_freq_dict, ptm_map = sample_data[sample]['freq'], sample_data[sample]['ptm']
        for prot in prot_freq_dict:
            domain_cov_ptm_csv(prot_freq_dict,ptm_map,info_dict,prot,folder+sample.replace('/', '_').replace('\u0394','') + '_' + prot)


    # PTM tables sample
    # js_code = open('F:/matrisomedb2.0/db_script/ptm_button.js').read()
    # sample_data = pickle.load(open('F:\matrisomedb2.0/sample.data','rb'))
    # output_base = r'F:\matrisomedb2.0/table_htmls_2/'
    # ptm_table_bokeh2(sample_data,protein_dict,protein_info_dict,output_base)
    # ptm_table_bokeh3(sample_data,output_base='F:/matrisomedb2.0/table_htmls/',js_code=js_code)

    # ptm table global

    # ptm_table_bokeh2(glob_ptm_map,output_base='F:/matrisomedb2.0/table_htmls/',js_code=js_code)

    # domain html generation
    """
    for sample in sample_prot_psm_dict:
        print (sample)
        for prot in sample_prot_psm_dict[sample]:
            if os.path.exists(r'F:\matrisomedb2.0\domain_htmls/'+sample+'_'+prot+'_domaincov.html'):
                print (sample,prot,'exists')
                continue

            else:
                start = time.time()
                print (prot)
                # new_protein_dict = {prot:protein_dict[prot],'XXX':'XXX'}
                sample_psm_list = sample_prot_psm_dict[sample][prot]
                sample_psm_dict = psmlist_todict(sample_psm_list)
                # protein_dict = {prot:protein_dict[prot],'XXX':'XXX','AAA':'AAA'}
                # global_psm_list = global_protein_psm_dict[prot]
                # glob_psm_dict = psmlist_todict(global_psm_list)

                sample_prot_freq_dict = peptide_map(sample_psm_dict,protein_dict)[0]
                # print (prot,len(sample_prot_freq_dict[prot]),len(protein_dict[prot]))
                # glob_prot_freq_dict = peptide_map(glob_psm_dict,protein_dict)[0]

                sample_ptm_map = ptm_map(sample_psm_list,protein_dict)[0]
                # glob_ptm_map = ptm_map(global_psm_list,protein_dict)[0]
                if len(sample_prot_freq_dict[prot]) != len(protein_dict[prot]):
                    print (prot, 'NON EQUAL LENGTH!')
                    continue
                else:
                    sample_seq_cov = seq_cov_gen(sample_prot_freq_dict[prot],sample_ptm_map[prot],protein_dict[prot])
                    glob_seq_cov = seq_cov_gen(glob_prot_freq_dict[prot],glob_ptm_map[prot],protein_dict[prot])

                    sample_domain_cov = domain_cov_ptm(sample_prot_freq_dict,sample_ptm_map,info_dict,prot,data_source='sample')
                    glob_domain_cov = domain_cov_ptm(glob_prot_freq_dict, glob_ptm_map, info_dict, prot,data_source='global')

                    new_html = html_tempalte_read.replace('<!-- COPY/PASTE domain coverage SCRIPT HERE -->',sample_domain_cov[0]).\
                        replace('<!-- COPY/PASTE domain coverage global SCRIPT HERE -->',glob_domain_cov[0]).\
                        replace('<!-- INSERT domain DIVS HERE -->',sample_domain_cov[1]).\
                        replace('<!-- INSERT domain DIVS global HERE -->', glob_domain_cov[1]).\
                        replace('<!-- COPY/PASTE seq coverage value HERE-->',str(sample_seq_cov[1])).\
                        replace('<!-- COPY/PASTE seq coverage value global HERE-->',str(glob_seq_cov[1])).\
                        replace('<!-- COPY/PASTE seq coverage str HERE-->',sample_seq_cov[0]).\
                        replace('<!-- COPY/PASTE seq coverage str global HERE-->',glob_seq_cov[0]).\
                        replace('<!-- UniprotID -->', protein_info_dict[prot][0] + '  (' + protein_info_dict[prot][
                                            1].rstrip(' ') + ')').replace('<!-- sample_type -->',sample).\
                        replace('<!-- 3d cov URL -->',sample+'_'+prot+'_3dcov.html').replace('<!-- 3d cov global URL -->',prot+'_3dcov.html').\
                        replace('<!-- PTM table URL -->', sample+'_'+prot+'_ptmtable.html').replace('<!-- PTM table global URL -->',prot+'_ptmtable.html')
                    with open(r'F:\matrisomedb2.0\domain_htmls/'+sample.replace('/','-')+'_'+prot+'_domaincov.html','w') as f_o:
                        f_o.write(new_html)
                    print (f'time for {prot}: {time.time()-start}')
    """

    # domain html one protein test
    """
    start = time.time()
    out_put = r'F:\matrisomedb2.0/test/Basement membrane of lens capsule_P08572.html'
    sample = 'Basement membrane of lens capsule'
    prot = 'P08572'
    smart_url = 'https://smart.embl.de/smart/show_motifs.pl?ID=' + prot
    # sample_psm_list = pickle.load(open(r'F:\matrisomedb2.0/data_for_test/sample_psm_list.p','rb'))
    # sample_psm_dict = psmlist_todict(sample_psm_list)
    # protein_dict = {prot:protein_dict[prot],'XXX':'XXX','AAA':'AAA'}
    # global_psm_list = global_protein_psm_dict[prot]
    # glob_psm_dict = psmlist_todict(global_psm_list)

    # sample_prot_freq_dict = peptide_map(sample_psm_dict, protein_dict)[0]
    sample_prot_freq_dict = sample_data[sample]['freq']
    # print (prot,len(sample_prot_freq_dict[prot]),len(protein_dict[prot]))
    # glob_prot_freq_dict = peptide_map(glob_psm_dict,protein_dict)[0]

    # sample_ptm_map = ptm_map(sample_psm_list, protein_dict)[0]
    sample_ptm_map = sample_data[sample]['ptm']
    # glob_ptm_map = ptm_map(global_psm_list,protein_dict)[0]

    sample_seq_cov = seq_cov_gen(sample_prot_freq_dict[prot], sample_ptm_map[prot], protein_dict[prot])
    glob_seq_cov = seq_cov_gen(glob_prot_freq_dict[prot], glob_ptm_map[prot], protein_dict[prot])

    sample_domain_cov = domain_cov_ptm(sample_prot_freq_dict, sample_ptm_map, info_dict, prot, data_source='sample')
    glob_domain_cov = domain_cov_ptm(glob_prot_freq_dict, glob_ptm_map, info_dict, prot, data_source='global')

    new_html = html_tempalte_read.replace('<!-- COPY/PASTE domain coverage SCRIPT HERE -->', sample_domain_cov[0]). \
        replace('<!-- COPY/PASTE domain coverage global SCRIPT HERE -->', glob_domain_cov[0]). \
        replace('<!-- INSERT domain DIVS HERE -->', sample_domain_cov[1]). \
        replace('<!-- INSERT domain DIVS global HERE -->', glob_domain_cov[1]). \
        replace('<!-- COPY/PASTE seq coverage value HERE-->', str(sample_seq_cov[1])). \
        replace('<!-- COPY/PASTE seq coverage value global HERE-->', str(glob_seq_cov[1])). \
        replace('<!-- COPY/PASTE seq coverage str HERE-->', sample_seq_cov[0]). \
        replace('<!-- COPY/PASTE seq coverage str global HERE-->', glob_seq_cov[0]). \
        replace('<!-- UniprotID -->', protein_info_dict['P21980'][0] + '  (' + protein_info_dict['P21980'][
        1].rstrip(' ') + ')').replace('<!-- sample_type -->', sample). \
        replace('<!-- 3d cov URL -->', sample + '_' + prot + '_3dcov.html').replace('<!-- 3d cov global URL -->',
                                                                                    prot + '_3dcov.html'). \
        replace('<!-- PTM table URL -->', sample + '_' + prot + '_ptmtable.html').replace(
        '<!-- PTM table global URL -->', prot + '_ptmtable.html').replace('<!-- SMART URL -->',smart_url)
    with open(out_put,'w') as f_o:
        f_o.write(new_html)
    print(f'time for {prot}: {time.time() - start}')
    """