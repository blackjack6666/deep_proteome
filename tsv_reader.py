"""
include functions to read search results (peptide,protein,psm.tsv) and some plot functions (venn diagram) and text logo

"""

from collections import defaultdict
import re
import pandas as pd

def peptide_info_reader(psm_path:str):
    """
    read peptide sequence as key, [(filename, spectram number,mass, charge,retention time),()] as key
    :param psm_path:
    :return:
    """
    info_dict = defaultdict(list)
    with open(psm_path,'r') as f:
        for i in range(1):
            next(f)
        for line in f:
            line_split = line.split('\t')
            pep_seq = line_split[1]
            file_name = line_split[0].split('.')[0]
            spec_no = int(line_split[0].split('.')[-2])
            mass = float(line_split[6])
            charge = int(line_split[4])
            ret_time = float(line_split[5])
            info_dict[pep_seq].append([file_name,spec_no,mass,charge,ret_time])
    return info_dict


def peptide_charger_reader(pep_tsv):
    """
    get charge of the peptide
    :param pep_tsv:
    :return:
    """
    peptide_charge_dict = {}
    with open(pep_tsv,'r') as file_open:
        next(file_open)
        for line in file_open:
            line_split = line.split('\t')
            peptide_seq = line_split[0]
            charge = [line_split[2]] if ', ' not in line_split[2] else line_split[2].split(', ')
            peptide_charge_dict[peptide_seq] = charge
        return peptide_charge_dict


def id_pep_from_peptsv(pep_tsv):

    import numpy as np
    id_pep_dict = defaultdict(set)
    df = pd.read_csv(pep_tsv, sep='\t')
    for tp in df.itertuples():
        pep_seq = tp[1]
        prot_list = [tp[9].split('|')[1]]
        mapped_proteins = tp[-1]

        # print (mapped_proteins)
        if type(mapped_proteins)==str:
            if ',' in mapped_proteins:
                for each in mapped_proteins.split(', '):
                    prot_list.append(each.split('|')[1])
            else:
                prot_list.append(mapped_proteins.split('|')[1])
        # print (prot_list)
        for prot in prot_list:
            id_pep_dict[prot].add(pep_seq)
    return id_pep_dict

def peptide_phospho_reader(peptide_tsv_file, mod=0.9840): # 79.9663 is the delta mass of phosphorylation on STY
    pep_phos_dict = defaultdict()
    mod_str = "%.4f" % mod # keep last 4 digits
    with open(peptide_tsv_file) as file_open:
        for i in range(1):
            next(file_open)
        for line in file_open:
            pep_seq = line.split('\t')[0]

            pattern = re.compile('\d+\w{1}\('+mod_str+'\)')
            regex = re.findall(pattern, line)
            #print (regex)
            for ele in regex:
                if ele != '':
                    pep_phos_dict[pep_seq]=regex
    return pep_phos_dict


def modified_peptide_from_psm(psm_path):
    psm_list = []
    with open(psm_path, 'r') as f_open:
        next(f_open)
        for line in f_open:
            line_split = line.split('\t')
            match = re.search('\w{1}\[\d+\.?\d+\]',line)
            if match:
                psm_list.append(line_split[3])
            else:
                psm_list.append(line_split[2])
    return psm_list

def psm_ptm_reader(psm_tsv_file_list,gene_set,mod=0.9840):
    total_psm = 0
    ptm_psm_dict = defaultdict(list)
    mod_str = "%.4f" % mod  # keep last 4 digits
    for psm_tsv_file in psm_tsv_file_list:
        print (psm_tsv_file)
        with open(psm_tsv_file,'r') as file_open:
            next(file_open)
            for line in file_open:
                line_split = line.split('\t')
                psm = line_split[2]
                gene = line_split[29]
                total_psm+=1
                if gene in gene_set:
                    pattern = re.compile('\w{1}\(' + mod_str + '\)')
                    regex = re.findall(pattern, line)
                    if regex:
                        for ele in regex:
                            ptm_psm_dict[ele].append(psm)


    return ptm_psm_dict,total_psm


def protein_tsv_reader(protein_tsv_file):

    with open(protein_tsv_file, 'r') as file_open:
        next(file_open)
        return [line.split("\t")[3] for line in file_open]


def protein_reader(protein_tsv):
    df = pd.read_csv(protein_tsv, sep="\t")
    df_indis = df.dropna(subset=['Indistinguishable Proteins'])
    indis_protein_list = [prot.split('|')[1] for each in df_indis['Indistinguishable Proteins'].tolist() for prot in
                          each.split(", ")]
    protein_list = df.dropna(subset=['Protein ID'])['Protein ID'].tolist()

    return set(indis_protein_list + protein_list)

def protein_tsv_reader_no_contam(protein_tsv_file):
    with open(protein_tsv_file, 'r') as file_open:
        next(file_open)
        return [line.split("\t")[3] for line in file_open if 'contaminant' not in line.split("\t")[3]]

def peptide_counting(peptide_tsv_file):

    with open(peptide_tsv_file, 'r') as file_open:
        next(file_open)

        peptide_list = [line.split("\t")[0] for line in file_open]
    return peptide_list

def psm_reader(psm_path,fragpipe_ver=13.0):
    pep_spec_count_dict = defaultdict(int)
    ret_pep_dict = {}
    with open(psm_path, 'r') as f:
        for i in range(1):
            next(f)
        for line in f:
            line_split = line.split('\t')
            pep_seq = line_split[2] if fragpipe_ver==13.0 else line_split[1]
            retention_time = float(line_split[5])/60  if fragpipe_ver==13.0 else float(line_split[4])/60 # in minute
            pep_spec_count_dict[pep_seq]+=1
            ret_pep_dict[retention_time] = pep_seq
    return pep_spec_count_dict, ret_pep_dict

def venn_diagram_gen(dictionary, title=''): # parameter could be a dictionary of proteins or peptides from different samples {'sample1': [], 'sample2': []}
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2, venn3

    value_list_of_sets = [set(l) for l in dictionary.values()]
    sample_name_list = [n for n in dictionary.keys()]
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))

    if len(dictionary) == 2:  # two samples for venn diagram
        out = venn2(value_list_of_sets, set_labels=sample_name_list)

    elif len(dictionary) == 3:  # 3 samples for venn diagram
        out = venn3(value_list_of_sets, set_labels=sample_name_list)

    else:
        print ('Error: only 2 or 3 comparison for venn diagram are accepted in this script.')


    for t in out.set_labels: t.set_fontsize(22)
    for t in out.subset_labels: t.set_fontsize(20)

    fig.suptitle(title, fontsize=22)
    plt.show()


def venn_diagram_gen2(dictionary, title=None, png_output=None): # parameter could be a dictionary of proteins or peptides from different samples {'sample1': [], 'sample2': []}
    import matplotlib.pyplot as plt
    import matplotlib
    import venn

    if len(dictionary)>=6:
        raise ValueError('maximum groups to compare is 5')

    else:
        matplotlib.rcParams.update({'font.size': 15})

        venn_dict = {each:set(dictionary[each]) for each in dictionary}
        ax = venn.venn(venn_dict)
        if title:
            ax.set_title(title)
    if png_output:
        plt.savefig(png_output,dpi=300)
    plt.show()


def logomaker_from_stop_codon(plot_seq_list):
    """
    plot a logo graph
    :param plot_seq_list: a list of equal_len_sequence
    :return:
    """
    import logomaker as lm
    from matplotlib import pyplot as plt

    counts_mat = lm.alignment_to_matrix(plot_seq_list)
    prob_mat = lm.transform_matrix(counts_mat, normalize_values=True)
    print (prob_mat)
    # call lm.list_font_names() for font list
    logo = lm.Logo(prob_mat,
                   fade_probabilities=True,
                   alpha=0.9,
                   font_name='Berlin Sans FB Demi',
                   color_scheme='skylign_protein',
                   )

    # styling
    logo.ax.set_xticks(range(len(plot_seq_list[0])))
    logo.ax.set_xticklabels(['P'+str(i) for i in range(8,0,-1)]+['P'+str(i)+"'" for i in range(1,8)])
    logo.ax.axvline(7.5, color='black',linewidth=2,linestyle=':')
    logo.ax.set_ylabel('probability')
    plt.show()
    return prob_mat


def heatmap_gen(matrix):
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    fig, ax = plt.subplots(1,1, figsize =(10,10))
    g = sns.heatmap(matrix,
                    ax=ax,
                    cbar_kws={'label': 'ratio','ticks': range(0,120,20)},vmin=0,vmax=100, cmap="YlGnBu")


    ax.set_xticks(np.arange(0.5,15.5,1))
    ax.set_xticklabels(['P'+str(i) for i in range(8,0,-1)]+['P'+str(i)+"'" for i in range(1,8)], fontsize=15)
    ax.set_xlabel('Position', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=15)
    plt.tight_layout()
    plt.show()
    return ax


def myplot(x, y, s, bins=1000):
    """
    usage:   img, extent = myplot(x, y, 16)
             ax.imshow(img, extent=extent, origin='lower', cmap=cm.jet)
             plt.show()
    :param x: x coordinates
    :param y: y coordinates
    :param s: gaussian_filter sigma
    :param bins:
    :return:
    """
    import matplotlib.cm as cm
    from scipy.ndimage.filters import gaussian_filter

    heatmap, xedges, yedges = np.histogram2d(x, y, bins=bins)
    heatmap = gaussian_filter(heatmap, sigma=s)

    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    return heatmap.T, extent


def combined_proteintsv_map(combined_protein_tsv):
    """
    map spectra count from combined protein tsv file to each file
    :param combined_protein_tsv:
    :return:
    """
    info_dict = {}
    import pandas as pd
    df = pd.read_csv(combined_protein_tsv,sep='\t',index_col=False)
    # print (df.head)
    protein_list = df['Protein ID']

    for each_column in df.columns:
        if 'Total Spectral Count' in each_column:
            file_name = each_column.split(' ')[0]
            spec_count_list = df[each_column]
            protein_spec_dict = {i:j for i,j in zip(protein_list,spec_count_list) if j != 0}
            info_dict[file_name] = protein_spec_dict
    return info_dict


def plot_prot_combined_tsv(combined_protein_tsv):
    """
    plot identified proteins from combined protein tsv file
    :param combined_protein_tsv:
    :return:
    """
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    file_path = 'D:/data/Naba_deep_matrisome/matrisome coverage.xlsx'
    df_cov = pd.read_excel(file_path).drop_duplicates()
    ecm_prot_list = df_cov['protein_id']

    df = pd.read_csv(combined_protein_tsv,delimiter='\t')
    filt_col_list = [each for each in df.columns.tolist()
                     if 'Total Spectral Count' in each and 'Summarized Total' not in each][7:]
    filt_new_col_list = [each.split(' ')[0] for each in filt_col_list]
    filt_col_list.insert(0,'Protein ID')
    print (filt_col_list)

    filt_df = df[filt_col_list]
    filt_df = filt_df[filt_df['Protein ID'].isin(ecm_prot_list)]
    print (filt_df.shape)
    filt_df = filt_df.set_index('Protein ID')
    filt_df.columns = filt_new_col_list
    print (filt_df.columns)
    filt_df = filt_df.reindex(['18_2B05','18_2B1','18_2B2','18_2B4','18_2B18','18_2B20','18_2A'],axis=1)


    filt_df = filt_df.apply(np.log10)
    # fill infinite number with 0, infinite number comes from log10 of zero
    filt_df[filt_df==-np.inf] = 0
    print (filt_df.head())

    # print(filt_df.loc[filt_df['18_2B05'].idxmax()])
    # print (filt_df.head())
    #
    filt_df.T.plot.line(legend=False,figsize=(15,8),ylabel='Log10 spectra count')
    # plt.xticks(rotation=30)
    plt.show()

def protein_info_from_combined(combined_protein_tsv):
    """
    get protein name, gene name, entry name, gene name and Description

    :param combined_protein_tsv:
    :return:
    """
    info_dict = {}
    with open(combined_protein_tsv, 'r') as f_open:
        next(f_open)
        for line in f_open:
            line_split = line.split('\t')
            protein_id = line_split[3]
            protein_name = line_split[4]
            gene_name = line_split[5]
            description = line_split[10]
            info_dict[protein_id] = (protein_name,gene_name,description)
    return info_dict

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

def map_psm_file(psm_tsv):
    """
    map each peptide to the file where it's from
    :param psm_tsv:
    :return:
    """

    file_peptide_dict = defaultdict(list)
    with open(psm_tsv,'r') as f:
        next(f)
        for line in f:
            line_split = line.split('\t')
            file_name,psm = line_split[0].split('.')[0],line_split[2]
            file_peptide_dict[file_name].append(psm)
    return file_peptide_dict


def convert_uniprot_tonum(uniprot_id):

    return sum([ord(i) for i in uniprot_id])-320+ord(uniprot_id[-1])/10+ord(uniprot_id[-2])/100


def pep_xml_info(pepxml_file):
    with open(pepxml_file,'r') as f_read:
        f_read = f_read.read()
        f_split = f_read.split('<spectrum_query start_scan=')
        head, tail = f_split[0], f_read.split('</spectrum_query>')[-1]
        info_list = [each.split('</spectrum_query>')[0] for each in f_split[1:]]

    return head,tail,info_list

def pepxml_peptide_getter(pepxml):
    psm_hit_list = []
    with open(pepxml,'r') as f_open:
        file_split = f_open.read().split('peptide="')
        for each in file_split[1:]:
            peptide = each.split('"')[0]
            psm_hit_list.append(peptide)
    return psm_hit_list


def pep_mass_dist(pep_tsv_list:list, plot='kde'):
    from calculations_and_plot import protein_mass
    import seaborn as sns
    import matplotlib.pyplot as plt
    from statistics import mean

    ylabel_dict = {'kde':'Density','hist':'Frequency'}

    for pep_tsv in pep_tsv_list:
        sample_name = pep_tsv.split('/')[-2]
        mass_list = [protein_mass(each) for each in peptide_counting(pep_tsv)]
        print (sample_name,'average mass %f' % mean(mass_list))
        # draw density plot
        if plot == 'kde':
            sns.kdeplot(mass_list,linewidth=2, alpha=.5,label=sample_name)
        else:
            sns.histplot(mass_list,bins=100, alpha=.5, fill=False,label=sample_name)
    plt.legend(prop={'size': 10}, title='Time points')
    plt.title('%s plot native digestion multiple time points' % plot)
    plt.xlabel('Mass in Da')
    plt.ylabel('%s' % ylabel_dict[plot])

    plt.show()


def get_unique_peptide(list_of_peptsv:list):
    """
    from pep tsv file only get unique peptides compared with previous ones, e.g. in 4 hour sample, filter out peptides
    in 1h,2h and only retain peptides uniquely identified in 4h
    :param list_of_peptide:
    :return:
    """

    unique_peptide_dict = {}
    peptide_list = []
    for idx, val in enumerate(list_of_peptsv):

        file_name = val.split("/")[-2]
        print (file_name)
        unique_peptide_list = [each for each in peptide_counting(val) if each not in peptide_list]

        peptide_list += unique_peptide_list

        unique_peptide_dict[file_name] = unique_peptide_list

    return unique_peptide_dict


def map_k_r(psm_list, protein_dict):
    """
    map the start and end of each tryptic peptide
    :param psm_list:
    :param protein_dict:
    :param regex_dict: {regex:HEX color}
    :return:
    """

    import time
    import multiprocessing_naive_algorithym
    from aho_corasick import automaton_matching,automaton_trie

    id_kr_mapp_dict = {}

    # aho mapping
    id_list, seq_list = multiprocessing_naive_algorithym.extract_UNID_and_seq(protein_dict)
    seq_line = multiprocessing_naive_algorithym.creat_total_seq_line(seq_list, sep="|")
    zero_line = multiprocessing_naive_algorithym.zero_line_for_seq(seq_line)
    separtor_pos_array = multiprocessing_naive_algorithym.separator_pos(seq_line)

    aho_result = automaton_matching(automaton_trie(psm_list), seq_line)
    for tp in aho_result:
        # matched_pep = tp[2]  # without ptm site
        zero_line[tp[0]-1]+=1
        zero_line[tp[1]]+=1

    time_start = time.time()
    for i in range(len(separtor_pos_array)-1):
        zero_line_slice = zero_line[separtor_pos_array[i]+1:separtor_pos_array[i+1]]
        if np.count_nonzero(zero_line_slice) != 0:
            id_kr_mapp_dict[id_list[i]] = zero_line_slice


    return id_kr_mapp_dict


def kr_calculate(id_kr_mapp_dict,protein_dict):
    """
    calculate the number of K and R being cut for each protein
    :param id_kr_mapp_dict:
    :param protein_dict:
    :return:
    """

    from collections import Counter
    id_kr_count_dict = {}
    k_sum, r_sum = 0,0
    for prot in id_kr_mapp_dict:
        kr_index = np.nonzero(id_kr_mapp_dict[prot])[0]
        prot_seq = protein_dict[prot]
        kr_count_dict = Counter([prot_seq[idx] for idx in kr_index])
        k_sum += kr_count_dict['K']
        r_sum += kr_count_dict['R']
        id_kr_count_dict[prot] = {'K': kr_count_dict['K'],'R':kr_count_dict['R']}
    return id_kr_count_dict,k_sum,r_sum

if __name__=="__main__":
    from glob import glob
    import numpy as np
    import pandas as pd
    from pandas import ExcelWriter
    import matplotlib.pyplot as plt
    import math
    # path = 'C:/uic/lab/data/naba/search_result/*/peptide.tsv'
    # file_list = glob(path)
    # print (file_list)
    # # for f in file_list:
    # #     print (f.split('\\')[-2], sum([v for v in psm_reader(f)[0].values()]))
    # protein_info_dict = protein_info_from_combined('D:/data/Naba_deep_matrisome/11_11_combined_search/combined_protein.tsv')
    # info_dict = combined_proteintsv_map('D:/data/Naba_deep_matrisome/11_11_combined_search/combined_protein.tsv')
    # fasta_path = 'D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_tr_isoforms.fasta'
    # fasta_info_dict = protein_info_from_fasta(fasta_path)
    #
    # base_path = 'D:/data/pats/results/'
    # psm_tsv = ['trypsin_4h_sp_only_default_close/psm.tsv',
    #            'trypsin_4h_sp_isoform_default_close/psm.tsv',
    #            'trypsin_4h_sp_tr_default_close/psm.tsv',
    #            'trypsin_4h_sp_tr_isoform_default_close/psm.tsv']
    #
    # peptide_tsv = ['hek_trypsin_4hour_sp_isoforms/peptide.tsv',
    #            'hek_trypsin_4hour_sp_only/peptide.tsv',
    #            'hek_trypsin_4hour_sp_tr/peptide.tsv',
    #            'hek_trypsin_4hour_sp_tr_isoforms/peptide.tsv']

    pepxml = 'D:/data/deep_proteome/non_specfic_search/ct_4h/CT_37C_4h.pepXML'
    info_list = pep_xml_info(pepxml)[-1][-1]
    print (info_list)

    # venn_dict = {}
    #     # for each in psm_tsv:
    #     #     file = '_'.join(each.split('/')[0].split('_')[2:])
    #     #     psm_path = base_path + each
    #     #     protein_path = base_path + each.replace('psm', 'protein')
    #     #     peptide_tsv_path = base_path + each.replace('psm', 'peptide')
    #     #     protein_set = protein_reader(protein_path)
    #     #     print ('protein', len(protein_set))
    #     #     gene_set = set([fasta_info_dict[prot][0] for prot in protein_set])
    #     #     print ('gene', len(gene_set))
    #     #     psm_dict = psm_reader(psm_path)[0]
    #     #
    #     #     peptide_list = peptide_counting(peptide_tsv_path)
    #     #     print(len(peptide_list))
    #     #     psm_list = [pep + '_' + str(i) for pep in psm_dict for i in range(psm_dict[pep])]
    #     #     print(file, len(psm_list))
    #     #     venn_dict[file] = gene_set
    #     # venn_diagram_gen2(venn_dict,title='gene ids')
    #     # gene_list = [v for v in venn_dict.values()]
    # print ([gene for gene in gene_list[2] if gene not in gene_list[-1] and gene not in gene_list[1]])
    # venn_dict = {'163_3_dec': [k for k in info_dict['163_3_dec']],'163_3_extr': [k for k in info_dict['163_3_extr']]}
    # venn_dict = {each_file.split('\\')[-2]:peptide_counting(each_file)
    #              for each_file in file_list[:2]}
    # venn_diagram_gen(venn_dict,title='peptide compare')

    # for each in info_dict:
    #     print (each,len(info_dict[each]))

    # with ExcelWriter('D:/data/Naba_deep_matrisome/11_11_protein_ids.xlsx') as writer:
    #     for each in info_dict:
    #         info_list = [[prot,info_dict[each][prot],protein_info_dict[prot][0],protein_info_dict[prot][1],protein_info_dict[prot][2]]
    #                      for prot in info_dict[each]]
    #         df = pd.DataFrame(info_list, columns=['protein id', 'spectra count', 'entry name', 'gene name', 'description'])
    #         df.to_excel(writer,'%s' % each)
    # from calculations_and_plot import miss_cleavage_identify
    # for each_file in file_list:
    #     print (each_file)
    #     pep_list = peptide_counting(each_file)
    #     miss_cleav_dict = miss_cleavage_identify(pep_list,regex_pattern=r'(?:F|W|Y)\w+')
    #     print (float(np.count_nonzero([each for each in miss_cleav_dict.values()]))/len(pep_list))

    # df = pd.read_csv('D:/data/Naba_deep_matrisome/02152021_1/dash_info.csv')
    # uniprot_num_list = [convert_uniprot_tonum(uni_id) for uni_id in df['protein_id']]
    # df['uniprot_num'] = uniprot_num_list
    # df['log_protein_len'] = [math.log2(each) for each in df['length']]
    # df.to_csv('D:/data/Naba_deep_matrisome/02152021_1/dash_info.csv')
    # df = df.drop('Unnamed: 0', axis=1)


