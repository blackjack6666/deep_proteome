import pickle as ppp
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from collections import defaultdict
import multiprocessing_naive_algorithym
from aho_corasick import automaton_trie, automaton_matching
import glob
from parameters import aa_mass_table


def protein_mass(protein_seq):
    mass = 0
    for aa in protein_seq:
        mass += aa_mass_table[aa]
    mass += 18.01528 # add the weight of water

    return mass

def coverage_calculation(zero_line, sep_pos_array, ID_list):
    coverage_list_ordered = [] # append the coverage starting from left of zeroline all the way to right
    len_list = np.array([])
    protein_coverage_dict = {}
    total_seq_len = len(zero_line)-len(sep_pos_array)+2 #artifically added two separator positions into sep_pos_array so need to plus 2
    total_non_zero = np.count_nonzero(zero_line)
    overall_percentage = float(total_non_zero)/total_seq_len*100
    for i in range(len(sep_pos_array)-1): # iterate from the left of zeroline all the way to the right
        non_zero_for_one_protein = np.count_nonzero(zero_line[sep_pos_array[i]+1:sep_pos_array[i+1]])
        total_len_for_one_protein = len(zero_line[sep_pos_array[i]+1:sep_pos_array[i+1]])
        len_list = np.append(len_list, total_len_for_one_protein)
        coverage_for_one_protein = float(non_zero_for_one_protein)/total_len_for_one_protein*100
        coverage_list_ordered.append(coverage_for_one_protein)
        protein_coverage_dict[ID_list[i]] = coverage_for_one_protein
    coverage_array = np.array(coverage_list_ordered)
    identified_coverage_array = [x for x in coverage_array if x != 0]
    average_coverage = np.mean(coverage_array)
    identified_average_coverage = np.mean(identified_coverage_array)

    return overall_percentage, identified_coverage_array, coverage_list_ordered, protein_coverage_dict, average_coverage, identified_average_coverage, len_list


def identified_proteome_cov(aho_result, protein_dict):
    """
    -----
    get average coverage of identified proteome, and a dictionary that has coverage of each protein
    -----
    :param aho_result:
    :param protein_dict:
    :return:
    """
    prot_cov_dict = {}
    ID_list, seq_list = multiprocessing_naive_algorithym.extract_UNID_and_seq(protein_dict)

    seq_line = multiprocessing_naive_algorithym.creat_total_seq_line(seq_list,sep='|')
    separtor_pos_array = multiprocessing_naive_algorithym.separator_pos(seq_line)

    zero_line = multiprocessing_naive_algorithym.zero_line_for_seq(seq_line)
    for pos in aho_result:
        zero_line[pos[0]:pos[1] + 1] += 1

    for i in range(len(separtor_pos_array)-1): # iterate from the left of zeroline all the way to the right
        non_zero_for_one_protein = np.count_nonzero(zero_line[separtor_pos_array[i]+1:separtor_pos_array[i+1]])
        total_len_for_one_protein = len(zero_line[separtor_pos_array[i] + 1:separtor_pos_array[i + 1]])
        coverage_for_one_protein = float(non_zero_for_one_protein) / total_len_for_one_protein * 100
        prot_cov_dict[ID_list[i]] = coverage_for_one_protein

    # delete unidentified proteins
    iden_prot_cov_dict = {each:prot_cov_dict[each] for each in prot_cov_dict if prot_cov_dict[each] != 0}
    return np.mean([v for v in iden_prot_cov_dict.values()]), iden_prot_cov_dict, prot_cov_dict


def my_replace(match_obj):
    match_obj = match_obj.group()
    return match_obj[0]  # gives back the first element of matched object as string


def freq_array_and_PTM_index_generator(peptide_list, protein_seq_string,regex_pat='\w{1}\[\d+\.?\d+\]'):
    """
    map single protein seq
    :param peptide_list:
    :param protein_seq_string:
    :param regex_pat:
    :return:
    """
    import re
    freq_array = np.zeros(len(protein_seq_string))
    PTM_sites_counting = defaultdict(int)
    PTM_loc_list = []

    # reformat the peptide with PTM numbers into characters only
    new_pep_list = [re.sub(regex_pat, my_replace, pep) for pep in peptide_list]
    PTM_list = [re.findall(regex_pat, pep) for pep in peptide_list]
    # print (PTM_list)
    # calculation

    for pep, new_pep, PTM in zip(peptide_list, new_pep_list,PTM_list):  # PTM_list is list of list
        if new_pep in protein_seq_string:
            # print(new_pep)
            start_pos = protein_seq_string.find(new_pep)
            end_pos = start_pos + len(new_pep) -1
            # print (start_pos,end_pos,new_pep)
            freq_array[start_pos:end_pos + 1] += 1
            if PTM:  # the peptide has ptm site
                for ele in PTM:

                    PTM_index = pep.find(ele)
                   #PTM_site = pep[PTM_index] # single amino acid
                    PTM_sites_counting[ele] += 1
                    PTM_loc_list.append(start_pos+PTM_index)
    # print (PTM_sites_counting, PTM_loc_list)

    return freq_array, PTM_loc_list, PTM_sites_counting


def one_d_covearge_bar(html_template,peptide_list,protein_seq,output_html_path=None,screenshot=None):
    """
    generate 1d coverage bar on html
    :param html_template:
    :param peptide_list: 2d list containing result from different sample
    :param protein_seq:
    :return:
    """


    freq_array_2d = [freq_array_and_PTM_index_generator(pep_list,protein_seq)[0] for pep_list in peptide_list]
    # print (freq_array_2d)
    # rgb color list for 4 time points
    rgb_list = ['rgba(156,0,252,1) ','rgba(236,142,56,1) ','rgba(155,255,119,1) ','rgba(71,182,221,1) ']
    color_assign_dict = {i:j for i,j in zip(range(len(freq_array_2d)),rgb_list)}
    replace_str = ''
    for i in range(len(protein_seq)):
        percen = f'{i/float(len(protein_seq))*100:.1f}%'
        percen_next = f'{(i+1)/float(len(protein_seq))*100:.1f}%'
        # print (i,i+1,len(protein_seq),percen,percen_next)
        cov = False
        for j in range(len(freq_array_2d)):
            if freq_array_2d[j][i] == 0:
                continue
            else:
                cov = True
                replace_str += color_assign_dict[j]
                replace_str += percen
                replace_str += ', '
                replace_str += color_assign_dict[j]
                replace_str += percen_next
                replace_str += ', '
                break  # make sure the only color the peptide that shows first time
        if cov == False:
            replace_str += 'rgba(191,191,191,1) '  # not covered part show as grey
            replace_str += percen
            replace_str += ', '
            replace_str += 'rgba(191,191,191,1) '
            replace_str += percen_next
            replace_str += ', '
    replace_str = replace_str[:-2]  # delete ', ' in the end
    if output_html_path:
        f_open = open(html_template, 'r')
        new_html = f_open.read().replace('#insert_replace_here',replace_str)
        f_write = open(output_html_path,'w')
        f_write.write(new_html)
        if screenshot:
            from html2image import Html2Image
            hti = Html2Image()
            hti.screenshot(url=output_html_path,save_as=screenshot,size=(800,600))
        f_open.close()
        f_write.close()
    return replace_str


def whole_proteome_cov(aho_result, seq_line):
    """
    -----
    get proteome coverage from aho-corasick result
    -----
    :param aho_result:
    :param seq_line: long sequence line, no separator "|"
    :return:
    """

    total_seq_len = len(seq_line)
    zero_line = np.zeros(total_seq_len)
    for pos in aho_result:
        zero_line[pos[0]:pos[1] + 1] += 1
    non_zero_len = np.count_nonzero(zero_line)
    # print (non_zero_len)

    return float(non_zero_len)/total_seq_len*100


def length_KR_ratio(protein_dict):
    """
    calculate the ratio of protein length and K/R ratio
    :param protein_dict:
    :return:
    """

    return {each:float(protein_dict[each].count('K')+protein_dict[each].count('R'))/len(protein_dict[each])
            for each in protein_dict}


def specific_aa_length_ratio(protein_dict,spec_aa='KR'):
    """
    calculate the ratio of specific aa numbers to protein length
    :param protein_dict:
    :param spec_aa:
    :return:
    """
    product_dict = {}
    for each in protein_dict:
        count = 0
        for aa in spec_aa:
            count += protein_dict[each].count(aa)
        product_dict[each] = float(count)/len(protein_dict[each])
    return product_dict


def miss_cleavage_identify(peptide_list,regex_pattern=r'(?:K|R)\w+'):
    """
    calculate the ratio of peptides that has least one missed cleavage to all peptides
    :param peptide_list:
    :param regex_pattern: regex pattern for mapping enzyme specificity. Default is trypsin
    :return:
    """
    from re import search
    return {each: bool(search(regex_pattern,each)) for each in peptide_list}

# the following function returns a dictionary with uniprotID as key and aa frequency as value
def ID_number_line(zero_line, sep_pos_array, ID_list):
    ID_number_line_dict = {}
    for i in range(len(sep_pos_array)-1):
        ID_number_line_dict[ID_list[i]] = zero_line[sep_pos_array[i]+1:sep_pos_array[i+1]]
    return ID_number_line_dict



def ratio_of_coverage_to_sum(zero_line, sep_pos_array, ID_list):
    ID_ratio_dict = {}
    for i in range(len(sep_pos_array)-1):
        non_zero_for_one_protein = np.count_nonzero(zero_line[sep_pos_array[i]+1:sep_pos_array[i+1]])
        if non_zero_for_one_protein == 0: # check if the protein is identified or not
            continue
        else:
            total_len_for_one_protein = len(zero_line[sep_pos_array[i]+1:sep_pos_array[i+1]])
            coverage_for_one_identified = float(non_zero_for_one_protein)/total_len_for_one_protein*100
            sum_for_one_identified_protein = np.sum(zero_line[sep_pos_array[i]+1:sep_pos_array[i+1]])
            ratio = coverage_for_one_identified/sum_for_one_identified_protein
            ID_ratio_dict[ID_list[i]] = ratio
    return ID_ratio_dict

def read_into_dataframe(protein_coverage_dict, identified_coverage_array):
    uniprot_ID_list = protein_coverage_dict.keys()
    coverage_list = protein_coverage_dict.values()
    df1 = pd.DataFrame(dict(ID=uniprot_ID_list, coverage=coverage_list))

    num_list = np.zeros(101)  # creat a list with 100 zeros, from 0 to 100%, 1% range
    for percentage in coverage_list:
        for i in range(0, 101):
            if i <= percentage < i + 1:
                num_list[i] += 1

    df2 = pd.DataFrame(dict(coverage=np.arange(0,101,1), number_of_proteins=num_list))

    identified_num_list = np.zeros(101)
    for percentage in identified_coverage_array:
        for i in range(0, 101):
            if i <= percentage < i + 1:
                identified_num_list[i] += 1

    df3 = pd.DataFrame(dict(coverage=np.arange(0,101,1), number_of_identified_proteins=identified_num_list))


    return df1, df2, df3

# the following function generates a dataframe that has protein_length and coverage as two columns
def protein_length_vs_coverage(len_list, coverage_list):
    coverage_array = np.array(coverage_list)
    len_array = np.array(len_list)
    zero_index = np.argwhere(coverage_array == 0)
    identified_array = np.delete(coverage_array, zero_index)
    identified_len_array = np.delete(len_array, zero_index)
    length_interval = range(0, 10010, 10)
    average_coverage_array = np.array([])
    for i in range(len(length_interval) - 1):
        flag = 0
        coverage_sum = 0
        for j in range(len(identified_len_array)):
            if length_interval[i] <= identified_len_array[j] < length_interval[i + 1]:
                flag += 1
                coverage_sum += identified_array[j]
        if flag == 0:
            average_coverage_array = np.append(average_coverage_array, 0)
        else:
            average_coverage = float(coverage_sum) / flag
            average_coverage_array = np.append(average_coverage_array, average_coverage)
    df = pd.DataFrame(dict(protein_length=np.arange(10,10010,10), coverage=average_coverage_array))


    return average_coverage_array, df

def len_cov_trend_line(average_cov_array, interval):
    average_array = np.array([])
    for i in range(len(average_cov_array)-interval):
        average = np.mean(average_cov_array[i:i+interval])
        average_array = np.append(average_array, average)
        df = pd.DataFrame(dict(length=np.arange(len(average_array)), coverage=average_array))

    return df

def length_vs_freq_plot(len_list, len_interval):
    #len_interval = np.arange(0,10010,10)
    freq_array = np.zeros(len(len_interval))
    for leng in len_list:
        for i in range(len(len_interval)-1):
            if len_interval[i]<=leng<len_interval[i+1]:
                freq_array[i] += 1
    df = pd.DataFrame(dict(protein_length=len_interval, frequency=freq_array))
    return df

# the following function use aho-corasick result and pos_ID_dict as parameter to generate a dict with pep sequence as key and uniprotID as value.
def creat_pep_ID_dict(aho_result, pos_ID_dict):
    pep_ID_dict = defaultdict(set)
    for i in aho_result:
        pep_ID_dict[i[2]].add(pos_ID_dict[i[0]])
    return pep_ID_dict

# a dictionary that has Uniport ID as key, identified peptides for that ID as value
def creat_ID_pep_dict(aho_result, pos_ID_dict):
    ID_pep_dict = defaultdict(set)
    for i in aho_result:
        ID_pep_dict[pos_ID_dict[i[0]]].add(i[2])
    return ID_pep_dict

def counting_peptide_in_one_protein(ID_pep_dict):
    ID_pep_num_dict = dict()
    for ID in ID_pep_dict:
        ID_pep_num_dict[ID] = len(ID_pep_dict[ID])
    return ID_pep_num_dict

def counting_unique_pep_in_one_protein(pep_id_dict):
    ID_unique_pep_dict = defaultdict(set) # UniprotID as key, set of unique peptides as value
    ID_number_of_unique_pep_dict = {} # UniprotID as key, number of unique peptides as value
    for pep in pep_id_dict:
        if len(pep_id_dict[pep]) == 1: # unique peptide only has one corresponding uniprot_id
            for id in pep_id_dict[pep]: # pep_id_dict[pep] is set, so iterate it to convert set back to string
                ID_unique_pep_dict[id].add(pep)

    for ID in ID_unique_pep_dict:
        ID_number_of_unique_pep_dict[ID] = len(ID_unique_pep_dict[ID])
    return ID_unique_pep_dict,ID_number_of_unique_pep_dict

def NSAF_calculation(protein_dict, ID_pep_num_dict):
    ID_NSAF_dict = dict()
    sum_si_over_Li = 0
    for ID in ID_pep_num_dict:
        Li = len(protein_dict[ID])
        si = int(ID_pep_num_dict[ID])
        sum_si_over_Li += float(si)/Li

    for ID in ID_pep_num_dict:
        Li = len(protein_dict[ID])
        si = ID_pep_num_dict[ID]
        si_over_Li = float(si)/Li
        NASF = si_over_Li/sum_si_over_Li
        ID_NSAF_dict[ID] = NASF
    return ID_NSAF_dict

def read_reverse_into_fasta(protein_dict, protein_description_dict, reverse_fasta_file): #  the usage of protein description dict is from protein_coverage.py
    with open(reverse_fasta_file, 'a') as file_open:

        for ID in protein_dict:
            reverse_seq = protein_dict[ID][::-1]
            if ('SWISS-PROT' or 'TREMBL' or 'ENSEMBL') not in ID: # read normal ID
                file_open.write('>Reverse|'+ID+'|'+protein_description_dict[ID]+'\n')
                interval = np.arange(0,len(reverse_seq),60)
                interval = np.append(interval, len(reverse_seq))
                for i in range(len(interval)-1):
                    file_open.write(reverse_seq[interval[i]:interval[i+1]]+'\n')
            else:# read contaminants
                file_open.write('>Reverse|'+ID+'|'+protein_description_dict[ID])
                interval = np.arange(0, len(reverse_seq), 60)
                interval = np.append(interval, len(reverse_seq))
                for i in range(len(interval) - 1):
                    file_open.write(reverse_seq[interval[i]:interval[i + 1]] + '\n')

    return reverse_fasta_file

if __name__=='__main__':
    # load information required to do calculations
    pfiles = ['number_line_5_5.p', 'seq_line_5_5.p', 'ID_list_5_5.p', 'seq_list_5_5.p', 'sep_pos_5_5.p', 'aho_result_5_5.p']
    zero_line = ppp.load(open(pfiles[0], 'rb'))
    sep_pos_array = ppp.load(open(pfiles[4], 'rb'))
    uniprot_list = ppp.load(open(pfiles[2], 'rb'))
    seq_list = ppp.load(open(pfiles[3], 'rb'))
    seq_line = ppp.load(open(pfiles[1], 'rb'))
    aho_result = ppp.load(open(pfiles[5], 'rb'))

    pep_ID_dict = creat_pep_ID_dict(aho_result, multiprocessing_naive_algorithym.read_position_ID_into_dict(uniprot_list, seq_list, zero_line))
    ID_number_of_unique_pep_dict = counting_unique_pep_in_one_protein(pep_ID_dict)[1]
    #print ID_number_of_unique_pep_dict
    print (len(zero_line), len(sep_pos_array), len(uniprot_list))

    overall_percentage, identified_coverage_array, coverage_list, protein_coverage_dict, average_coverage, identified_average_coverage, length_list = coverage_calculation(zero_line, sep_pos_array, uniprot_list)
    print ('The overall overage is: %s\nThe average coverage for all of proteins is: %s\nThe average coverage for identified proteins is: %s' \
    % (overall_percentage, average_coverage, identified_average_coverage))
    '''
    for ID in protein_coverage_dict:
        if 20>protein_coverage_dict[ID]:
            print ID, protein_coverage_dict[ID], ID_number_line(zero_line, sep_pos_array, uniprot_list)[ID]
    '''



    ID_ratio_dict = ratio_of_coverage_to_sum(zero_line, sep_pos_array, uniprot_list)
    identified_ID_list = ID_ratio_dict.keys()
    print (len(identified_ID_list))
    ratio_list = ID_ratio_dict.values()
    average_cov_array, df1 = protein_length_vs_coverage(length_list, coverage_list)

    df4 = pd.DataFrame(dict(Uniprot_ID=identified_ID_list, ratio_of_coverage_to_freq_sum=ratio_list))
    #df4.to_excel('protein_seletion.xlsx', sheet_name='sheet1', index=False)
    df5 = pd.DataFrame(dict(protein_length=length_list, sequence_coverage=coverage_list))
    sns.set()

    protein_graph, identified_protein_graph = read_into_dataframe(protein_coverage_dict, identified_coverage_array)[1:]
    sns.lineplot(x='coverage', y='number_of_proteins', data=protein_graph)
    plt.show()

    sns.barplot(x='coverage', y='number_of_identified_proteins', data=identified_protein_graph)
    plt.xlabel('sequence coverage')
    plt.ylabel('number of identified protein isoforms')
    plt.xticks(range(0,105,5), range(0,105,5))
    plt.show()
    '''
    #plt.savefig('coverage_vs_iden_protein_number.png', dpi=400)
    len_intervl = np.arange(0,10010,10)
    df_freq = length_vs_freq_plot(length_list, len_intervl)
    sns.lineplot(x='protein_length', y='frequency', data=df_freq)
    plt.show()
    pd.set_option('display.max_rows', 1000)
    df = len_cov_trend_line(average_cov_array, 10)
    sns.lineplot(x='length', y='coverage', data=df)
    plt.show()
    '''

    plt.ylim((0,100))
    sns.scatterplot(x='protein_length', y='coverage', data=df1)
    plt.xlabel('protein length')
    plt.ylabel('sequence coverage')
    plt.show()



