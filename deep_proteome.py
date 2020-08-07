import tsv_reader
import time
import pandas as pd
import numpy as np
from protein_coverage import fasta_reader, read_fasta_into_dict
from MS_tools_in_silico_digestion import peptide_generator


def coverage_calculating(list_of_peptide, protein_seq):
    freq_array = np.zeros(len(protein_seq))

    for peptide in list_of_peptide:
        if peptide in protein_seq:
            start_pos = protein_seq.find(peptide)
            end_pos = start_pos + len(peptide) - 1
            freq_array[start_pos:end_pos + 1] += 1

    non_zero = np.count_nonzero(freq_array)
    coverage = float(non_zero) / len(freq_array) * 100

    return coverage


def coverage_from_tsv(protein_tsv_path):
    coverage_dict = {}  # protein ID is key, coverage is value
    with open(protein_tsv_path, 'r') as f:
        for i in range(1):
            next(f)
        for line in f:
            line_split = line.split('\t')
            protein_ID = line_split[3]
            cov = float(line_split[7])
            coverage_dict[protein_ID] = cov
    return coverage_dict


def KR_ratio_calculator(protein_seq):
    KR_counting = 0
    for aa in protein_seq:
        if aa == 'K' or aa == 'R':
            KR_counting += 1
        else:
            continue
    return float(KR_counting)/len(protein_seq)

def trypsin_chymotrypsin_seq_cov(trypsin_protein_path, chymotrypsin_protein_path, trypsin_peptide_path,chymotrypsin_peptide_path, fasta_path,out_excel):

    trypsin_protein_list = tsv_reader.protein_tsv_reader(trypsin_protein_path)
    chymo_protein_list = tsv_reader.protein_tsv_reader(chymotrypsin_protein_path)
    intersection = [p for p in trypsin_protein_list if p in chymo_protein_list and not p.startswith('CT')]  # list of protein Uniprot IDs
    print (intersection)

    trypsin_pep_list = tsv_reader.peptide_counting(trypsin_peptide_path)
    chymo_pep_list = tsv_reader.peptide_counting(chymotrypsin_peptide_path)
    combined_pep_list = list(set(trypsin_pep_list+chymo_pep_list))

    info_list = []

    protein_dict = fasta_reader(fasta_path)
    for protein in intersection:
        sequence = protein_dict[protein]

        # calculating trypsin coverage
        trypsin_coverage = coverage_calculating(trypsin_pep_list,sequence)

        # calculating chymotrypsin coverage
        chymotrypsin_coverage = coverage_calculating(chymo_pep_list,sequence)

        # combined coverage
        combined_coverage = coverage_calculating(combined_pep_list,sequence)

        # getting trypsin inslico coverage
        inslico_coverage = coverage_calculating(peptide_generator((protein,sequence),'trypsin'),sequence)

        info_list.append((protein,trypsin_coverage,chymotrypsin_coverage,combined_coverage,inslico_coverage))

    protein_list, trypsin_cov_list, chymo_cov_list,comb_cov_list, inslico_cov_list= [list(t) for t in zip(*info_list)]

    df = pd.DataFrame(dict(Uniprot_ID=protein_list,
                           trypsin_coverage=trypsin_cov_list,
                           chymotrypsin_coverage=chymo_cov_list,
                           combined_coverage=comb_cov_list,
                           inslico_coverage=inslico_cov_list),
                      columns=['Uniprot_ID','trypsin_coverage','chymotrypsin_coverage','combined_coverage','inslico_coverage'])
    df.to_excel(out_excel)

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

def output_differential_digestion_coverage(list_of_protein_tsv, list_of_peptide_tsv, protein_fasta):

    """
    list of protein tsv and list of peptide tsv are corresponding with each other, aho-corasick implementation
    :param list_of_proteintsv:
    :param list_of_peptide_tsv:
    :param protein_fasta:
    :param out_excel:
    :return:
    """
    import aho_corasick
    from collections import defaultdict
    from calculations_and_plot import coverage_calculation
    import multiprocessing_naive_algorithym
    from itertools import combinations


    all_file_protein_list = [tsv_reader.protein_tsv_reader(protein_path) for protein_path in list_of_protein_tsv]  # 2D array

    #common_protein_list = list(set.intersection(*map(set, all_file_protein_list)))  # get the proteins IDs they all have in common, 1D list

    #common_protein_list = [pro for pro in common_protein_list if not pro.startswith('CT')]
    all_protein_list = []
    for path in list_of_protein_tsv:
        all_protein_list += tsv_reader.protein_tsv_reader(path)
    common_protein_list = [pro for pro in set(all_protein_list) if not pro.startswith('CT')]
    print (len(common_protein_list))
    protein_dict = fasta_reader(protein_fasta)

    #print (len(common_protein_list))
    # get the long sequence line ready for aho implementation
    seq_line = '|'.join([protein_dict[pro] for pro in common_protein_list])
    #print (seq_line)
    sep_pos_array = multiprocessing_naive_algorithym.separator_pos(seq_line)
    #print (sep_pos_array)
    df = pd.DataFrame()
    df['Protein ID'] = common_protein_list

    all_file_peptide_list = [tsv_reader.peptide_counting(peptide_path) for peptide_path in list_of_peptide_tsv]

    # for each peptide list, calculate the protein coverage list
    time_start = time.clock()
    zero_line_list = []  # 2D array with zero_lines of all matching
    for pep_list, pep_path in zip(all_file_peptide_list, list_of_peptide_tsv):

        automaton = aho_corasick.automaton_trie(pep_list)
        zero_line = multiprocessing_naive_algorithym.zero_line_for_seq(seq_line)
        result = aho_corasick.automaton_matching(automaton, seq_line)
        #print (result)
        for pos in result:
            zero_line[pos[0]:pos[1] + 1] += 1

        zero_line_list.append(zero_line)
        coverage_list = coverage_calculation(zero_line, sep_pos_array, common_protein_list)[2]

        df[pep_path.split('/')[-2]] = coverage_list
    print (time.clock()-time_start)

    # calculate the overlap sequence ratio
    combination_of_two_zero_line = list(combinations(zero_line_list,2))  # [([],[]), ([],[]), ([],[])]
    combination_of_two_path_list = list(combinations(list_of_peptide_tsv,2))  # [('tsv1', 'tsv2'), (), ()...]
    flag = 0
    for zero_line_tuple, tsv_comb in zip(combination_of_two_zero_line,combination_of_two_path_list):
        overlap_zero_line = multiprocessing_naive_algorithym.zero_line_for_seq(seq_line)

        overlap_pos = [ind for ind,val in enumerate(zero_line_tuple[0]) if val !=0 and zero_line_tuple[1][ind] != 0]
        for pos in overlap_pos:
            overlap_zero_line[pos] += 1
        overlap_coverage_list = coverage_calculation(overlap_zero_line, sep_pos_array,common_protein_list)[2]
        df[tsv_comb[0].split('/')[-2]+'_'+tsv_comb[1].split('/')[-2]] = overlap_coverage_list
        flag += 1
        print (flag)
    return df


def getting_overlap_sequence(p1, p2, protein_seq):  # p1 and p2 are identified in protein_seq
    """
    get the overlapping sequence between two peptides, ex: p1 from typsin digestion, p2 from chymotrypsin digestion,
    p1 and p2 are identified from the same protein

    :param p1:
    :param p2:
    :return:
    """
    import difflib
    p1_pos, p2_pos = protein_seq.find(p1), protein_seq.find(p2)
    if abs(p1_pos - p2_pos) < len(p1) or abs(p1_pos - p2_pos) < len(p2):  # make sure p1 and p2 are overlapped
        s = difflib.SequenceMatcher(None, p1, p2)
        pos_a, pos_b, size = s.find_longest_match(0, len(p1), 0, len(p2))
        print(pos_a, pos_b, size)
        overlap = p1[pos_a:pos_a + size]
        print(p2_pos, p1_pos, len(p1), len(p2), len(overlap))
        return overlap
    # if (pos_a+size == len(p1) or pos_b+size == len(p2)) and abs(p2_pos-p1_pos) == len(p1)+len(p2)-len(overlap):  #


def find_aho_compare(peptide_list, protein_dict):
    """
    caculate the run time of .find for peptide matching based on increasing number of proteins and fixed number of peptides,
    returns two array, number of proteins per run, and corresponding run time.
    :param peptide_list:
    :param protein_dict:
    :return:
    """
    import time
    import matplotlib.pyplot as plt

    # traditional .find method
    flag = 1
    result_dict = {}
    number_of_proteins = []  # x-axis
    time_runing = []  # y-axis
    protein_ID_list = [ID for ID in protein_dict.keys()]

    while flag <= len(protein_dict):
        number_of_proteins.append(flag)
        print('matching %i proteins' % (flag))
        start = time.clock()
        find_result = [(protein_dict[ID].find(pep), ID, pep) for ID in
                       protein_ID_list[:flag] for pep in peptide_list]  # find implementation
        time_runing.append(time.clock() - start)
        #result_dict[flag] = find_result
        #print(time.clock() - start)
        flag += 1
    return number_of_proteins, time_runing


def aho_run_time(peptide_list, protein_dict):
    '''
    similar as the last function, only implementing ahocorasick for peptide matching
    :param peptide_list:
    :param protein_dict:
    :return:
    '''
    import time
    import aho_corasick
    import matplotlib.pyplot as plt
    # aho_implementation
    aho_result_dict = {}
    number_of_proteins = []  # x-axis
    aho_time_running = []
    protein_seq_list = [seq for seq in protein_dict.values()]
    A = aho_corasick.automaton_trie(peptide_list)
    flag2 = 20000
    while flag2 <= len(protein_seq_list):
        number_of_proteins.append(flag2)
        seq_line = ''.join(protein_seq_list[:flag2])
        print ('aho_matching %s protein' % (flag2))
        aho_start = time.clock()
        aho_result = aho_corasick.automaton_matching(A, seq_line)
        aho_time_running.append(time.clock()-aho_start)
        flag2 += 1

    return number_of_proteins, aho_time_running


def protein_peptide_match(peptide_list, protein_list, fasta_path):
    """
    returns a dictionary with each matched protein as key, a list of matched peptides as value
    :param peptide_list:
    :param protein_list:
    :param fasta_path:
    :return:
    """
    import aho_corasick
    import multiprocessing_naive_algorithym
    from collections import defaultdict
    protein_matched_peptide_dict = defaultdict(list)

    protein_dict = fasta_reader(fasta_path)
    seq_list = [protein_dict[pro] for pro in protein_list]
    seq_line = '|'.join(seq_list)
    # automaton trie building and iterate the sequence line
    A = aho_corasick.automaton_trie(peptide_list)
    aho_result = aho_corasick.automaton_matching(A, seq_line)
    # generate the protein ID for each position on the sequence line
    pos_ID_dict = multiprocessing_naive_algorithym.read_position_ID_into_dict(protein_list,seq_list,seq_line)
    for each in aho_result:
        protein_matched_peptide_dict[pos_ID_dict[each[0]]].append(each[2])

    return protein_matched_peptide_dict


def protein_pep_pos_gen(peptide_list, protein_list, fasta_path):
    '''
    return a dictionary with each protein as key, a list of matched peptides position sets as value
    :param peptide_list:
    :param protein_list:
    :param fasta_path:
    :return:
    '''
    import aho_corasick
    import re
    import multiprocessing_naive_algorithym
    from collections import defaultdict

    protein_pep_pos_dict = defaultdict(list)
    protein_dict = fasta_reader(fasta_path)
    seq_list = [protein_dict[pro] for pro in protein_list]
    seq_line = '|'.join(seq_list)
    # automaton trie building and iterate the sequence line
    A = aho_corasick.automaton_trie(peptide_list)
    aho_result = aho_corasick.automaton_matching(A, seq_line)
    sep_pos_array = np.array([m.start() for m in re.finditer('\|', seq_line)])
    sep_pos_array = np.insert(sep_pos_array, 0, 0)
    # generate the protein ID for each position on the sequence line
    pos_ID_dict = multiprocessing_naive_algorithym.read_position_ID_into_dict(protein_list, seq_list, seq_line)

    for each in aho_result:
        matched_protein_ID = pos_ID_dict[each[0]]
        matched_protein_index = protein_list.index(matched_protein_ID)
        protein_initial_index = sep_pos_array[matched_protein_index]
        protein_pep_pos_dict[matched_protein_ID].append((each[0]-protein_initial_index, each[1]-protein_initial_index,each[2]))
    return protein_pep_pos_dict


if __name__ == "__main__":

    '''
    trypsin_protein_path = 'D:/data/Deep_proteome/10_29/trypsin_result/seq_grade_try/protein.tsv'
    chymo_protein_path = 'D:/data/Deep_proteome/10_29/chymotrypsin/chymotrypsin/protein.tsv'
    trypsin_peptide_path = 'D:/data/Deep_proteome/10_29/trypsin_result/seq_grade_try/peptide.tsv'
    chymo_peptide_path = 'D:/data/Deep_proteome/10_29/chymotrypsin/chymotrypsin/peptide.tsv'
    fasta_path = 'D:/data/proteome_fasta/uniprot-proteome_UP000005640.fasta'
    trypsin_chymotrypsin_seq_cov(trypsin_protein_path,chymo_protein_path,trypsin_peptide_path,chymo_peptide_path,fasta_path,'2019_10_29.xlsx')
    
    import matplotlib.pyplot as plt
    pep_tsv = '/home/xshao/lab/wang/10_04_19_second_search/SC/peptide.tsv'
    protein_tsv_path = '/home/xshao/lab/wang/10_04_19_second_search/SC_phos/protein.tsv'
    fasta_path = '/home/xshao/lab/fasta_files/proteome_fasta/UniProt_Mouse_04-17-2019_reversed.fasta'
    from protein_coverage import fasta_reader
    import tsv_reader
    peptide_list = tsv_reader.peptide_counting(pep_tsv)
    total_protein_dict = fasta_reader(fasta_path)

    protein_dict = {p:total_protein_dict[p] for p in [protein for protein in total_protein_dict.keys()][:50000]}

    #protein_dict_from_tsv = {p:total_protein_dict[p] for p in tsv_reader.protein_tsv_reader(protein_tsv_path)}
    #number_of_proteins1, find_time_running = find_aho_compare(peptide_list, protein_dict)
    number_of_proteins, aho_time_running = aho_run_time(peptide_list,protein_dict)
    plt.plot(number_of_proteins, aho_time_running, 'b.')
    plt.xlabel('number of proteins')
    plt.ylabel('time processed')
    plt.title('ahocorasick run time vs number of proteins')
    plt.show()
    '''
    fasta_path = '/home/xshao/lab/fasta_files/proteome_fasta/uniprot-proteome_UP000005640.fasta'
    tryp_5min = '/home/xshao/lab/Deep_proteome/11_14/tryp_5min/protein.tsv'
    tryp_10min = '/home/xshao/lab/Deep_proteome/11_14/tryp_10min/protein.tsv'
    tryp_9h = '/home/xshao/lab/Deep_proteome/11_14/tryp_9h/protein.tsv'
    tryp_19h = '/home/xshao/lab/Deep_proteome/11_14/tryp_19h/protein.tsv'
    tryp_24h = '/home/xshao/lab/Deep_proteome/11_14/tryp_24h/protein.tsv'

    CT_5min = '/home/xshao/lab/Deep_proteome/11_14/CT_5min/protein.tsv'
    CT_10min = '/home/xshao/lab/Deep_proteome/11_14/CT_10min/protein.tsv'
    CT_9h = '/home/xshao/lab/Deep_proteome/11_14/CT_9h/protein.tsv'
    CT_19h = '/home/xshao/lab/Deep_proteome/11_14/CT_19h/protein.tsv'
    CT_24h = '/home/xshao/lab/Deep_proteome/11_14/CT_24h/protein.tsv'

    tryp_5min_pep = '/home/xshao/lab/Deep_proteome/11_14/tryp_5min/peptide.tsv'
    tryp_10min_pep = '/home/xshao/lab/Deep_proteome/11_14/tryp_10min/peptide.tsv'
    tryp_9h_pep = '/home/xshao/lab/Deep_proteome/11_14/tryp_9h/peptide.tsv'
    tryp_19h_pep = '/home/xshao/lab/Deep_proteome/11_14/tryp_19h/peptide.tsv'
    tryp_24h_pep = '/home/xshao/lab/Deep_proteome/11_14/tryp_24h/peptide.tsv'

    CT_5min_pep = '/home/xshao/lab/Deep_proteome/11_14/CT_5min/peptide.tsv'
    CT_10min_pep = '/home/xshao/lab/Deep_proteome/11_14/CT_10min/peptide.tsv'
    CT_9h_pep = '/home/xshao/lab/Deep_proteome/11_14/CT_9h/peptide.tsv'
    CT_19h_pep = '/home/xshao/lab/Deep_proteome/11_14/CT_19h/peptide.tsv'
    CT_24h_pep = '/home/xshao/lab/Deep_proteome/11_14/CT_24h/peptide.tsv'

    list_of_protein_path = [tryp_5min, tryp_10min, tryp_9h,tryp_19h, tryp_24h,CT_5min, CT_10min, CT_9h, CT_19h, CT_24h]
    list_of_peptide_path = [tryp_5min_pep, tryp_10min_pep, tryp_9h_pep, tryp_19h_pep, tryp_24h_pep,
                            CT_5min_pep, CT_10min_pep, CT_9h_pep, CT_19h_pep, CT_24h_pep]

    df = output_differential_digestion_coverage(list_of_protein_path,list_of_peptide_path,fasta_path)
    df.to_excel('eee1.xlsx')