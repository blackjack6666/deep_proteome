from collections import defaultdict
import numpy as np
import glob
import pickle as ppp
import multiprocessing as mp
import time

import re


def extract_UNID_and_seq(protein_dict):
    UNID_list = [key for key in protein_dict.keys()]
    seq_list = [value for value in protein_dict.values()]
    return UNID_list, seq_list

def read_seq_length_into_nparray(seq_list):
    seq_length_array = np.array([])
    for each_seq in seq_list:
        seq_length_array = np.append(seq_length_array, len(each_seq))
    return seq_length_array

def creat_total_seq_line(seq_list, sep=None):
    seq_line = '|'.join(seq_list) if sep == '|' else ''.join(seq_list)
    return seq_line

def zero_line_for_seq(seq_line):
    zero_line = np.zeros(len(seq_line))
    return zero_line

# the following function create a dictionary that read each position in long sequence line as key, corresponding UNIPORT ID as value.
def read_position_ID_into_dict(UNID_list, seq_list, zero_line):
    m = 0
    j = 0
    seq_line_ID_dict = dict()
    for i in range(len(zero_line)):
        if j < len(seq_list[m]):
            seq_line_ID_dict[i] = UNID_list[m]
            j += 1
        else:
            j = 0

            m += 1
    return seq_line_ID_dict


def creat_ID_pep_dict(aho_result, pos_ID_dict):
    ID_pep_dict = defaultdict(set)
    for i in aho_result:
        ID_pep_dict[pos_ID_dict[i[0]]].add(i[2])
    return ID_pep_dict


def creat_pep_ID_dict(aho_result, pos_ID_dict):
    pep_ID_dict = defaultdict(set)
    for i in aho_result:
        pep_ID_dict[i[2]].add(pos_ID_dict[i[0]])
    return pep_ID_dict


def create_unique_id_peptide_dict(pep_id_dict):
    """
    get a dictionary with unique peptides for each protein
    :param pep_id_dict:
    :return:
    """
    unique_id_peptide_dict = defaultdict(set)
    unique_id_peptide_count_dict = defaultdict(int)
    unique_pep_id_dict = {pep:prot for pep in pep_id_dict for prot in pep_id_dict[pep]
                          if len(pep_id_dict[pep])==1}
    for pep in unique_pep_id_dict:
        unique_id_peptide_dict[unique_pep_id_dict[pep]].add(pep)

    for id in unique_id_peptide_dict:
        unique_id_peptide_count_dict[id]=len(unique_id_peptide_dict[id])

    return unique_id_peptide_dict, unique_id_peptide_count_dict



def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def find_all(a_str, sub):  # a generator
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start +=  len(sub)  # use start += 1 to find overlapping matches

def find_pep_start(param):
    seq_line, peptide_list=param
    res_dict={}
    for peptide in peptide_list:
        res_dict[peptide]=[m for m in find_all(seq_line,peptide)]
    return res_dict

# the following function returns a dictionary with each peptide as key and corresponding start position list as value.
def start_end_pos_dict(res_dicts):
    start_pos_dict = {}
    end_pos_dict = {}
    for res_dict in res_dicts:
        for peptide in res_dict:
            start_pos_dict[peptide] = res_dict[peptide]
            end_pos_dict[peptide] = [i + len(peptide) - 1 for i in res_dict[peptide]]
    return start_pos_dict, end_pos_dict

# the following function returns a number_line after matching peptides to long_seq_line.
def adding_numbers_to_zero_line(zero_line, start_dict, end_dict):
    for peptide in start_dict:
        start_list_for_each = start_dict[peptide]
        end_list_for_each = end_dict[peptide]
        for i, j in zip(start_list_for_each, end_list_for_each):
            zero_line[i:j] += 1
    return zero_line

def separator_pos(seq_line):
    sep_pos_array = np.array([m.start() for m in re.finditer('\|', seq_line)])
    sep_pos_array = np.insert(sep_pos_array, 0, 0)
    sep_pos_array = np.append(sep_pos_array, len(seq_line))
    return sep_pos_array


def find_protein_start(id_list,seq_list):
    """
    find the start of each protein on the long seq
    :param id_list:
    :param seq_list:
    :return:
    """
    prot_start_pos_dict = {}
    pos = 0
    for each, seq in zip(id_list, seq_list):
        prot_start_pos_dict[each] = pos
        pos += len(seq)+1

    return prot_start_pos_dict


def map_peptide_index(protein_dict,peptide_list):
    """
    get index of peptide on each protein
    :return:
    """
    from aho_corasick import automaton_trie,automaton_matching
    peptide_start_end_dict = defaultdict(set)

    # load variables
    id_list, seq_list = extract_UNID_and_seq(protein_dict)
    protein_start_dict = find_protein_start(id_list,seq_list)

    seq_line = creat_total_seq_line(seq_list,sep='|')
    pos_id_dict = read_position_ID_into_dict(id_list,seq_list,seq_line)

    # ahocorasick mapping
    aho_result = automaton_matching(automaton_trie(peptide_list),seq_line)

    # find start position for each peptide
    for each in aho_result:
        start, end, pep = each
        protein_id = pos_id_dict[start]
        protein_start_pos = protein_start_dict[protein_id]
        start = start-protein_start_pos  # absolute start of peptide - the start of protein
        end = start+len(pep)-1
        peptide_start_end_dict[pep].add((protein_id,start,end))

    return peptide_start_end_dict

if __name__ == "__main__":
    from protein_coverage import fasta_reader
    from tsv_reader import peptide_counting
    fasta_path = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
    protein_dict = fasta_reader(fasta_path)
    protein_ids = [key for key in protein_dict]
    sub_protein_dict={each:protein_dict[each] for each in protein_ids[:3]}
    print (sub_protein_dict.items())
    pep_tsv = 'D:/data/Naba_deep_matrisome/07232021_secondsearch/GFP_120D/peptide.tsv'
    peptide_list = peptide_counting(pep_tsv)
    # peptide_list = ['ISTVEFNYSGDLLATGDKGGRVVIFQR','PEQPVTLRSVLRDEKGALFRAHARYRADSHGELDLARVPALGGSFSG']
    map_peptide_index(protein_dict,peptide_list)

    """
    #chunk_list = chunks(peptide_list, 10)
    chunk_list = chunks(peptide_list_unique, 10)
    parm_list = [(seq_line, p) for p in chunk_list]
    start = time.clock()
    pool = mp.Pool(processes = mp.cpu_count()-4)
    res_dicts = pool.map(find_pep_start, parm_list)
    pool.close()
    pool.join()

    start_pos_dict, end_pos_dict = start_end_pos_dict(res_dicts)
    zero_line = adding_numbers_to_zero_line(zero_line, start_pos_dict, end_pos_dict)
    print (time.clock()-start)

    sep_pos_array = separator_pos(seq_line)
    """
    # trie implementation
    '''
    zero_line_trie = in_trie(make_trie(peptide_list_unique), seq_line)
    total_seq_len_trie = len(zero_line_trie)-len(sep_pos_array) + 2
    total_non_zero_trie = np.count_nonzero(zero_line_trie)
    overall_percentage_trie = float(total_non_zero_trie)/total_seq_len_trie*100
    print total_non_zero_trie
    zero_line_naive = ppp.load(open('zero_line.p', 'rb'))
    print np.count_nonzero(zero_line_naive)
    '''

    """
    print (time.clock()-start)
    total_seq_len = len(zero_line) - len(sep_pos_array) + 2  # artifically added two separator positions into sep_pos_array so need to plus 2
    total_non_zero = np.count_nonzero(zero_line)
    overall_percentage = float(total_non_zero) / total_seq_len * 100
    print (overall_percentage)
    #print len(uniprot_ID_list), len(sep_pos_array)
    ppp.dump(zero_line, open('zero_line.p', 'wb'), protocol=-1)
    ppp.dump(sep_pos_array, open('separator_pos.p', 'wb'), protocol=-1)
    ppp.dump(uniprot_ID_list, open('uniprotID.p', 'wb'), protocol=-1)
        #ppp.dump(uniprot_ID_list, pf)

    #print protein_dict, seq_line, peptide_list, sep_pos_array, zero_line_trie
    """