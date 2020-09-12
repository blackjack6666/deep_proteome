'''
import numpy as np
import re

zero_line = np.zeros(100)

sep_pos_list = [5,9,15,23,35,68]

zero_line[10:20] += 1
zero_line[28:36] += 1
print zero_line
for i in range(len(sep_pos_list)-1):
    print np.count_nonzero(zero_line[sep_pos_list[i]+1:sep_pos_list[i+1]])

seq_line = 'adsfasdgsadfasdf|asdgsadgasd|fasdfasdf|asdfwaef|'
sep_pos_list = [m.start() for m in re.finditer('\|', seq_line)]
print sep_pos_list
'''
'''
import numpy as np
import pickle as p
A = np.array([1,2,31,2123,12,31,24,12,3,12,31,23,12,3,12,3])
p.dump(A, open('a.p', 'wb'))
B = [12,31,2,41,2,31,4,12,3,12,4,1,23,1,24,124,12,3,12,4,1,1,23,13]
B = np.array(B)
'''
'''
import numpy as np
a = np.zeros(100)
a[40:60] += 1
index = np.argwhere(a > 0)
b = np.delete(a, index)
print len(b)
'''

from protein_coverage import read_fasta_info_dict2

#from test4 import read_peptide
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

def creat_total_seq_line(seq_list):
    seq_line = '|'.join(seq_list)
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


if __name__ == "__main__":

    filename = 'C:/uic/lab/data/xinhao_data1/uniprot-proteome_UP000005640.fasta'
    test_filename = 'C:/uic/lab/data/TEST/test_fasta.txt'
    protein_dict = read_fasta_into_dict(filename)[0]
    #print protein_dict, protein_dict.values()
    uniprot_ID_list, seq_list = extract_UNID_and_seq(protein_dict)
    seq_line = creat_total_seq_line(seq_list)
    zero_line = zero_line_for_seq(seq_line)
    seq_line_ID_dict = read_position_ID_into_dict(uniprot_ID_list, seq_list, zero_line)
    #ppp.dump(seq_line_ID_dict, open('ID_position_dict.P'), protocol=-1) # read it into a pickle file
    path = 'C:/uic/lab/data/xinhao_data1/'
    test_path = 'C:/uic/lab/data/TEST/'
    dtafiles = glob.glob(path+'*.dta')
    start = time.clock()
    peptide_list = read_peptide(dtafiles)
    peptide_list_set = set(peptide_list)
    peptide_list_unique = list(peptide_list_set)
    print (time.clock()-start, len(peptide_list), len(peptide_list_unique))

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