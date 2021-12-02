import ahocorasick
import numpy as np
from collections import defaultdict
import re


def peptide_counting(peptide_tsv_file):

    with open(peptide_tsv_file, 'r') as file_open:
        next(file_open)

        return [line.split("\t")[0] for line in file_open]


def fasta_reader(fasta_file_path):

    with open(fasta_file_path, 'r') as file_open:
        file_split = file_open.read().split('\n>')

    return {each.split('\n')[0].split('|')[1]: ''.join(each.split('\n')[1:]) for each in file_split}


def automaton_trie(peptide_list):
    A = ahocorasick.Automaton()
    for idx, peptide in enumerate(peptide_list):
        A.add_word(peptide, (idx, peptide))
    A.make_automaton()
    return A


def automaton_matching(A, seq_line):
    result = []
    for end_idx, (insert_order, original_value) in A.iter(seq_line):
        start_idx = end_idx - len(original_value) + 1
        result.append((start_idx, end_idx, original_value))
        assert seq_line[start_idx:start_idx+len(original_value)] == original_value
    return result


def extract_UNID_and_seq(protein_dict):
    UNID_list = [key for key in protein_dict.keys()]
    seq_list = [value for value in protein_dict.values()]
    return UNID_list, seq_list


def creat_total_seq_line(seq_list, sep=None):
    seq_line = '|'.join(seq_list) if sep == '|' else ''.join(seq_list)
    return seq_line


def zero_line_for_seq(seq_line):
    zero_line = np.zeros(len(seq_line))
    return zero_line


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


def separator_pos(seq_line):
    sep_pos_array = np.array([m.start() for m in re.finditer('\|', seq_line)])
    sep_pos_array = np.insert(sep_pos_array, 0, 0)
    sep_pos_array = np.append(sep_pos_array, len(seq_line))
    return sep_pos_array


def main(fasta_file:str, peptide_list:list):

    # load variables
    prot_cov_dict = {}
    protein_dict = fasta_reader(fasta_file)
    ID_list, seq_list = extract_UNID_and_seq(protein_dict)

    seq_line = creat_total_seq_line(seq_list, sep='|')
    zero_line = zero_line_for_seq(seq_line)
    separtor_pos_array = separator_pos(seq_line)

    # ahocorasick mapping
    aho_result = automaton_matching(automaton_trie(peptide_list),seq_line)

    # calculate coverage for each protein
    for pos in aho_result:
        zero_line[pos[0]:pos[1] + 1] += 1

    for i in range(len(separtor_pos_array)-1): # iterate from the left of zeroline all the way to the right
        non_zero_for_one_protein = np.count_nonzero(zero_line[separtor_pos_array[i]+1:separtor_pos_array[i+1]])
        total_len_for_one_protein = len(zero_line[separtor_pos_array[i] + 1:separtor_pos_array[i + 1]])
        coverage_for_one_protein = float(non_zero_for_one_protein) / total_len_for_one_protein * 100
        prot_cov_dict[ID_list[i]] = coverage_for_one_protein

    # delete unidentified proteins
    iden_prot_cov_dict = {each:prot_cov_dict[each] for each in prot_cov_dict if prot_cov_dict[each] != 0}
    return iden_prot_cov_dict


if __name__ == '__main__':
    peptide_tsv = 'D:/data/native_protein_digestion/11182021/search_result_XS/0060min_XS/peptide.tsv'
    peptide_list = peptide_counting(peptide_tsv)
    fasta_file = ''
    cov_dict = main(fasta_file,peptide_list)
