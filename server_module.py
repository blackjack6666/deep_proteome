import ahocorasick
import numpy as np
from collections import defaultdict
import re


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


def creat_total_seq_line(seq_list):
    seq_line = '|'.join(seq_list)
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


def peptide_phospho_reader(peptide_tsv_file, mod='\d+\.\d+'): # 79.9663 is the delta mass of phosphorylation on STY
    pep_phos_dict = defaultdict()
    with open(peptide_tsv_file) as file_open:
        for i in range(1):
            next(file_open)
        for line in file_open:
            pep_seq = line.split('\t')[0]

            pattern = re.compile('\d+\w{1}\('+str(mod)+'\)')
            regex = re.findall(pattern, line)
            #print (regex)
            for ele in regex:
                if ele != '':
                    pep_phos_dict[pep_seq]=regex
    return pep_phos_dict


def peptide_dataframe(peptide_list, fasta_file, peptide_ptm_dict):
    """
    get sequence coverage, mapped array, ptm information into a dictionary
    :param peptide_list:
    :param fasta_file:
    :param peptide_ptm_dict: a dictionary of identified peptides with list of ptms
    :return:
    """
    # get some variables ready
    prot_dict = fasta_reader(fasta_file)
    ID_list, seq_list = extract_UNID_and_seq(prot_dict)
    seq_line = creat_total_seq_line(seq_list)
    separtor_pos_array = separator_pos(seq_line)
    zero_line = zero_line_for_seq(seq_line)
    pos_id_dict = read_position_ID_into_dict(ID_list,seq_list,zero_line)

    # peptide mapping
    A = automaton_trie(peptide_list)
    aho_result = automaton_matching(A,seq_line)
    id_pep_dict = creat_ID_pep_dict(aho_result,pos_id_dict)

    # add to array
    for pos in aho_result:
        zero_line[pos[0]:pos[1] + 1] += 1

    # calculations
    result_dict = {}
    for i in range(len(separtor_pos_array)-1): # iterate from the left of zeroline all the way to the right

        mapped_array_one_prot = zero_line[separtor_pos_array[i]+1:separtor_pos_array[i+1]]
        non_zero_for_one_protein = np.count_nonzero(mapped_array_one_prot)
        if non_zero_for_one_protein != 0:  # filter unseen protein
            protein_id = ID_list[i]
            total_len_for_one_protein = len(mapped_array_one_prot)
            coverage_for_one_protein = float(non_zero_for_one_protein) / total_len_for_one_protein * 100
            ptm_dict_one_prot = {each:peptide_ptm_dict[each] for each in id_pep_dict[protein_id] if each in peptide_ptm_dict}
            result_dict[protein_id] = (coverage_for_one_protein,
                                       mapped_array_one_prot,
                                       ptm_dict_one_prot)
        else:
            continue

    return result_dict



if __name__ == '__main__':
    from tsv_reader import peptide_counting
    tsv_file = 'D:/data/Naba_deep_matrisome/01102021/163_3B18/peptide.tsv'
    pep_ptm_dict = peptide_phospho_reader(tsv_file)
    peptide_list = peptide_counting(tsv_file)
    fasta_path = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
    result_dict = peptide_dataframe(peptide_list,fasta_path,pep_ptm_dict)
    print (result_dict['Q9EQH3'])