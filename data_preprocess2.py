"""
mapping peptide to in-silico cleavage site
"""
from MS_tools_parameters import expasy_rules
import re
from protein_coverage import fasta_reader2
from collections import defaultdict
import numpy as np


def cleavage_site_identify(protein_list,proteome_seq_dict,enzyme):
    """
    return a dictionary with each protein id in protein list as key, a list of insilico cleavage site index as value
    :param protein_list: a list of protein ids
    :param proteome_seq_dict: whole proteome seq dictionary
    :param enzyme: protease enzyme
    :return:
    """
    regex_cut_rule = expasy_rules[enzyme]
    return {prot:[m.end()-1 for m in re.finditer(regex_cut_rule, proteome_seq_dict[prot])] for prot in protein_list}


def protein_id_peplist_dict_getter(proteome_dict,peptide_list):
    """
    get a dictionary with protein id as key, a set of identified peps as value.
    :param proteome_dict:
    :param peptide_list:
    :return:
    """
    import multiprocessing_naive_algorithym
    import aho_corasick

    ID_list, seq_list = multiprocessing_naive_algorithym.extract_UNID_and_seq(proteome_dict)
    seq_line = multiprocessing_naive_algorithym.creat_total_seq_line(seq_list)
    automaton = aho_corasick.automaton_trie(peptide_list)
    aho_result = aho_corasick.automaton_matching(automaton,seq_line)
    pos_id_dict = multiprocessing_naive_algorithym.read_position_ID_into_dict(ID_list,seq_list,seq_line)
    id_pep_dict = multiprocessing_naive_algorithym.creat_ID_pep_dict(aho_result,pos_id_dict)

    return id_pep_dict, [prot for prot in id_pep_dict]


def polymer_miss_cleav_charc(cleavage_site_dict, proteome_seq_dict):
    """
    output a dictionary of dictionary, {protein:{int(cleavage site index): 'polymer that includes cleavage site'}}
    :param cleavage_site_dict: returned from cleavage_site_identify
    :param proteome_seq_dict:
    :return:
    """
    total_dict = {}
    for prot in cleavage_site_dict:
        seq = proteome_seq_dict[prot]
        cleav_index_polymer_dict= {}
        for cleav_index in cleavage_site_dict[prot]:

            if cleav_index-15 >=0 and cleav_index +15 < len(seq):  # 31mer is within the protein sequence
                polymer = seq[cleav_index-15:cleav_index+16]

            elif cleav_index-15 < 0 and cleav_index + 15 < len(seq): # adding Zs to make up the N-terminal
                polymer = 'Z'*(15-cleav_index)+seq[:cleav_index+16]

            elif cleav_index-15 < 0 and cleav_index + 15 >= len(seq): # adding Zs to both N and C terminal
                polymer = 'Z'*(15-cleav_index)+seq+'Z'*(16-(len(seq)-cleav_index))

            elif cleav_index-15 >=0 and cleav_index+15 > len(seq): # adding Zs to make up C-terminal
                polymer = seq[cleav_index-15:]+'Z'*(16-(len(seq)-cleav_index))


            cleav_index_polymer_dict[cleav_index] = polymer

        total_dict[prot]=cleav_index_polymer_dict

    return total_dict


def protein_polymer_convert(protein_cleav_polymer_dict):
    """
    convert {protein:{cleavage site1:polymer1, cleavage site2: polymer2}} into {protein:{polymer1,polymer2}}
    :param polymer_cleav_polymer_dict: returned by polymer_miss_cleav_charc
    :return:
    """
    protein_polymer_dict = defaultdict(set)
    for prot in protein_cleav_polymer_dict:
        for each_cleav in protein_cleav_polymer_dict[prot]:
            protein_polymer_dict[prot].add(protein_cleav_polymer_dict[prot][each_cleav])
    #print (protein_polymer_dict)
    return protein_polymer_dict


def map_to_cleavage(id_pep_dict,
                    protein_cleavage_dict,
                    protein_polymer_dict,
                    proteome_dict,
                    psm_dict):
    """
    map identified peptide to cleavage sites
    :param id_pep_dict:
    :param protein_cleavage_dict:
    :param
    :return:
    """
    import time
    protein_polymer_sc_dict = {}
    print("mapping peptides to 31mer...")
    start = time.time()
    for prot in id_pep_dict:
        # print (prot)
        scn,scc,scm = defaultdict(int),defaultdict(int),defaultdict(int)
        if prot in protein_cleavage_dict:
            pep_set = id_pep_dict[prot]
            cleavage_array = np.array(protein_cleavage_dict[prot])
            cleav_polymer_dict = protein_polymer_dict[prot]
            seq = proteome_dict[prot]
            for pep in pep_set:
                pep_loc = seq.find(pep)
                pep_end_loc = pep_loc+len(pep)
                # print (pep,pep_loc,pep_end_loc)
                # get missed cleavage location within peptide
                missed_cleav_array = cleavage_array[np.where((cleavage_array>=pep_loc)&(cleavage_array<pep_end_loc-1))]
                # print (missed_cleav_array)
                if pep_end_loc-1 not in cleavage_array and pep_loc-1 in cleavage_array: # protein end
                    pep_start_polymer = cleav_polymer_dict[pep_loc - 1]
                    scc[pep_start_polymer] += psm_dict[pep]
                elif pep_end_loc-1 in cleavage_array and pep_loc-1 not in cleavage_array: # protein start
                    pep_end_polymer = cleav_polymer_dict[pep_end_loc - 1]
                    scn[pep_end_polymer] += psm_dict[pep]
                elif pep_end_loc-1 not in cleavage_array and pep_loc-1 not in cleavage_array: # ex. Q92614 GIVPLAK
                    continue
                else:

                    pep_start_polymer = cleav_polymer_dict[pep_loc - 1]
                    scc[pep_start_polymer] += psm_dict[pep]
                    pep_end_polymer = cleav_polymer_dict[pep_end_loc - 1]
                    scn[pep_end_polymer] += psm_dict[pep]


                if missed_cleav_array.size >0:
                    for each_miss in missed_cleav_array:
                        # print (cleav_polymer_dict[each_miss])
                        scm[cleav_polymer_dict[each_miss]]+=psm_dict[pep]
        protein_polymer_sc_dict[prot] = (scn,scc,scm)
    print("mapping end, %fs" % (time.time() - start))

    return protein_polymer_sc_dict


def cleavage_site_label(protein_polymer_sc_dict,polymer_dict):
    """
    -----
    the cleavage site should be labeled as 1 if SCn or SCc was at least 1 and SCm was zero.
    labeled as 0 if both SCn and SCc were zero and SCm was at least 1.
    -----
    :param protein_polymer_sc_dict: returned by last function spec_count_polymer,
    {protein_id:(SCn_dict,SCc_dict,SCm_dict)}
    :param polymer_dict: in-silico generated protein-polymer dict, {protein1:{polymer1, polymer2...}}, returned by
    protein_polymer_convert
    :return:
    """
    polymer_label_dict = {}  # {polymer:1 or 0}
    protein_poly_dict = defaultdict(set)
    uncertain_polymer_no = 0
    for prot in polymer_dict:
        for polymer in polymer_dict[prot]:
            scn,scc,scm = protein_polymer_sc_dict[prot][0][polymer],\
                          protein_polymer_sc_dict[prot][1][polymer],\
                          protein_polymer_sc_dict[prot][2][polymer]
            # print (scn,scc,scm)

            if (scm == 0 and scn >= 1) or (scm == 0 and scc >= 1):
                polymer_label_dict[polymer] = 1
                protein_poly_dict[prot].add(polymer)
            elif scn == 0 and scc == 0 and scm >= 1:
                polymer_label_dict[polymer] = 0
                protein_poly_dict[prot].add(polymer)
            else:
                uncertain_polymer_no+=1
                protein_poly_dict[prot].add(polymer)
                # print(each_polymer)
                # continue
    return polymer_label_dict, protein_poly_dict,uncertain_polymer_no

if __name__ == '__main__':
    from tsv_reader import peptide_counting,psm_reader,protein_tsv_reader_no_contam
    import pickle as ppp
    from collections import Counter
    protein_tsv_path = "D:/data/deep_proteome/20200915_tryp_37C_1440min/protein.tsv"
    peptide_tsv_path = "D:/data/deep_proteome/20200915_tryp_37C_1440min/peptide.tsv"
    psm_tsv_path = "D:/data/deep_proteome/20200915_tryp_37C_1440min/psm.tsv"

    fasta_path = 'D:/data/proteome_fasta/uniprot-proteome_UP000005640.fasta'
    proteome_dict = fasta_reader2(fasta_path)
    pep_list = peptide_counting(peptide_tsv_path)
    id_pep_dict,protein_list = protein_id_peplist_dict_getter(proteome_dict, pep_list)
    # id_pep_dict = {'P24539':id_pep_dict['P24539']}
    # print (id_pep_dict)
    # protein_list = ['P24539']
    protein_miss_clea_loc_dict = cleavage_site_identify(protein_list, proteome_dict, 'trypsin')
    # protein_miss_clea_loc_dict = {'P24539':protein_miss_clea_loc_dict['P24539']}
    psm_dict = psm_reader(psm_tsv_path)[0]
    polymers_dict = polymer_miss_cleav_charc(protein_miss_clea_loc_dict, proteome_dict)
    # polymers_dict = {'P24539':polymers_dict['P24539']}


    protein_polymer_sc_dict = map_to_cleavage(id_pep_dict,protein_miss_clea_loc_dict,polymers_dict,proteome_dict,psm_dict)
    cleavage_site_label_dict, protein_poly_dict, uncertain_polymer_no = cleavage_site_label(protein_polymer_sc_dict, protein_polymer_convert(polymers_dict))
    print('number of proteins with polymers reported: %i' % len(protein_poly_dict))
    print(Counter([v for v in cleavage_site_label_dict.values()]), len(cleavage_site_label_dict))
    print('uncertain ploymer number: %i' % uncertain_polymer_no)
    # print (cleavage_site_label_dict)
    ppp.dump(cleavage_site_label_dict, open('D:/data/deep_proteome/pickle_file/20200915_tryp_37C_1440min_new.p', 'wb'))
    print(len(cleavage_site_label_dict), len(protein_poly_dict))
    for each in cleavage_site_label_dict:
        if (len(each)) !=31:
            print (each)