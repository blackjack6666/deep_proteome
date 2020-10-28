"""
label in-silico cleavage site as 0 or 1 based on spectral counts of cleavage site N, C terminal and spectral counts of
peptide that covers cleavage site, data preparation for downstream statistical training.
"""

from MS_tools_parameters import expasy_rules
import re
from protein_coverage import fasta_reader2
from collections import defaultdict


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


def miss_clea_adjacent_loc(protein_miss_clea_loc_dict):
    """
    map the index of cleavage site adjacent amino acids to corresponding missed cleavage
    :param protein_miss_clea_loc_dict: {protein_id: missed cleavage location list}
    :return: dictionary of dictionary, {each_protein: {each_aa_loc:set(mapped miss_cleav_index)}}
    """
    loc_dict = {}
    for each_prot in protein_miss_clea_loc_dict:
        # print (each_prot)
        miss_cleav_list = protein_miss_clea_loc_dict[each_prot]
        if len(miss_cleav_list)>=1:
            # print (miss_cleav_list)
            aa_miss_map_dict = defaultdict(set)

            first_iter_range = range(0,miss_cleav_list[0]+15) # first missed cleavage
            for each_loc in first_iter_range:
                aa_miss_map_dict[each_loc].add(miss_cleav_list[0])

            for each_miss_cleav_loc in miss_cleav_list[1:]:
                iter_range = range(each_miss_cleav_loc-15,each_miss_cleav_loc+15) if each_miss_cleav_loc-15>0 else \
                    range(0,each_miss_cleav_loc+15)
                for each_loc in iter_range:
                    aa_miss_map_dict[each_loc].add(each_miss_cleav_loc)
            loc_dict[each_prot]=aa_miss_map_dict

    return loc_dict


def ploymer_miss_cleav_charc(cleavage_site_dict, proteome_seq_dict):
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

            elif cleav_index-15 < 0 and cleav_index + 15 <= len(seq): # adding Zs to make up the N-terminal
                polymer = 'Z'*(15-cleav_index)+seq[:cleav_index+16]

            elif cleav_index-15 < 0 and cleav_index + 15 > len(seq): # adding Zs to both N and C terminal
                polymer = 'Z'*(15-cleav_index)+seq+'Z'*(16-(len(seq)-cleav_index))

            elif cleav_index-15 >=0 and cleav_index+15 > len(seq): # adding Zs to make up C-terminal
                polymer = seq[cleav_index-15:]+'Z'*(16-(len(seq)-cleav_index))

            cleav_index_polymer_dict[cleav_index] = polymer

        total_dict[prot]=cleav_index_polymer_dict

    return total_dict


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


def spec_count_polymer(cleavage_site_dict,
                        proteome_seq_dict,
                        id_pep_dict,
                        aa_index_cleav_map_dict,
                        cleavage_polymer_dict,
                        psm_dict):
    """
    -----
    stats on missed cleavage polymers with SCn, SCc, SCm, outputs a dictionary {protein_id:(SCn_dict,SCc_dict,SCm_dict)}
    SCn: Spectral count of peps in missed cleavage n terminal
    SCc: spectral counts of peps in missed cleavage c terminal
    SCm: spectral counts of peps that contains missed cleavage site.
    -----
    :param proteome_seq_dict:
    :param id_pep_dict:
    :param aa_index_cleav_map_dict:
    :return:
    """
    protein_polymer_sc_dict = {}
    for prot_id in id_pep_dict:
        if prot_id in aa_index_cleav_map_dict:
            seq = proteome_seq_dict[prot_id]  # protein sequence
            #cleavage_site_list = cleavage_site_dict[prot_id]  # the cleavage index list for current protein

            SCn_count_dict, SCc_count_dict, SCm_count_dict = defaultdict(int),defaultdict(int),defaultdict(int)

            for pep in id_pep_dict[prot_id]:  # for each matched peptide in one protein
                spec_count = psm_dict[pep]
                pep_ind = seq.find(pep)  # find location of the peptide
                cleav_ind_set = aa_index_cleav_map_dict[prot_id][pep_ind]  # get cleavage sites that are near the peptide

                for each_cleav_ind in cleav_ind_set:
                    polymer = cleavage_polymer_dict[prot_id][each_cleav_ind]
                    SCn, SCc, SCm = 0,0,0  # or array = np.zeros(3)
                    if pep_ind > each_cleav_ind:  # peptide is identified at C terminal of cleavage site
                        SCc += spec_count
                    elif pep_ind <= each_cleav_ind:  # peptide is identified at N terminal of cleavage site
                        if pep_ind+len(pep) <= each_cleav_ind+1:
                            SCn += spec_count
                        else:
                            SCm += spec_count

                    # add spec count to polymer
                    SCn_count_dict[polymer] += SCn
                    SCc_count_dict[polymer] += SCc
                    SCm_count_dict[polymer] += SCm
            protein_polymer_sc_dict[prot_id]=(SCn_count_dict,SCc_count_dict,SCm_count_dict)
    return protein_polymer_sc_dict


def cleavage_site_label(protein_polymer_sc_dict):
    """
    -----
    the cleavage site should be labeled as 1 if SCn or SCc was at least 1 and SCm was zero.
    labeled as 0 if both SCn and SCc were zero and SCm was at least 1.
    -----
    :param protein_polymer_sc_dict: returned by last function spec_count_polymer,
    {protein_id:(SCn_dict,SCc_dict,SCm_dict)}
    :return:
    """
    polymer_label_dict = {}  # {polymer:1 or 0}
    protein_poly_dict = defaultdict(set)
    for prot in protein_polymer_sc_dict:
        dict_tuple = protein_polymer_sc_dict[prot]

        for each_polymer in dict_tuple[0]:
            SCn,SCc,SCm = dict_tuple[0][each_polymer],dict_tuple[1][each_polymer],dict_tuple[2][each_polymer]
            if SCm == 0 and SCn >= 1 and SCc >= 1:
                polymer_label_dict[each_polymer] = 1
                protein_poly_dict[prot].add(each_polymer)
            elif SCn == 0 and SCc == 0 and SCm >= 1:
                polymer_label_dict[each_polymer] = 0
                protein_poly_dict[prot].add(each_polymer)
            else:
                # print(each_polymer)
                continue
    return polymer_label_dict, protein_poly_dict


if __name__ == '__main__':
    from tsv_reader import peptide_counting,psm_reader,protein_tsv_reader_no_contam
    import pickle as ppp
    protein_tsv_path = "D:/data/deep_proteome/20200915_ct_37C_240min/protein.tsv"
    peptide_tsv_path = "D:/data/deep_proteome/20200915_ct_37C_240min/peptide.tsv"
    psm_tsv_path = "D:/data/deep_proteome/20200915_ct_37C_240min/psm.tsv"

    fasta_path = 'D:/data/proteome_fasta/uniprot-proteome_UP000005640.fasta'
    proteome_dict = fasta_reader2(fasta_path)
    pep_list = peptide_counting(peptide_tsv_path)
    id_pep_dict,protein_list = protein_id_peplist_dict_getter(proteome_dict, pep_list)
    print (len(id_pep_dict))
    # protein_list = ['Q6QHF9-11']
    # protein_list = protein_tsv_reader_no_contam(protein_tsv_path)

    # print (proteome_dict['Q6QHF9-11'])
    protein_miss_clea_loc_dict = cleavage_site_identify(protein_list,proteome_dict,'chymotrypsin high specificity')

    # print (protein_miss_clea_loc_dict['Q6QHF9-11'])

    aa_miss_cleavage_dict = miss_clea_adjacent_loc(protein_miss_clea_loc_dict)
    # print (miss_clea_adjacent_loc(protein_miss_clea_loc_dict)['Q6QHF9-11'])
    polymers_dict = ploymer_miss_cleav_charc(protein_miss_clea_loc_dict,proteome_dict)
    #print (polymers_dict)

    # id_pep_dict = {'Q6QHF9-11':{'MESTGSVGEAPGGGHGPR', 'RGPHPLGALLR','GGGGRALDPWALPG'}}

    # psm_dict = {'MESTGSVGEAPGGGHGPR':3,'RGPHPLGALLR':5,'GGGGRALDPWALPG':2}
    psm_dict = psm_reader(psm_tsv_path)[0]
    protein_polymer_sc_dict = spec_count_polymer(protein_miss_clea_loc_dict,proteome_dict,id_pep_dict,aa_miss_cleavage_dict,polymers_dict,psm_dict)
    cleavage_site_label_dict,protein_poly_dict = cleavage_site_label(protein_polymer_sc_dict)
    print (cleavage_site_label_dict)
    ppp.dump(cleavage_site_label_dict,open('ct_37C_240min_cleavage_label.p','wb'))
    print (len(cleavage_site_label_dict), len(protein_poly_dict))
