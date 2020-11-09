"""
mapping identified peptides to cleaved and missed cleaved sites without using cleavage rule
"""

from collections import defaultdict


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


def polymer_gen(protein_seq, aa_loc):
    """
    generate 31mer based on peptide location, peptide location is at center of 31 mer
    :param protein_seq:
    :param aa_loc:
    :return:
    """
    if aa_loc - 15 >= 0 and aa_loc + 15 < len(protein_seq):  # 31mer is within the protein sequence
        polymer = protein_seq[aa_loc - 15:aa_loc + 16]

    elif aa_loc - 15 < 0 and aa_loc + 15 < len(protein_seq):  # adding Zs to make up the N-terminal
        polymer = 'Z' * (15 - aa_loc) + protein_seq[:aa_loc + 16]

    elif aa_loc - 15 < 0 and aa_loc + 15 >= len(protein_seq):  # adding Zs to both N and C terminal
        polymer = 'Z' * (15 - aa_loc) + protein_seq + 'Z' * (16 - (len(protein_seq) - aa_loc))

    elif aa_loc - 15 >= 0 and aa_loc + 15 >= len(protein_seq):  # adding Zs to make up C-terminal
        polymer = protein_seq[aa_loc - 15:] + 'Z' * (16 - (len(protein_seq) - aa_loc))

    else:
        print (protein_seq,aa_loc)
        raise ValueError('No polymer being generated')
    return polymer


def cleavage_map(id_pep_dict,proteome_dict):
    """
    map identified peptide to cleavage/missed cleaved sites
    :param id_pep_dict:
    :param proteome_dict:
    :return:
    """
    import time
    protein_polymer_sc_dict = {}
    print("mapping peptides to 31mer...")
    start = time.time()
    for prot in id_pep_dict:
        # print (prot)
        scn, scc, scm = defaultdict(int), defaultdict(int), defaultdict(int)
        pep_set = id_pep_dict[prot]
        seq = proteome_dict[prot]
        for pep in pep_set:
            pep_loc = seq.find(pep)
            pep_range = range(pep_loc,pep_loc+len(pep))
            pep_start_31mer,pep_end_31mer = polymer_gen(seq,pep_loc-1),polymer_gen(seq,pep_range[-1])
            scc[pep_start_31mer]+=1
            scn[pep_end_31mer]+=1
            # print (pep_start_31mer,pep_end_31mer)
            # exclude two cleavage sites, count rest aa as missed cleavage
            pep_miss_31mer_list = [polymer_gen(seq,each) for each in pep_range[:-1]]
            for each_missmer in pep_miss_31mer_list:

                scm[each_missmer]+=1
        protein_polymer_sc_dict[prot]=(scn,scc,scm)
    print("mapping end, %fs" % (time.time() - start))

    return protein_polymer_sc_dict


def label_cleavage(protein_polymer_sc_dict):
    """
    the cleavage site should be labeled as 1 if SCn or SCc was at least 1 and SCm was zero.
    labeled as 0 if both SCn and SCc were zero and SCm was at least 1.
    :param protein_polymer_sc_dict: returned from cleavage_map
    :return:
    """
    polymer_label_dict = {}  # {polymer:1 or 0}
    uncertain_polymer_no = 0
    protein_poly_dict = defaultdict(set)
    for prot in protein_polymer_sc_dict:
        scn_dict, scc_dict, scm_dict = protein_polymer_sc_dict[prot]
        total_polymer = set([k for k in scn_dict]+[k for k in scc_dict]+[k for k in scm_dict])
        for polymer in total_polymer:
            scn, scc, scm = protein_polymer_sc_dict[prot][0][polymer], \
                            protein_polymer_sc_dict[prot][1][polymer], \
                            protein_polymer_sc_dict[prot][2][polymer]
            if (scm == 0 and scn >= 1) or (scm == 0 and scc >= 1):
                polymer_label_dict[polymer] = 1
                protein_poly_dict[prot].add(polymer)
            elif scn == 0 and scc == 0 and scm >= 1:
                polymer_label_dict[polymer] = 0
                protein_poly_dict[prot].add(polymer)
            else:
                uncertain_polymer_no += 1
                protein_poly_dict[prot].add(polymer)
                # print(each_polymer)
                # continue
    return polymer_label_dict, protein_poly_dict, uncertain_polymer_no


if __name__=='__main__':
    from tsv_reader import peptide_counting, psm_reader, protein_tsv_reader_no_contam
    import pickle as ppp
    from collections import Counter
    from protein_coverage import fasta_reader2

    # protein_tsv_path = "D:/data/deep_proteome/20200915_tryp_37C_1440min/protein.tsv"
    peptide_tsv_path = "D:/data/non_specific_search/peptide.tsv"
    # psm_tsv_path = "D:/data/deep_proteome/20200915_tryp_37C_1440min/psm.tsv"

    fasta_path = 'D:/data/proteome_fasta/uniprot-proteome_UP000000558_ecoli.fasta'
    proteome_dict = fasta_reader2(fasta_path)
    # protein_seq = proteome_dict['P22234']

    # print (protein_seq[410])
    pep_list = peptide_counting(peptide_tsv_path)
    id_pep_dict,protein_list = protein_id_peplist_dict_getter(proteome_dict, pep_list)
    protein_polymer_sc_dict = cleavage_map(id_pep_dict,proteome_dict)
    polymer_label_dict, protein_poly_dict, uncertain_polymer_no = label_cleavage(protein_polymer_sc_dict)
    print('number of proteins with polymers reported: %i' % len(protein_poly_dict))
    print(Counter([v for v in polymer_label_dict.values()]), len(polymer_label_dict))
    print('uncertain ploymer number: %i' % uncertain_polymer_no)
    ppp.dump(polymer_label_dict, open('D:/data/non_specific_search/ecoli_non_specific_search_poly_dict.p', 'wb'))
