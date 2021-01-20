import pandas as pd
from multiprocessing_naive_algorithym import extract_UNID_and_seq,creat_total_seq_line,creat_ID_pep_dict,\
    read_position_ID_into_dict
from tsv_reader import psm_reader,peptide_counting
from calculations_and_plot import identified_proteome_cov
import aho_corasick
import os
from protein_coverage import fasta_reader
import matplotlib.pyplot as plt
import statistics


def dash_dataframe(pep_path_list, psm_path_list, protein_dict, ecm_prot_list, ecm_info_dict):
    file_name_number_dict = {pep_tsv.split('/')[-2]: i for pep_tsv, i in zip(pep_path_list, range(len(pep_path_list)))}
    id_list, seq_list = extract_UNID_and_seq(protein_dict)
    seq_line = creat_total_seq_line(seq_list)
    pos_id_dict = read_position_ID_into_dict(id_list, seq_list, seq_line)

    info_list = []
    for pep_tsv, psm_tsv in zip(pep_path_list, psm_path_list):
        file_name = pep_tsv.split('/')[-2]
        print(file_name)
        pep_list = peptide_counting(pep_tsv)
        automaton = aho_corasick.automaton_trie(pep_list)
        aho_result = aho_corasick.automaton_matching(automaton, seq_line)
        coverage_dict = identified_proteome_cov(aho_result, protein_dict)[1]
        id_pep_dict = creat_ID_pep_dict(aho_result, pos_id_dict)
        psm_dict = psm_reader(psm_tsv)[0]
        prot_spec_dict = {}
        for id in id_pep_dict:
            spec = 0
            for pep in id_pep_dict[id]:
                spec += psm_dict[pep]
            prot_spec_dict[id] = spec
        protein_id_ls = [k for k in id_pep_dict if k in ecm_prot_list]

        file_list = [[prot,
                      len(protein_dict[prot]),
                      coverage_dict[prot],
                      ecm_info_dict[prot][0],
                      prot_spec_dict[prot],
                      ecm_info_dict[prot][2],
                      file_name,
                      file_name_number_dict[file_name]] for prot in protein_id_ls]
        info_list += file_list

    info_df = pd.DataFrame(info_list,
                           columns=['protein_id', 'length', 'coverage', 'gene', 'spec_count', 'ecm_class', 'file_name',
                                    'file_number'])
    info_df.to_csv('D:/data/Naba_deep_matrisome/01102021/dash_info_new.csv')


if __name__=='__main__':

    # file_path = 'D:/data/Naba_deep_matrisome/matrisome coverage.xlsx'
    # df = pd.read_excel(file_path)
    # df = df.drop_duplicates()
    # ecm_prot_list = df['protein_id'].tolist()
    # df = df.set_index('protein_id')
    # ecm_info_dict = {each:(df.loc[each,'gene_id'], df.loc[each, 'loc'], df.loc[each,'category']) for each in ecm_prot_list}
    #
    base_path = 'D:/data/Naba_deep_matrisome/01102021/'
    folders = [f for f in os.listdir(base_path) if '.' not in f]
    psm_path_list = [base_path + each + '/psm.tsv' for each in folders]
    pep_path_list = [base_path + each + '/peptide.tsv' for each in folders]
    peptide_163_05_20 = list(set([pep for pep_file in pep_path_list[8:] for pep in peptide_counting(pep_file)]))
    # peptide_163_3A = peptide_counting(pep_path_list[7])
    print (pep_path_list[8:])

    fasta_path = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
    protein_dict = fasta_reader(fasta_path)

    #
    # dash_dataframe(pep_path_list,psm_path_list,protein_dict,ecm_prot_list,ecm_info_dict)

    from tsv_reader import venn_diagram_gen
    file_path = 'D:/data/Naba_deep_matrisome/matrisome coverage.xlsx'

    df = pd.read_excel(file_path)
    df = df.drop_duplicates()
    ecm_prot_list = df['protein_id'].tolist()
    new_protein_dict = {prot:protein_dict[prot] for prot in ecm_prot_list if prot in protein_dict}

    df = df.set_index('protein_id')

    ecm_info_dict = {each: (df.loc[each, 'gene_id'], df.loc[each, 'loc'], df.loc[each, 'category']) for each in
                     ecm_prot_list}
    id_list, seq_list = extract_UNID_and_seq(new_protein_dict)
    seq_line = creat_total_seq_line(seq_list)
    automaton = aho_corasick.automaton_trie(peptide_163_05_20)
    aho_result = aho_corasick.automaton_matching(automaton,seq_line)
    coverage_dict = identified_proteome_cov(aho_result, new_protein_dict)[1]
    print (len(coverage_dict))
    coverage_array = [v for v in coverage_dict.values()]
    print (statistics.mean(coverage_array))
    plt.hist(coverage_array,bins=15, color='black')
    plt.xlim(0,100)
    plt.xlabel('sequence coverage%')
    plt.ylabel('frequency')
    plt.title('ECM sequence coverage distribution in time lapsed 18_2B digestion')
    plt.show()
    #
    # df = pd.read_csv('D:/data/Naba_deep_matrisome/01102021/dash_info_new.csv')
    # protein_18_2A_list = df[df['file_number']==0]['protein_id'].tolist()
    # protein_1805_1820_list = list(set([prot for i in range(1,7) for prot in df[df['file_number']==i]['protein_id'].tolist()]))
    #
    # venn_dict = {'163_3_time_lapsed_digestion': protein_1805_1820_list,'163_3_20hour_digestion': protein_18_2A_list}
    # venn_diagram_gen(venn_dict)
    # print ([ecm_info_dict[prot] for prot in protein_1805_1820_list if prot not in protein_18_2A_list])

