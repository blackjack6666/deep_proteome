import pandas as pd
from multiprocessing_naive_algorithym import extract_UNID_and_seq,creat_total_seq_line,creat_ID_pep_dict,\
    read_position_ID_into_dict
from tsv_reader import psm_reader,peptide_counting
from calculations_and_plot import identified_proteome_cov
import aho_corasick
import os
from protein_coverage import fasta_reader


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

    file_path = 'D:/data/Naba_deep_matrisome/matrisome coverage.xlsx'
    df = pd.read_excel(file_path)
    df = df.drop_duplicates()
    ecm_prot_list = df['protein_id'].tolist()
    df = df.set_index('protein_id')
    ecm_info_dict = {each:(df.loc[each,'gene_id'], df.loc[each, 'loc'], df.loc[each,'category']) for each in ecm_prot_list}

    base_path = 'D:/data/Naba_deep_matrisome/01102021/'
    folders = [f for f in os.listdir(base_path) if '.' not in f]
    psm_path_list = [base_path + each + '/psm.tsv' for each in folders]
    pep_path_list = [base_path + each + '/peptide.tsv' for each in folders]
    fasta_path = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
    protein_dict = fasta_reader(fasta_path)

    dash_dataframe(pep_path_list,psm_path_list,protein_dict,ecm_prot_list,ecm_info_dict)