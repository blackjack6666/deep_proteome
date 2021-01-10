import pandas as pd
from protein_coverage import fasta_reader, read_description
import aho_corasick
import numpy as np
from glob import glob
from tsv_reader import peptide_counting
from calculations_and_plot import whole_proteome_cov, identified_proteome_cov
import multiprocessing_naive_algorithym

#custom_db build

# file_path = 'C:/uic/lab/data/naba/matrisome coverage.xlsx'
#
# df = pd.read_excel(file_path)
# df = df.drop_duplicates()
# print(df.shape)
#
# ecm_protein_id_set = set(df['protein_id'].to_list())
# print (len(ecm_protein_id_set))
#
# fasta_path = 'D:/data/proteome_fasta/uniprot-proteome_UP000000589_mouse.fasta'
# protein_dict = fasta_reader(fasta_path)
# protein_descript_dict = read_description(fasta_path)
#
# count=0
# with open('D:/data/Naba_deep_matrisome/mouse_ecm_costom_proteome_db.fasta', 'w', newline='\n') as file_w:
#     for each in ecm_protein_id_set:
#         if each in protein_dict:
#            count+=1
#            seq = protein_dict[each]
#            prefix,description = protein_descript_dict[each]
#            file_w.write('>'+str(prefix)+'|'+str(each)+'|'+str(description)+'\n')
#            bins = range(0,len(seq)+60,60)
#            for i in range(len(bins)-1):
#                file_w.write(seq[bins[i]:bins[i+1]]+'\n')
#
# print (count,'fasta done')

# coverage calculation

# fasta_path = 'D:/data/Naba_deep_matrisome/mouse_ecm_costom_proteome_db.fasta'
# protein_dict = fasta_reader(fasta_path)
# ID_list, seq_list = multiprocessing_naive_algorithym.extract_UNID_and_seq(protein_dict)
# seq_line = multiprocessing_naive_algorithym.creat_total_seq_line(seq_list)
#
# trypsin_peptide_tsv = 'D:/data/Naba_deep_matrisome/09_10_2020_T_search_result/peptide.tsv'
# chymot_peptide_tsv = 'D:/data/Naba_deep_matrisome/09_10_2020_chymoT_search_result/peptide.tsv'
#
# tryp_peptide_list,chymo_peptide_list = peptide_counting(trypsin_peptide_tsv),peptide_counting(chymot_peptide_tsv)
# combined = tryp_peptide_list+chymo_peptide_list
#
# all_protein_cov_dict_protease = {} # different protease as key
#
# for each_pep_list,protease in zip([tryp_peptide_list,chymo_peptide_list,combined],['trypsin','chymotrypsin','T_CT_combined']):
#     automaton = aho_corasick.automaton_trie(each_pep_list)
#     aho_result = aho_corasick.automaton_matching(automaton,seq_line)
#     proteome_coverage = whole_proteome_cov(aho_result, protein_dict)
#     identifed_coverage, idenfied_cov_dict,all_protein_cov_dict = identified_proteome_cov(aho_result,protein_dict)
#     all_protein_cov_dict_protease[protease] = all_protein_cov_dict
#     print (protease, proteome_coverage, identifed_coverage )
#
# trypsin_automaton = aho_corasick.automaton_trie(tryp_peptide_list)
# trypsin_aho_result = aho_corasick.automaton_matching(trypsin_automaton,seq_line)
# trypsin_idenfied_dict = identified_proteome_cov(trypsin_aho_result,protein_dict)[1]
#
# combined_automaton = aho_corasick.automaton_trie(combined)
# combined_aho_result = aho_corasick.automaton_matching(combined_automaton,seq_line)
# combined_identified_dict = identified_proteome_cov(combined_aho_result,protein_dict)[1]
# new_identified_dict = {each:combined_identified_dict[each] for each in trypsin_idenfied_dict}
#
#
# cov_list = [[id,all_protein_cov_dict_protease['trypsin'][id],
#              all_protein_cov_dict_protease['chymotrypsin'][id],
#              all_protein_cov_dict_protease['T_CT_combined'][id]]for id in protein_dict]
#
# df = pd.DataFrame(cov_list, columns=['protein id', 'trypsin cov', 'chymotrypsin cov', 'T_CT_cov'])
# df.to_excel('D:/data/Naba_deep_matrisome/coverage_9_11_2020.xlsx')

# collagen length distribution

# import matplotlib.pyplot as plt
#
# file_path = 'C:/uic/lab/data/naba/matrisome coverage.xlsx'
# human_fasta = 'C:/uic/lab/data/proteome_fasta/uniprot-proteome_UP000005640.fasta'
# mouse_fasta = 'C:/uic/lab/data/proteome_fasta/uniprot-proteome_UP000000589_mouse.fasta'
#
# human_dict,mouse_dict = fasta_reader(human_fasta), fasta_reader(mouse_fasta)
#
# df = pd.read_excel(file_path)
# df = df.drop_duplicates()
# print (df.shape)
# df_human = df[df['protein_id'].isin(human_dict)]
# print (df_human.shape)
# df_mouse = df[df['protein_id'].isin(mouse_dict)]
# print (df_mouse.shape)
#
# human_collagen_list = df_human[df_human['category']=='Collagens']['protein_id'].tolist()
# mouse_collagen_list = df_mouse[df_mouse['category']=='Collagens']['protein_id'].tolist()
#
# human_collagen_length_list = [len(human_dict[each]) for each in human_collagen_list]
# mouse_collagen_length_list = [len(mouse_dict[each]) for each in mouse_collagen_list]
# print (np.mean(human_collagen_length_list),np.mean(mouse_collagen_length_list))
#
#
#
# human_core_list = [len(human_dict[each]) for each in df_human[df_human['loc']=='Core matrisome']['protein_id'].tolist()]
# mouse_core_list = [len(mouse_dict[each]) for each in df_mouse[df_mouse['loc']=='Core matrisome']['protein_id'].tolist()]
# print (np.mean(human_core_list),np.mean(mouse_core_list))
#
# fig,(ax0,ax1) = plt.subplots(1,2,figsize=(15,10))
# n_bins = 50
#
# ax0.hist(human_core_list,bins=n_bins)
# ax0.set_title('human core matrisome length dist')
# ax1.hist(mouse_core_list,bins=n_bins)
# ax1.set_title('mouse core matrisome length dist')
#
# plt.show()

from glob import glob
from collections import defaultdict
import os
import aho_corasick
from multiprocessing_naive_algorithym import extract_UNID_and_seq, creat_total_seq_line, creat_ID_pep_dict,read_position_ID_into_dict
from tsv_reader import peptide_counting, map_psm_file, protein_info_from_combined, protein_info_from_fasta, psm_reader
from protein_coverage import fasta_reader
from pandas import ExcelWriter

fasta_path = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
protein_dict = fasta_reader(fasta_path)
id_list,seq_list = extract_UNID_and_seq(protein_dict)
seq_line = creat_total_seq_line(seq_list)
pos_id_dict = read_position_ID_into_dict(id_list,seq_list,seq_line)
protein_info_dict = protein_info_from_fasta(fasta_path)
print ('done')


# write protein IDs from different exp into excel
base_path = 'D:/data/Naba_deep_matrisome/01102021/'
folders = [f for f in os.listdir(base_path) if '.' not in f]
psm_path_list = [base_path+each+'/psm.tsv' for each in folders]
pep_path_list = [base_path+each+'/peptide.tsv' for each in folders]
protein_path_list = [base_path+each+'/protein.tsv' for each in folders]



with ExcelWriter('D:/data/Naba_deep_matrisome/011021_protein_ids.xlsx') as writer:
    for pep_path in pep_path_list:
        file_name = pep_path.split('/')[-2]
        print (file_name)
        pep_list = peptide_counting(pep_path)
        automaton = aho_corasick.automaton_trie(pep_list)
        aho_result = aho_corasick.automaton_matching(automaton, seq_line)
        id_pep_dict = creat_ID_pep_dict(aho_result, pos_id_dict)
        protein_id_ls = [k for k in id_pep_dict]
        info_list = [[prot,
                      protein_info_dict[prot][0],
                      protein_info_dict[prot][1]]
                      for prot in protein_id_ls]
        df = pd.DataFrame(info_list, columns=['protein ID', 'gene name', 'description'])
        df.to_excel(writer, '%s' % file_name)

# with ExcelWriter('D:/data/Naba_deep_matrisome/1218_protein_ids.xlsx') as writer:
#     for each_psm_path in psm_path_list:
#         file_pep_dict = map_psm_file(each_psm_path)
#         # print ([(k,len(file_pep_dict[k])) for k in file_pep_dict])
#         for each_file in file_pep_dict:
#             pep_list = list(file_pep_dict[each_file])
#             automaton = aho_corasick.automaton_trie(pep_list)
#             aho_result = aho_corasick.automaton_matching(automaton,seq_line)
#             id_pep_dict = creat_ID_pep_dict(aho_result,pos_id_dict)
#             protein_id_ls = [k for k in id_pep_dict]
#
#             info_list = [[prot,
#                           protein_info_dict[prot][0],
#                           protein_info_dict[prot][1]]
#                          for prot in protein_id_ls]
#             df = pd.DataFrame(info_list, columns=['protein ID', 'gene name', 'description'])
#             df.to_excel(writer,'%s' % each_file)



