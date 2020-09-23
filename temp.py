import aho_corasick
from protein_coverage import fasta_reader
from glob import glob
from tsv_reader import peptide_counting, psm_reader
from calculations_and_plot import whole_proteome_cov, identified_proteome_cov
import multiprocessing_naive_algorithym
import pandas as pd
import numpy as np

tryp_path = 'D:/data/deep_proteome/20200915_tryp*/peptide.tsv'
ct_path = 'D:/data/deep_proteome/20200915_ct*/peptide.tsv'
tryp_file_list = glob(tryp_path)
ct_file_list = glob(ct_path)

tryp_psm_path = 'D:/data/deep_proteome/20200915_tryp*/psm.tsv'
ct_psm_path = 'D:/data/deep_proteome/20200915_ct*/psm.tsv'
tryp_psm_file_list = glob(tryp_psm_path)
ct_psm_file_list = glob(ct_psm_path)

fasta_path = 'D:/data/proteome_fasta/uniprot-proteome_UP000005640.fasta'
protein_dict = fasta_reader(fasta_path)
print (len(protein_dict))
ID_list, seq_list = multiprocessing_naive_algorithym.extract_UNID_and_seq(protein_dict)
seq_line = multiprocessing_naive_algorithym.creat_total_seq_line(seq_list)
coverage_dict = {}
df = pd.DataFrame()
pd.set_option('display.max_rows', None)

# whole proteome coverage calculation, singe file and combined
# for each_tryp_file in tryp_file_list:
#     for each_ct_file in ct_file_list:
#         print (each_tryp_file, each_ct_file)
#
#         enzyme_tm_time_tryp, enzyme_tm_time_ct = '_'.join(each_tryp_file.split('\\')[-2].split('_')[1:]), \
#                                                  '_'.join(each_ct_file.split('\\')[-2].split('_')[1:])
#         tryp_peplist, ct_peplist = peptide_counting(each_tryp_file), peptide_counting(each_ct_file)
#         #combined_peplist = list(set(tryp_peplist+ct_peplist))
#         for each_peplist, cond in zip([tryp_peplist, ct_peplist],
#                                       [enzyme_tm_time_tryp, enzyme_tm_time_ct]):
#         #for each_peplist, cond in zip([tryp_peplist,ct_peplist,combined_peplist],[enzyme_tm_time_tryp, enzyme_tm_time_ct, enzyme_tm_time_tryp+'_'+enzyme_tm_time_ct]):
#             automaton = aho_corasick.automaton_trie(each_peplist)
#             aho_result = aho_corasick.automaton_matching(automaton,seq_line)
#             proteome_coverage = whole_proteome_cov(aho_result,protein_dict)
#             print (cond, proteome_coverage)
#             df.loc[cond,'cov'] = proteome_coverage

# whole proteome coverage caluculation, single file
# total_file_list = tryp_file_list+ct_file_list
# for each_file in total_file_list:
#     print (each_file)
#
#     ez_tm = '_'.join(each_file.split('\\')[-2].split('_')[1:3])
#     time = each_file.split('\\')[-2].split('_')[-1]
#     print(ez_tm,time)
#     pep_list = peptide_counting(each_file)
#     automaton = aho_corasick.automaton_trie(pep_list)
#     aho_result = aho_corasick.automaton_matching(automaton,seq_line)
#     proteome_coverage = whole_proteome_cov(aho_result,protein_dict)
#     print (ez_tm,time, proteome_coverage)
#     df.loc[ez_tm,time] = proteome_coverage

# identified protein coverage, single file
# total_file_list = tryp_file_list+ct_file_list
# for each_file in total_file_list:
#     print (each_file)
#     cond = '_'.join(each_file.split('\\')[-2].split('_')[1:])
#     # ez_tm = '_'.join(each_file.split('\\')[-2].split('_')[1:3])
#     # time = each_file.split('\\')[-2].split('_')[-1]
#     print(cond)
#     pep_list = peptide_counting(each_file)
#     automaton = aho_corasick.automaton_trie(pep_list)
#     aho_result = aho_corasick.automaton_matching(automaton,seq_line)
#     identified_proteome_cov_dict = identified_proteome_cov(aho_result,protein_dict)[1]
#     for each_prot in identified_proteome_cov_dict:
#         df.loc[each_prot,cond] = identified_proteome_cov_dict[each_prot]
# df = df.fillna(0)

#psm counting
# total_file_list = tryp_psm_file_list+ct_psm_file_list
# for each_file in total_file_list:
#     print (each_file)
#     psm_count=0
#     ez_tm = '_'.join(each_file.split('\\')[-2].split('_')[1:3])
#     time = each_file.split('\\')[-2].split('_')[-1]
#     print(ez_tm,time)
#     psm_dict = psm_reader(each_file,fragpipe_ver=13.0)[0]
#     for each in psm_dict:
#         psm_count+=psm_dict[each]
#     df.loc[ez_tm,time] = psm_count

# compare the same identified proteome between files
# for each_tryp_file in tryp_file_list:
#     print (each_tryp_file)
#     ez_tm_time_tryp = '_'.join(each_tryp_file.split('\\')[-2].split('_')[1:])
#     tryp_peplist = peptide_counting(each_tryp_file)
#     automaton_tryp = aho_corasick.automaton_trie(tryp_peplist)
#     aho_result_tryp = aho_corasick.automaton_matching(automaton_tryp,seq_line)
#     average_proteome_coverage_tryp, identified_proteome_cov_dict_tryp = identified_proteome_cov(aho_result_tryp,protein_dict)[:2]
#     print (average_proteome_coverage_tryp)
#     df.loc[ez_tm_time_tryp,'cov'] = average_proteome_coverage_tryp
#     for each_ct_file in ct_file_list:
#         print (each_ct_file)
#         ez_tm_time_ct = '_'.join(each_ct_file.split('\\')[-2].split('_')[1:])
#         ct_peplist = peptide_counting(each_ct_file)
#         combined_peplist= list(set(tryp_peplist+ct_peplist))
#         automaton_comb = aho_corasick.automaton_trie(combined_peplist)
#         aho_result_comb = aho_corasick.automaton_matching(automaton_comb,seq_line)
#         identified_proteome_cov_dict_comb = identified_proteome_cov(aho_result_comb,protein_dict)[1]
#         new_id_tryp_cov_dict = {each:identified_proteome_cov_dict_comb[each] for each in identified_proteome_cov_dict_tryp}
#         new_tryp_average_cov = np.mean([v for v in new_id_tryp_cov_dict.values()])
#         print (new_tryp_average_cov)
#         df.loc[ez_tm_time_tryp+'_'+ez_tm_time_ct,'cov'] = new_tryp_average_cov



# df.to_excel('D:/data/deep_proteome/9_20_coverage.xlsx')
# df = pd.read_excel('D:/data/deep_proteome/9_20_coverage.xlsx')
# print (df)
# cols = df.columns.tolist()
# new_cols = ['Unnamed: 0', '0min', '10min','20min','30min','60min', '120min', '240min', '1440min', '3180min', '4320min']
# df = df.reindex(columns=new_cols)
# df.to_excel('D:/data/deep_proteome/9_20_identified_cov_total.xlsx')

# use the protein order from total table and calculate coverage after data combined