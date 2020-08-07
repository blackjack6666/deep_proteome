import aho_corasick
from protein_coverage import fasta_reader
from glob import glob
from tsv_reader import peptide_counting
from calculations_and_plot import whole_proteome_cov, identified_proteome_cov
import multiprocessing_naive_algorithym
import pandas as pd
import numpy as np

tryp_path = 'D:/data/deep_proteome/20200806_tryp*/peptide.tsv'
ct_path = 'D:/data/deep_proteome/20200806_CT*/peptide.tsv'
tryp_file_list = glob(tryp_path)
ct_file_list = glob(ct_path)

fasta_path = 'D:/data/proteome_fasta/uniprot-proteome_UP000005640.fasta'
protein_dict = fasta_reader(fasta_path)
ID_list, seq_list = multiprocessing_naive_algorithym.extract_UNID_and_seq(protein_dict)
seq_line = multiprocessing_naive_algorithym.creat_total_seq_line(seq_list)
coverage_dict = {}
df = pd.DataFrame()
pd.set_option('display.max_rows', None)


for each_tryp_file in tryp_file_list:
    for each_ct_file in ct_file_list:
        print (each_tryp_file, each_ct_file)

        enzyme_tm_time_tryp, enzyme_tm_time_ct = '_'.join(each_tryp_file.split('\\')[-2].split('_')[1:]), \
                                                 '_'.join(each_ct_file.split('\\')[-2].split('_')[1:])
        tryp_peplist, ct_peplist = peptide_counting(each_tryp_file), peptide_counting(each_ct_file)
        combined_peplist = list(set(tryp_peplist+ct_peplist))
        for each_peplist, cond in zip([tryp_peplist,ct_peplist,combined_peplist],[enzyme_tm_time_tryp, enzyme_tm_time_ct, enzyme_tm_time_tryp+'_'+enzyme_tm_time_ct]):
            automaton = aho_corasick.automaton_trie(each_peplist)
            aho_result = aho_corasick.automaton_matching(automaton,seq_line)
            proteome_coverage = whole_proteome_cov(aho_result,protein_dict)
            print (cond, proteome_coverage)
            df.loc[cond,'cov'] = proteome_coverage

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



df.to_excel('8_7_whole_proteome_cov.xlsx')
