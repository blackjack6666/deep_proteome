from tsv_reader import peptide_counting, id_pep_from_peptsv
from protein_coverage import fasta_reader2, fasta_reader
import re
import pandas as pd
from collections import defaultdict
import numpy as np
from data_preprocess_3 import polymer_gen
import matplotlib.pyplot as plt
import multiprocessing

fasta_path = 'D:/data/proteome_fasta/uniprot-proteome_UP000005640_reverse_beta.fasta'
protein_dict = fasta_reader(fasta_path)
# pep_tsv = 'D:/data/deep_proteome/non_specfic_search/ct_4h/peptide.tsv'
psm_tsv = 'D:/data/deep_proteome/non_specfic_search/tryps_4h/psm.tsv'
#

# id_peptide_dict = id_pep_from_peptsv(pep_tsv)
df_psm = pd.read_csv(psm_tsv,sep='\t')
psm_score_dict = defaultdict(list)  # psm and hyperscore
for each_row in df_psm.itertuples():
    psm_score_dict[each_row[3]].append(each_row[16])
pep_average_score_dict = {each:np.mean(psm_score_dict[each]) for each in psm_score_dict}

# check if pep fits customized regex cleavage rule
# rule = r'([KR](?=[^P]))'
#
# matched_peptide = 0
# no_matched_peptide = 0
# matched_peptide_score_dict, no_matched_peptide_score_dict = {}, {}
# for prot in id_peptide_dict:
#     prot_seq = protein_dict[prot]
#     for pep in id_peptide_dict[prot]:
#         pep_loc = prot_seq.find(pep) - 1
#         if pep_loc!=-2:
#             pep_end_loc = pep_loc+len(pep)
#             start_31mer, end_31_mer = polymer_gen(prot_seq,pep_loc), polymer_gen(prot_seq,pep_end_loc)
#             start_match, end_match = [m.end() for m in re.finditer(rule,start_31mer)], [m.end() for m in re.finditer(rule,end_31_mer)]
#             if 16 in start_match and 16 in end_match:  # start and end position of a peptide matches with cleavage rule
#                 matched_peptide+=1
#                 matched_peptide_score_dict[pep] = pep_average_score_dict[pep]
#             elif (16 not in start_match) and (16 not in end_match):
#             #
#                 print (prot,pep, start_31mer, start_match,end_31_mer, end_match)
#                 no_matched_peptide+=1
#                 no_matched_peptide_score_dict[pep] = pep_average_score_dict[pep]
#
#
# matched_hyper_values = [v for v in matched_peptide_score_dict.values()]
# no_matched_hyper_values = [v for v in no_matched_peptide_score_dict.values()]
#
# print (np.mean(matched_hyper_values), np.mean(no_matched_hyper_values))

# plt.hist(matched_hyper_values, bins=50, alpha=0.5, label='regex-matched psm, average score=20.6',log=True)
# plt.hist(no_matched_hyper_values,bins=50, alpha=0.5, color='pink', label='regex-filtered psm, average score=16.2',log=True)
# plt.xlabel('Hyperscore')
# plt.ylabel('Frequency')
# plt.legend()
# plt.show()


# from scipy import stats
# print (stats.ttest_ind(matched_hyper_values, no_matched_hyper_values))

from tsv_reader import pep_xml_info
import pickle as ppp

pepxml = 'D:/data/Mankin/api_05_search_data/XS_Shura_Ribo_Api_0_5.pepXML'
head, tail, info_list = pep_xml_info(pepxml)
print (len(info_list))
custom_tryptic_peptide_dict_ecoli = ppp.load(open('D:/data/deep_proteome/non_specfic_search/custom_trypsin_human4h_1.p','rb'))
custom_pep_set = {pep for v in custom_tryptic_peptide_dict_ecoli.values() for pep in v}
print (len(custom_pep_set))
count = 0

fit_regex_score_list = [score for each in psm_score_dict for score in psm_score_dict[each] if each in custom_pep_set]
no_fit_regex_score_list = [score for each in psm_score_dict for score in psm_score_dict[each] if each not in custom_pep_set]
print (len(fit_regex_score_list),len(no_fit_regex_score_list))
print (np.mean(fit_regex_score_list), np.mean(no_fit_regex_score_list))
plt.hist(fit_regex_score_list, bins=100, alpha=0.5, label='regex-matched psm, average score=20.6')
plt.hist(no_fit_regex_score_list,bins=100, alpha=0.5, color='pink', label='regex-filtered psm, average score=16.2')
plt.xlabel('Hyperscore')
plt.ylabel('Frequency')
plt.legend()
# plt.boxplot([fit_regex_score_list,no_fit_regex_score_list],labels=['fit_regex', 'no_fit'])
plt.show()
"""
# write pepXML files
pepxml_fitrule = 'D:/data/Mankin/api_05_search_data_regexfit/XS_Shura_Ribo_Api_0_5_fit.pepXML'
pepxml_nofit = 'D:/data/Mankin/api_05_search_data_noregexfit/XS_Shura_Ribo_Api_0_5_no_fit.pepXML'

# create two pepxml files, write head first
f_fitrule = open(pepxml_fitrule,'w',newline='\n')
f_no_fit = open(pepxml_nofit,'w',newline='\n')
f_fitrule.write(head)
f_no_fit.write(head)
#if peptide hit from original pepxml file is in custom peptide set, write it to fit fitrule file,otherwise write it to no fit rule file
for each in info_list:
    peptide = each.split('peptide="')[1].split('"')[0]
    if peptide in custom_pep_set:
        f_fitrule.write('<spectrum_query start_scan='+each+'</spectrum_query>\n')
    else:
        f_no_fit.write('<spectrum_query start_scan='+each+'</spectrum_query>\n')
# fit_rule_count = 0
# for each in info_list:
#     peptide = each.split('peptide="')[1].split('"')[0]
#     protein_id = each.split('|')[1]
#     prot_seq = protein_dict[protein_id]
#     pep_loc = prot_seq.find(peptide) - 1
#     if pep_loc != -2:
#         pep_end_loc = pep_loc + len(peptide)
#         start_31mer, end_31_mer = polymer_gen(prot_seq, pep_loc), polymer_gen(prot_seq, pep_end_loc)
#         start_match, end_match = [m.end() for m in re.finditer(rule, start_31mer)], [m.end() for m in
#                                                                                      re.finditer(rule, end_31_mer)]
#         if pep_loc == -1  and  16 in end_match:
#             fit_rule_count+=1
#         elif pep_end_loc == len(prot_seq)-1 and 16 in start_match:
#             fit_rule_count+=1
#         elif 16 in start_match and 16 in end_match:
#             fit_rule_count+=1
#
#         else:
#             print (peptide, protein_id)
# print (fit_rule_count)


f_fitrule.write(tail.lstrip('\n'))
f_no_fit.write(tail.lstrip('\n'))
f_fitrule.close()
f_no_fit.close()

"""