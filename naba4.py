from tsv_reader import psm_reader, venn_diagram_gen
import pickle as p
from glob import glob
from collections import defaultdict
from dash_dataframe_naba import dash_dataframe,combined_coverage
import os
import pandas as pd
from protein_coverage import fasta_reader
import matplotlib.pyplot as plt
from math import log10
from matplotlib.lines import Line2D
import numpy as np
from pandas import ExcelWriter

folders_182 = ['18_2B05','18_2B1','18_2B2','18_2B4','18_2B18','18_2B20']
folders_163 = ['163_3B05','163_3B1','163_3B2','163_3B4','163_3B18','163_3B20']
path = 'D:/data/Naba_deep_matrisome/02152021_1/'

# psm_dict_182 = defaultdict(int)
# for each in folders_163:
#     psm_tsv = path+each+'/psm.tsv'
#     psm_dict = psm_reader(psm_tsv)[0]
#     for pep in psm_dict:
#         psm_dict_182[pep]+=psm_dict[pep]



# p.dump(psm_reader('D:/data/Naba_deep_matrisome/02152021_1/163_3A/psm.tsv')[0],open('D:/data/Naba_deep_matrisome/02152021_1/163A_psm_dict.p','wb'))

#getting coverage and other info matrix
fasta_path = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
protein_dict = fasta_reader(fasta_path)

file_path = 'D:/data/Naba_deep_matrisome/matrisome coverage.xlsx'
df = pd.read_excel(file_path)
df = df.drop_duplicates()
ecm_prot_list = df['protein_id'].tolist()
df = df.set_index('protein_id')
ecm_info_dict = {each:(df.loc[each,'gene_id'], df.loc[each, 'loc'], df.loc[each,'category']) for each in ecm_prot_list}
#
base_path = 'D:/data/Naba_deep_matrisome/05142021_secondsearch/SNEDC/'
folders = [f for f in os.listdir(base_path) if '.' not in f]
psm_path_list = [base_path + each + '/psm.tsv' for each in folders]
pep_path_list = [base_path + each + '/peptide.tsv' for each in folders]

info_dict = dash_dataframe(pep_path_list,psm_path_list,protein_dict,ecm_prot_list,ecm_info_dict,'D:/data/Naba_deep_matrisome/03122021/SNEDC/SNEDC_ecm_info_matrix.csv')
with ExcelWriter('D:/data/Naba_deep_matrisome/05142021_secondsearch/SNEDC//protein_spec_SNEDC.xlsx') as writer:
    for file in info_dict:
        info_list = [[prot, info_dict[file][prot], ecm_info_dict[prot][0], ecm_info_dict[prot][1],ecm_info_dict[prot][2]] for prot in info_dict[file]]
        df = pd.DataFrame(info_list, columns=['protein id', 'psm with count', 'gene name', 'loc', 'category'])
        df.to_excel(writer, '%s' % file)

# coverage line plot each category
# category_protein_id_dict=p.load(open('D:/data/Naba_deep_matrisome/02152021_1/ecm_category_protein_dict.p','rb'))
# accumulated_cov_163_3_df = pd.read_excel('D:/data/Naba_deep_matrisome/02152021_1/163_accumulated_norm.xlsx',index_col=0)
# print (accumulated_cov_163_3_df)
# protein_idx = accumulated_cov_163_3_df.columns.tolist()
# fig,axs = plt.subplots(6,1, figsize=(10,20))
# for each_cat,ax in zip(category_protein_id_dict,axs):
#     protein_list = [prot for prot in category_protein_id_dict[each_cat] if prot in protein_idx]
#     sub_df = accumulated_cov_163_3_df[protein_list]
#
#     # x = range(len(sub_df.index))
#     x = [log10(0.5),log10(1),log10(2),log10(4),log10(18),log10(20)]
#     for column in sub_df.columns:
#         y = sub_df[column].values
#         y_log10 = np.log10(y)
#         print (y)
#         if y[-1]-y[0] == 0: # if line is flat, make it transparent
#             ax.plot(x,y_log10, color='lightcoral', linestyle='-',alpha=0.1)
#         else:
#             ax.plot(x,y_log10, color='lightcoral', linestyle='-')
#
#     ax.set_xticks([log10(0.5),log10(1),log10(2),log10(4),log10(18),log10(20)])
#     ax.set_xticklabels(['0.5h','1h','2h','4h','18h', '20h'], fontsize=10)
#     ax.set_title(each_cat, fontsize=12)
#     ax.set_xlabel('log time point', fontsize=12)
#     ax.set_ylabel('log %coverage',fontsize=12)
#     ax.set(yticklabels=[])
#
# accumulated_cov_18_2_df = pd.read_csv('D:/data/Naba_deep_matrisome/02152021_1/18_2_accumulated_norm.csv',
#                                              index_col=0)
# protein_18_2_index = accumulated_cov_18_2_df.columns.tolist()
# for each_cat,ax in zip(category_protein_id_dict,axs):
#     protein_list = [prot for prot in category_protein_id_dict[each_cat] if prot in protein_18_2_index]
#     sub_df = accumulated_cov_18_2_df[protein_list]
#
#     # x = range(len(sub_df.index))
#     x = [log10(0.5),log10(1),log10(2),log10(4),log10(18),log10(20)]
#     for column in sub_df.columns:
#         y = sub_df[column].values
#         y_log10 = np.log10(y)
#         print (y)
#         if y[-1]-y[0] == 0: # if line is flat, make it transparent
#             ax.plot(x,y_log10, color='lightseagreen', linestyle='-',alpha=0.1)
#         else:
#             ax.plot(x,y_log10, color='lightseagreen', linestyle='-')
#
# colors, labels, lines= ['lightcoral','lightseagreen'], ['SNED1 OE', 'SNED1 KO'], [Line2D([0], [0], color=c, linewidth=3, linestyle='-') for c in ['lightcoral','lightseagreen']]
#
# fig.legend(lines, labels, framealpha=0.5, loc='upper right', fontsize=15)
#     # sub_df.to_csv('D:/data/Naba_deep_matrisome/02152021_1/'+'163_3_'+str(each_cat)+'_cov.csv')
#     # print (each_cat,sub_df)
#
# frame1 = plt.gca()
#
# frame1.axes.yaxis.set_ticklabels([])
# plt.show()

# df = pd.read_csv('D:/data/Naba_deep_matrisome/02152021_1/dash_info.csv')
# df_coverage = pd.read_csv('D:/data/Naba_deep_matrisome/02152021_1/coverage_boxplot_dataset.csv',index_col=0)
# print (df_coverage.index)
#
# sub_df_182 = df[(df['file_name']=='18_2B05')|(df['file_name']=='18_2B1')|(df['file_name']=='18_2B2')|
# (df['file_name']=='18_2B4')|(df['file_name']=='18_2B18')|(df['file_name']=='18_2B20')]
#
# spec_count_dict_182 = defaultdict(int)
# for prot, spec_count in zip(sub_df_182['protein_id'],sub_df_182['spec_count']):
#     spec_count_dict_182[prot] += spec_count
#
# sub_df_163 = df[(df['file_name']=='163_3B05')|(df['file_name']=='163_3B1')|(df['file_name']=='163_3B2')|
# (df['file_name']=='163_3B4')|(df['file_name']=='163_3B18')|(df['file_name']=='163_3B20')]
#
# spec_count_dict_163 = defaultdict(int)
# for prot, spec_count in zip(sub_df_163['protein_id'], sub_df_163['spec_count']):
#     spec_count_dict_163[prot]+=spec_count
#
# df_coverage['sned_ko_time_series_total_spec_count'] = [spec_count_dict_182[prot] for prot in df_coverage.index]
# df_coverage['sned_oe_time_series_total_spec_count'] = [spec_count_dict_163[prot] for prot in df_coverage.index]
#
# df_coverage.to_excel('D:/data/Naba_deep_matrisome/02152021_1/coverage_boxplot_dataset_1.xlsx')

# combined coverage output
# base_path = 'D:/data/Naba_deep_matrisome/03122021/SNEDC/'
# folders = [f for f in os.listdir(base_path) if '.' not in f]
#
# pep_path_list = [base_path + each + '/peptide.tsv' for each in folders]
# ko_05_4_peplist = [pep_path_list[0]]+[pep_path_list[2]]+[pep_path_list[3]]
# print (ko_05_4_peplist)
# ko_18_peplist = [pep_path_list[1]]
#
"""
file_path = 'D:/data/Naba_deep_matrisome/matrisome coverage.xlsx'
df = pd.read_excel(file_path)
df = df.drop_duplicates()
print (len(df['gene_id'].unique()))
ecm_prot_list = df['protein_id'].tolist()
df = df.set_index('protein_id')
ecm_info_dict = {each:(df.loc[each,'gene_id'], df.loc[each, 'loc'], df.loc[each,'category']) for each in ecm_prot_list}
fasta_path = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
protein_dict = fasta_reader(fasta_path)
#
#
# ko_05_4_cov_dict = combined_coverage(ko_05_4_peplist,protein_dict,ecm_prot_list)
# ko_18_cov_dict = combined_coverage(ko_18_peplist,protein_dict,ecm_prot_list)
#
# info_list = [[each,ko_05_4_cov_dict[each],ecm_info_dict[each][0],ecm_info_dict[each][-1]] for each in ko_05_4_cov_dict]
#
# df_info = pd.DataFrame(info_list, columns=['protein_id','coverage','gene','category'])
# df_info.to_excel('D:/data/Naba_deep_matrisome/03122021/SNEDC/SNEDC_05_to_4_coverage.xlsx')
# print (len([prot for prot in ko_05_4_cov_dict if prot in ko_18_cov_dict]))
# print (len(ko_05_4_cov_dict),len(ko_18_cov_dict))
# print (np.mean([v for v in ko_05_4_cov_dict.values()]), np.mean([v for v in ko_18_cov_dict.values()]))

# peptide length
# id_peptide_dict_18_2 = p.load(open('18_2_id_pepdict_0215.p','rb'))
# id_peptide_dict = id_peptide_dict_18_2['18_2B4']
# ko_05_peptide_len = [len(pep) for prot in id_peptide_dict for pep in id_peptide_dict[prot]]
# print (len(ko_05_peptide_len),np.mean(ko_05_peptide_len))
# plt.hist(ko_05_peptide_len)
# plt.show()

# from coverage_nsaf import protein_dataframe
# from tsv_reader import peptide_counting
# file_nsafdf_dict = {}
# base_path = 'D:/data/Naba_deep_matrisome/02152021_1/'
# fasta_path = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
# folders = [f for f in os.listdir(base_path) if '.' not in f]
# pep_path_list = [base_path + each + '/peptide.tsv' for each in folders]
# for each in pep_path_list:
#     file_name = each.split('/')[-2]
#     print(file_name)
#     file_nsafdf_dict[file_name] = protein_dataframe(peptide_counting(each),fasta_path)
# p.dump(file_nsafdf_dict,open('D:/data/Naba_deep_matrisome/02152021_1/file_nsafdf_dict.p','wb'))

df_info = pd.read_csv('D:/data/Naba_deep_matrisome/02152021_1/dash_info.csv')
ecm_prot_list = set(df_info['protein_id'].tolist())
file_nsaf_dict = p.load(open('D:/data/Naba_deep_matrisome/02152021_1/file_nsafdf_dict.p','rb'))
cluster_df=pd.DataFrame()
# for prot in ecm_prot_list:
#     for file,time in zip(['18_2B05','18_2B1','18_2B2','18_2B4','18_2B18','18_2B20'], ['0.5h','1h','2h','4h','18h','20h']):
#         nsaf_df = file_nsaf_dict[file]
#         protein_list = nsaf_df['protein_id'].tolist()
#         if prot in protein_list:
#             nsaf = nsaf_df.loc[nsaf_df['protein_id']==prot,'NSAF'].values[0]
#             cluster_df.loc[prot,time] = nsaf
#         else:
#             cluster_df.loc[prot,time] = 0
# cluster_df['category'] = [ecm_info_dict[prot][-1] for prot in cluster_df.index]
# cluster_df.to_csv('D:/data/Naba_deep_matrisome/02152021_1/cluster_df_snedKO.csv')

accumulated_cov_182_df = pd.read_csv('D:/data/Naba_deep_matrisome/02152021_1/18_2_accumulated_cov.csv',index_col=0)
accumulated_cov_182_df_transpose = accumulated_cov_182_df.transpose()
accumulated_cov_182_df_transpose['category'] = [ecm_info_dict[ind][-1] for ind in accumulated_cov_182_df_transpose.index]
accumulated_cov_182_df_transpose.to_csv('D:/data/Naba_deep_matrisome/02152021_1/18_2_accumulated_cov_transpose.csv')
"""