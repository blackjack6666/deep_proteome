import pandas as pd
import pickle as p
from tsv_reader import venn_diagram_gen
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import re

matrisome_cov_csv = 'D:/data/Naba_deep_matrisome/matrisome coverage.xlsx'
df_mat = pd.read_excel(matrisome_cov_csv)
df_mat = df_mat.drop_duplicates()
ecm_proteins = df_mat['protein_id'].tolist()
# ecm_category = df_mat['category'].tolist()
# category_protein_id_dict = defaultdict(list)
# for prot,category in zip(ecm_proteins,ecm_category):
#     category_protein_id_dict[category].append(prot)


# p.dump(category_protein_id_dict,open('ecm_category_protein_dict.p','wb'))


p_path = '18_2_id_pepdict_0215.p'
file_id_pep_dict = p.load(open(p_path,'rb'))
psm_dict_182B = p.load(open('D:/data/Naba_deep_matrisome/02152021_1/182B_psm_dict.p','rb'))
psm_dict_182A = p.load(open('D:/data/Naba_deep_matrisome/02152021_1/182A_psm_dict.p','rb'))
# get time-series accumulated covearge for ECM proteins
# df = pd.read_csv('D:/data/Naba_deep_matrisome/02152021_1/dash_info.csv')
"""
#
# df = df.drop(columns=['Unnamed: 0', 'Unnamed: 0.1'])
# df_standard = df[df['file_name']=='18_2A']
# df = df[(df['file_name']=='18_2B05')|(df['file_name']=='18_2B1')|(df['file_name']=='18_2B2')|
# (df['file_name']=='18_2B4')|(df['file_name']=='18_2B18')|(df['file_name']=='18_2B20')]
#
# df_163_standard = df[df['file_name']=='163_3A']
# df_163_timelapse = df[(df['file_name']=='163_3B05')|(df['file_name']=='163_3B1')|(df['file_name']=='163_3B2')|
# (df['file_name']=='163_3B4')|(df['file_name']=='163_3B18')|(df['file_name']=='163_3B20')]

time_lapse_file = ['18_2B05','18_2B1','18_2B2','18_2B4','18_2B18','18_2B20']

time_lapse_peptide_list = [pep for file in time_lapse_file for prot in file_id_pep_dict[file] for pep in file_id_pep_dict[file][prot] if prot in ecm_proteins]
time_lapse_psm_list = [pep+'_'+str(i) for pep in time_lapse_peptide_list for i in range(psm_dict_182B[pep])]

standard_ecm_peptide_list = [pep for prot in file_id_pep_dict['18_2A'] for pep in file_id_pep_dict['18_2A'][prot] if prot in ecm_proteins]
standard_ecm_psm_list = [pep+'_'+str(i) for pep in standard_ecm_peptide_list for i in range(psm_dict_182A[pep])]
# venn_diagram_gen({'time-lapse digestion':time_lapse_peptide_list,'20-hour standard digestion':standard_ecm_peptide_list}, title='Identified ECM peptides in 18_2')

second_p_path = '163_3_id_pepdict_0215.p'
file_id_pep_dict_2 = p.load(open(second_p_path,'rb'))
psm_dict_163B = p.load(open('D:/data/Naba_deep_matrisome/02152021_1/163B_psm_dict.p','rb'))
psm_dict_163A = p.load(open('D:/data/Naba_deep_matrisome/02152021_1/163A_psm_dict.p','rb'))

time_lapse_file2 = ['163_3B05','163_3B1','163_3B2','163_3B4','163_3B18','163_3B20']
time_lapse_peptide_list2 = [pep for file in time_lapse_file2 for prot in file_id_pep_dict_2[file] for pep in file_id_pep_dict_2[file][prot] if prot in ecm_proteins]
time_lapse_psm_list2 = [pep+'_'+str(i) for pep in time_lapse_peptide_list2 for i in range(psm_dict_163B[pep])]

standard_ecm_peptide_list2 = [pep for prot in file_id_pep_dict_2['163_3A'] for pep in file_id_pep_dict_2['163_3A'][prot] if prot in ecm_proteins]
standard_ecm_psm_list2 = [pep+'_'+str(i) for pep in standard_ecm_peptide_list2 for i in range(psm_dict_163A[pep])]
# venn_diagram_gen({'SNED1 KO time-series':time_lapse_psm_list,'':standard_ecm_psm_list,}, title='Identified ECM PSMs in SNED1 KO')
"""

# ptm analysis
from tsv_reader import peptide_phospho_reader

base_path = 'D:/data/Naba_deep_matrisome/02152021_1/18_2B18/peptide.tsv'
ecm_peptide = [pep for prot in file_id_pep_dict['18_2B18'] for pep in file_id_pep_dict['18_2B18'][prot] if prot in ecm_proteins]
print ('number of ecm peptides:',len(set(ecm_peptide)))
ptm = -17.0265
ptm_dict = peptide_phospho_reader(base_path,mod=ptm)

ptm_dict = {pep:ptm_dict[pep] for pep in ptm_dict if pep in ecm_peptide}
ptm_peptide_dict = defaultdict(set)
for pep in ptm_dict:
    for mod in ptm_dict[pep]:
        #print (mod)
        aa = re.search('[A-Z]', mod).group(0)
        ptm = re.search('\(.+\)',mod).group(0)[1:-1]
        ptm_peptide_dict[aa+'_'+ptm].add(pep)
for ptm in ptm_peptide_dict:
    print (ptm, len(ptm_peptide_dict[ptm]))