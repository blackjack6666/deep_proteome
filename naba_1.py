import pandas as pd
import pickle as p
from tsv_reader import venn_diagram_gen

matrisome_cov_csv = 'D:/data/Naba_deep_matrisome/matrisome coverage.xlsx'
df_mat = pd.read_excel(matrisome_cov_csv)
df_mat = df_mat.drop_duplicates()
ecm_proteins = df_mat['protein_id'].tolist()

p_path = '18_2_id_pepdict_0215.p'
file_id_pep_dict = p.load(open(p_path,'rb'))

# get time-series accumulated covearge for ECM proteins
# df = pd.read_csv('D:/data/Naba_deep_matrisome/02152021_1/dash_info.csv')
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

standard_ecm_peptide_list = [pep for prot in file_id_pep_dict['18_2A'] for pep in file_id_pep_dict['18_2A'][prot] if prot in ecm_proteins]

# venn_diagram_gen({'time-lapse digestion':time_lapse_peptide_list,'20-hour standard digestion':standard_ecm_peptide_list}, title='Identified ECM peptides in 18_2')

second_p_path = '163_3_id_pepdict_0215.p'
file_id_pep_dict_2 = p.load(open(second_p_path,'rb'))
time_lapse_file2 = ['163_3B05','163_3B1','163_3B2','163_3B4','163_3B18','163_3B20']
time_lapse_peptide_list2 = [pep for file in time_lapse_file2 for prot in file_id_pep_dict_2[file] for pep in file_id_pep_dict_2[file][prot] if prot in ecm_proteins]
standard_ecm_peptide_list2 = [pep for prot in file_id_pep_dict_2['163_3A'] for pep in file_id_pep_dict_2['163_3A'][prot] if prot in ecm_proteins]
venn_diagram_gen({'18_2':time_lapse_peptide_list,'163_3':time_lapse_peptide_list2}, title='Identified ECM peptides in time-lapse digestion')
