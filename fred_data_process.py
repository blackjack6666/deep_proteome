from tsv_reader import combined_proteintsv_map, protein_info_from_combined,protein_reader,protein_info_from_fasta
from pandas import ExcelWriter
import pandas as pd
from dash_dataframe_naba import dash_dataframe
import os
from protein_coverage import fasta_reader
from parameters import protein_mass_calculator

fasta_path = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
protein_dict = fasta_reader(fasta_path)
protein_info_dict = protein_info_from_fasta(fasta_path)
folder_path = 'D:/data/Naba_deep_matrisome/05142021_secondsearch/'
combined_protein_file_list = [folder_path+each+'/combined_protein.tsv' for each in ['KOB','KOC', 'SNEDB', 'SNEDC']]

total_protein_list = [prot for each in combined_protein_file_list for prot in protein_reader(each)]
total_protein_set = set(total_protein_list)

protein_mass_dict = protein_mass_calculator(total_protein_set,protein_dict)



base_path = 'D:/data/Naba_deep_matrisome/05142021_secondsearch/'
folders = [base_path+each+folder for each in ['KOB/','KOC/', 'SNEDB/', 'SNEDC/'] for folder in os.listdir(base_path+each) if '.' not in folder]
print (folders)
psm_path_list = [each+'/psm.tsv' for each in folders]
pep_path_list = [each+'/peptide.tsv' for each in folders]

file_list = ['KO_B_0o5', 'KO_B_2','KO_B_4', 'KO_B_18',
            'SNED_B_0o5', 'SNED_B_2','SNED_B_4', 'SNED_B_18',
             'KO_C_0o5', 'KO_C_2', 'KO_C_4', 'KO_C_18',
             'SNED_C_0o5', 'SNED_C_2', 'SNED_C_4', 'SNED_C_18']

file_protein_cov_dict, file_id_pep_dict,file_unique_id_pep_dict, file_prot_spec_dict,file_unique_id_pep_count_dict,file_id_pep_count_dict  = dash_dataframe(pep_path_list,psm_path_list,protein_dict,psm_path_list,psm_path_list,psm_path_list)
file_protein_info_dict = {file:protein_info_dict for file in file_list}


file_path = 'D:/data/Naba_deep_matrisome/matrisome coverage.xlsx'
df = pd.read_excel(file_path)
df = df.drop_duplicates()
ecm_prot_list = df['protein_id'].tolist()

# combined_protein_path = 'D:/data/Naba_deep_matrisome/05142021_secondsearch/SNEDC/combined_protein.tsv'
# file_protein_spec_dict = combined_proteintsv_map(combined_protein_path)
# protein_info_dict = protein_info_from_combined(combined_protein_path)
column =['gene', 'MW']+[file+'_'+i for file in file_list for i in ['total_spec', 'total peptides identified', 'total peptides count','coverage','unique peptides identified', 'unique peptides count']]

print (column)
df_info = pd.DataFrame(columns=column,index=list(total_protein_set))

for prot in total_protein_set:
    df_info.at[prot,'gene'] = protein_info_dict[prot][0]
    df_info.at[prot,'MW'] = protein_mass_dict[prot]
    for file in file_list:
        for i,j in zip(['total_spec', 'total peptides identified', 'total peptides count','coverage','unique peptides identified', 'unique peptides count'],
                       [file_prot_spec_dict,file_id_pep_dict,file_id_pep_count_dict,file_protein_cov_dict, file_unique_id_pep_dict, file_unique_id_pep_count_dict]):
            if prot in file_prot_spec_dict[file]:
                # print (prot,file)
                df_info.at[prot,file+'_'+i]=j[file][prot]
            else:
                df_info.at[prot,file+'_'+i]=0
df_info.to_excel('D:/data/Naba_deep_matrisome/05142021_secondsearch/6_14_summary_B_C_Xinhao.xlsx')



# with ExcelWriter('D:/data/Naba_deep_matrisome/05142021_secondsearch/SNEDC/mat_protein_spec_SNEDC.xlsx') as writer:
#     for file in file_protein_spec_dict:
#         if file != 'Summarized':
#             info_list = [[prot,file_protein_spec_dict[file][prot], file_protein_cov_dict[file][prot], protein_info_dict[prot][1], protein_info_dict[prot][2]]
#                          for prot in file_protein_spec_dict[file] if prot in ecm_prot_list and file_protein_cov_dict[file][prot]!=0]
#             df = pd.DataFrame(info_list, columns=['protein','total_spec_count','coverage', 'gene_name', 'description'])
#             df.to_excel(writer, '%s' % file)
