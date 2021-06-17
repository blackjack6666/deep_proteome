from tsv_reader import combined_proteintsv_map, protein_info_from_combined,protein_reader,protein_info_from_fasta
from pandas import ExcelWriter
import pandas as pd
from dash_dataframe_naba import dash_dataframe
import os
from protein_coverage import fasta_reader
from parameters import protein_mass_calculator
import numpy as np


fasta_path = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
protein_dict = fasta_reader(fasta_path)

# combine B and C sample and output a single spreadsheet

"""
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

"""

# get time-series aggregated data
"""
df = pd.read_excel('D:/data/Naba_deep_matrisome/05142021_secondsearch/6_14_summary_B_C_Xinhao.xlsx',index_col=0)  # manually delete
samples = ['KO_B','SNED_B','KO_C', 'SNED_C']
time_points = ['0o5', '2', '4', '18']
df_aggregated = pd.DataFrame()

for each in df.index:
    print (each)
    protein_seq = protein_dict[each]
    gene = df.at[each,'gene']
    MW = df.at[each,'MW']
    df_aggregated.at[each, 'gene'] = gene
    df_aggregated.at[each, 'MW_kDa'] = MW

    for sample in samples:

        # aggreted_peptide_set = set([pep for time in time_points for pep in df.at[each,sample+'_'+time+'_total peptides identified']
        #                             if df.at[each,sample+'_'+time+'_total peptides identified']!=0])

        for ind, time in enumerate(time_points):
            np_array = np.zeros(len(protein_seq))
            indx_list = time_points[:ind+1]
            aggreted_peptide_set = set()
            for each_time in indx_list:
                total_peptides = df.at[each,sample+'_'+each_time+'_total peptides identified']
                if total_peptides!=0:
                    if ', ' in total_peptides:
                        for pep in total_peptides.split(', '):
                            aggreted_peptide_set.add(pep)
                    else:
                        aggreted_peptide_set.add(total_peptides)
        # print (aggreted_peptide_set)
        #     aggregated_unique_pep_count = sum([df.at[each,sample+'_'+time+'_unique peptides count'] for time in indx_list])
            aggre_unique_pepset = set()
            for each_time in indx_list:
                unique_peptides = df.at[each,sample+'_'+each_time+'_unique peptides identified']
                if unique_peptides!='set()' and unique_peptides!=0:
                    if ', ' in unique_peptides:
                        for pep in unique_peptides.split(', '):
                            aggre_unique_pepset.add(pep)
                    else:
                        aggre_unique_pepset.add(unique_peptides)
            aggregated_unique_pep_count = len(aggre_unique_pepset)
            aggregated_pep_count = len(aggreted_peptide_set)
            for peptide in aggreted_peptide_set:
                pep_loc = protein_seq.find(peptide)
                pep_end_loc = pep_loc+len(peptide)
                np_array[pep_loc:pep_end_loc]+=1
            aggregated_cov = np.count_nonzero(np_array)/len(np_array)*100

            df_aggregated.at[each, sample+'_'+time+'_aggre_coverage'] = aggregated_cov
            df_aggregated.at[each, sample+'_'+time+'_aggre_unique_pep_count'] = aggregated_unique_pep_count

df_aggregated.to_excel('D:/data/Naba_deep_matrisome/05142021_secondsearch/6_15_summary_aggregated_B_C.xlsx')

"""

# get ECM category info with gene

df_ecm = pd.read_excel('D:/data/Naba_deep_matrisome/matrisome coverage_norepeat.xlsx')
ecm_gene_category_dict = {gene:category for gene,category in zip(df_ecm['gene_id'], df_ecm['category'])}
"""
df_summary = pd.read_excel('D:/data/Naba_deep_matrisome/05142021_secondsearch/6_15_summary_aggregated_B_C.xlsx',index_col=0)
df_summary_slice = df_summary[df_summary['gene'].isin(ecm_gene_category_dict)]
"""

# average coverage value from
df_summary_slice=pd.read_excel('D:/data/Naba_deep_matrisome/05142021_secondsearch/6_15_ecm_aggregated_B_C.xlsx', index_col=0)
samples, replicates, times = ['KO','SNED'], ['B', 'C'], ['0o5','2','4','18']

df_aggre_coverage = pd.DataFrame()
for prot in df_summary_slice.index:
    gene = df_summary_slice.at[prot,'gene']
    df_aggre_coverage.at[prot,'gene'] = gene
    df_aggre_coverage.at[prot, 'category'] = ecm_gene_category_dict[gene]
    for sample in samples:
        for time in times:
            aver_cov = np.mean([df_summary_slice.at[prot,sample+'_B_'+time+'_aggre_coverage'],df_summary_slice.at[prot,sample+'_C_'+time+'_aggre_coverage']])
            df_aggre_coverage.at[prot,sample+'_'+time+'_ave_aggre_cov']=aver_cov

df_aggre_coverage.to_excel('D:/data/Naba_deep_matrisome/05142021_secondsearch/6_15_ecm_aggregated_B_C_average.xlsx')