from tsv_reader import combined_proteintsv_map, protein_info_from_combined,protein_reader,protein_info_from_fasta
import pandas as pd
from dash_dataframe_naba import dash_dataframe
from collections import defaultdict
from protein_coverage import fasta_reader
from parameters import protein_mass_calculator
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import ttest_rel,ttest_1samp
import seaborn as sns
from statannot import add_stat_annotation
import pickle as ppp
import matplotlib
import os

fasta_path = 'D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_only.fasta'
protein_dict = fasta_reader(fasta_path)
### process tsv files
"""


protein_info_dict = protein_info_from_fasta(fasta_path)

total_protein_set = protein_reader('D:/data/native_protein_digestion/combined_protein.tsv')

base_path = 'D:/data/native_protein_digestion/'
folders = [base_path+folder for folder in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, folder))]

print (folders)
psm_path_list = [each+'/psm.tsv' for each in folders]
pep_path_list = [each+'/peptide.tsv' for each in folders]



file_list = ['1h_1_native','2h_1_native','4h_1_native','18h_1_native','30min_1_native_parallel','1h_1_native_parallel',
             '2h_1_native_parallel','4h_1_native_parallel','8h_1_native_parallel','18h_1_native_parallel']


file_protein_cov_dict, file_id_pep_dict,file_unique_id_pep_dict, file_prot_spec_dict,file_unique_id_pep_count_dict,file_id_pep_count_dict  = dash_dataframe(pep_path_list,psm_path_list,protein_dict,psm_path_list,psm_path_list,psm_path_list)
file_protein_info_dict = {file:protein_info_dict for file in file_list}

protein_mass_dict = protein_mass_calculator(total_protein_set,protein_dict)

# combined_protein_path = 'D:/data/Naba_deep_matrisome/05142021_secondsearch/SNEDC/combined_protein.tsv'
# file_protein_spec_dict = combined_proteintsv_map(combined_protein_path)
# protein_info_dict = protein_info_from_combined(combined_protein_path)
column =['gene','MW']+[file+'_'+i for file in file_list for i in ['total_spec', 'total peptides identified', 'total peptides count','coverage','unique peptides identified', 'unique peptides count']]

print (column)
df_info = pd.DataFrame(columns=column,index=list(total_protein_set))

for prot in total_protein_set:
    df_info.at[prot,'gene'] = protein_info_dict[prot][0]
    df_info.at[prot, 'MW'] = protein_mass_dict[prot]
    for file in file_list:
        for i,j in zip(['total_spec', 'total peptides identified', 'total peptides count','coverage','unique peptides identified', 'unique peptides count'],
                       [file_prot_spec_dict,file_id_pep_dict,file_id_pep_count_dict,file_protein_cov_dict, file_unique_id_pep_dict, file_unique_id_pep_count_dict]):
            if prot in file_prot_spec_dict[file]:
                # print (prot,file)
                df_info.at[prot,file+'_'+i]=j[file][prot]
            else:
                df_info.at[prot,file+'_'+i]=0
df_info.to_excel('D:/data/native_protein_digestion/raw_result.xlsx')

"""
### aggregate coverage from dialysis cassette digestion
"""
df = pd.read_excel('D:/data/native_protein_digestion/raw_result.xlsx',index_col=0)  # manually delete

time_points = ['1h','2h','4h','18h']

df_aggregated = pd.DataFrame()
protein_aggre_peptide_dict = {}
for each in df.index:
    print (each)
    protein_seq = protein_dict[each]
    gene = df.at[each,'gene']
    MW = df.at[each,'MW']
    df_aggregated.at[each, 'gene'] = gene
    df_aggregated.at[each, 'MW_kDa'] = MW

    for ind, time in enumerate(time_points):
        np_array = np.zeros(len(protein_seq))
        indx_list = time_points[:ind + 1]

        aggreted_peptide_set = set()
        for each_time in indx_list:
            total_peptides = df.at[each, each_time + '_1_native_total peptides identified']
            if total_peptides != 0:
                if ', ' in total_peptides:
                    for pep in total_peptides.split(', '):
                        aggreted_peptide_set.add(pep)
                else:
                    aggreted_peptide_set.add(total_peptides)
        # print (aggreted_peptide_set)
        #     aggregated_unique_pep_count = sum([df.at[each,sample+'_'+time+'_unique peptides count'] for time in indx_list])
        # aggre_unique_pepset = set()
        # for each_time in indx_list:
        #     unique_peptides = df.at[each, sample[:-1] + each_time + sample[-1] + '_total peptides identified']
        #     if unique_peptides != 'set()' and unique_peptides != 0:
        #         if ', ' in unique_peptides:
        #             for pep in unique_peptides.split(', '):
        #                 aggre_unique_pepset.add(pep)
        #         else:
        #             aggre_unique_pepset.add(unique_peptides)
        # aggregated_unique_pep_count = len(aggre_unique_pepset)
        aggregated_pep_count = len(aggreted_peptide_set)
        for peptide in aggreted_peptide_set:
            pep_loc = protein_seq.find(peptide)
            pep_end_loc = pep_loc + len(peptide)
            np_array[pep_loc:pep_end_loc] += 1
        aggregated_cov = np.count_nonzero(np_array) / len(np_array) * 100

        df_aggregated.at[each, time + '_1_native_aggre_coverage'] = aggregated_cov
        # df_aggregated.at[each, sample + '_' + time + '_aggre_unique_pep_count'] = aggregated_unique_pep_count

df_aggregated.to_excel('D:/data/native_protein_digestion/dialysis_casset_aggre_cov.xlsx')
"""
df_ecm = pd.read_excel('D:/data/Naba_deep_matrisome/matrisome coverage_norepeat.xlsx')
df_native = pd.read_excel('D:/data/native_protein_digestion/dialysis_casset_aggre_cov.xlsx',index_col=0)
protein_list = df_native.index.values
ecm_protein_list = df_ecm['protein_id'].values

for i in protein_list:
    if i in ecm_protein_list:
        print (i)