from tsv_reader import combined_proteintsv_map, protein_info_from_combined,protein_reader,protein_info_from_fasta
from pandas import ExcelWriter
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
from glob import glob

### get a summary table
"""
folders = glob('F:/sned1_biotinalytion/03112022/seach_result/*/')
print (folders)
file_list = [each.split('\\')[-2] for each in folders]
print (file_list)
fasta_path = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
protein_dict = fasta_reader(fasta_path)
protein_info_dict = protein_info_from_fasta(fasta_path)

total_protein_set = protein_reader('F:/sned1_biotinalytion/03112022/seach_result/combined_protein.tsv')
psm_path_list = [each+'psm.tsv' for each in folders]
pep_path_list = [each+'peptide.tsv' for each in folders]

file_protein_cov_dict, file_id_pep_dict,file_unique_id_pep_dict, file_prot_spec_dict,file_unique_id_pep_count_dict,file_id_pep_count_dict  = dash_dataframe(pep_path_list,psm_path_list,protein_dict,psm_path_list,psm_path_list,psm_path_list)

column =['gene']+[file+'_'+i for file in file_list for i in ['total_spec', 'total peptides identified', 'total peptides count','coverage','unique peptides identified', 'unique peptides count']]
df_info = pd.DataFrame(columns=column,index=list(total_protein_set))

for prot in total_protein_set:
    df_info.at[prot,'gene'] = protein_info_dict[prot][0]

    for file in file_list:
        for i,j in zip(['total_spec', 'total peptides identified', 'total peptides count','coverage','unique peptides identified', 'unique peptides count'],
                       [file_prot_spec_dict,file_id_pep_dict,file_id_pep_count_dict,file_protein_cov_dict, file_unique_id_pep_dict, file_unique_id_pep_count_dict]):
            if prot in file_prot_spec_dict[file]:
                # print (prot,file)
                df_info.at[prot,file+'_'+i]=j[file][prot]
            else:
                df_info.at[prot,file+'_'+i]=0
df_info.to_excel('F:/sned1_biotinalytion/03112022/seach_result/summary.xlsx')
"""
### compare public data for ECM non fractionation and 15 fractionation
matrisome_nonrepeat = 'D:/data/Naba_deep_matrisome/matrisome coverage_norepeat.xlsx'
matrisome_protein_list = pd.read_excel(matrisome_nonrepeat,index_col=1).index.tolist()

ecm_15frac_df = pd.read_csv('F:/proposal/public_data/ECM_15frac/protein.tsv',sep='\t',index_col=['Protein ID'])
ecm_15_df_filtered = ecm_15frac_df[ecm_15frac_df.index.isin(matrisome_protein_list)]
ecm_15_fra_protein = ecm_15_df_filtered.index.tolist()

ecm_nonfrac_df = pd.read_csv('F:/proposal/public_data/ECM_nonfrac/protein.tsv',sep='\t',index_col=['Protein ID'])
ecm_nonfrac_df_filtered = ecm_nonfrac_df[ecm_nonfrac_df.index.isin(matrisome_protein_list)]
ecm_nonfrac_protein_list = ecm_nonfrac_df_filtered.index.tolist()

total_ecm_list = set(ecm_15_fra_protein+ecm_nonfrac_protein_list)

df_plot = pd.DataFrame(columns=['coverage_in_nonfrac','coverage_in_15frac'], index=list(total_ecm_list))

for each in total_ecm_list:
    df_plot.loc[each,'coverage_in_nonfrac'] = ecm_nonfrac_df_filtered.loc[each,'Percent Coverage'] if each in ecm_nonfrac_protein_list else 0
    df_plot.loc[each,'coverage_in_15frac'] = ecm_15_df_filtered.loc[each,'Percent Coverage'] if each in ecm_15_fra_protein else 0

df_plot_sort = df_plot.sort_values(by=['coverage_in_nonfrac'])
print (df_plot_sort.columns)
fig,ax = plt.subplots()
x = range(df_plot_sort.shape[0])
y_nonfrac = df_plot_sort['coverage_in_nonfrac'].tolist()
y_15frac = df_plot_sort['coverage_in_15frac'].tolist()

# ax.set_ylim(0,100)
# ax.scatter(x,y_nonfrac,c='blue',s=10,label='non fractionation')
# ax.scatter(x,y_15frac,c='red',s=10,label='15 fractionations')
# ax.set_xlabel('ECM protein groups',fontsize=15)
# ax.set_ylabel('Sequence coverage %',fontsize=15)
# ax.legend()
df_new_plot = pd.DataFrame()
df_new_plot['sample_type'] = ['coverage_in_nonfrac']*df_plot_sort.shape[0]+['coverage_in_15frac']*df_plot_sort.shape[0]
df_new_plot['coverage'] = df_plot_sort['coverage_in_nonfrac'].tolist()+df_plot_sort['coverage_in_15frac'].tolist()
print (df_new_plot)

custom_palette = {'coverage_in_nonfrac':'blue','coverage_in_15frac':'red'}
sns.boxplot(data=df_new_plot,x='sample_type',y='coverage',linewidth=2.5,palette=custom_palette)
plt.show()