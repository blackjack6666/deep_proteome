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