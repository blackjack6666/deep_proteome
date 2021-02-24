# temp script doing matrisome analysis

import pandas as pd
import pickle as p
from protein_coverage import fasta_reader
import numpy as np
from tsv_reader import venn_diagram_gen, combined_proteintsv_map,protein_info_from_combined,protein_info_from_fasta
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pandas import ExcelWriter


# get spec count from combined for all identified proteins
"""
combined_prot_tsv = 'D:/data/Naba_deep_matrisome/01102021/combined_protein.tsv'
info_dict = combined_proteintsv_map(combined_prot_tsv)
protein_info_dict = protein_info_from_combined(combined_prot_tsv)

with ExcelWriter('D:/data/Naba_deep_matrisome/01102021/all_protein_spec_1_29.xlsx') as writer:
        for each in info_dict:
            info_list = [[prot,info_dict[each][prot],protein_info_dict[prot][0],protein_info_dict[prot][1],protein_info_dict[prot][2]]
                         for prot in info_dict[each]]
            df = pd.DataFrame(info_list, columns=['protein id', 'spectra count', 'entry name', 'gene name', 'description'])
            df.to_excel(writer,'%s' % each)
"""
matrisome_cov_csv = 'D:/data/Naba_deep_matrisome/matrisome coverage.xlsx'
df_mat = pd.read_excel(matrisome_cov_csv)
df_mat = df_mat.drop_duplicates()


p_path = '18_2_id_pepdict_0215.p'
file_id_pep_dict = p.load(open(p_path,'rb'))

fasta_path = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
protein_dict = fasta_reader(fasta_path)
protein_info_dict = protein_info_from_fasta(fasta_path)

# get time-series accumulated covearge for ECM proteins
df = pd.read_csv('D:/data/Naba_deep_matrisome/02152021_1/dash_info.csv')

df = df.drop(columns=['Unnamed: 0', 'Unnamed: 0.1'])
df = df[(df['file_name']=='18_2B05')|(df['file_name']=='18_2B1')|(df['file_name']=='18_2B2')|
(df['file_name']=='18_2B4')|(df['file_name']=='18_2B18')|(df['file_name']=='18_2B20')|(df['file_name']=='18_2A')]

print (df.shape)

protein_id_set = df['protein_id'].unique()
# print (protein_id_set)
file_set = df['file_name'].unique()
print (df['ecm_class'].unique())


print (len(protein_id_set))

# print ([prot for prot in protein_id_set if prot not in id_pep_dict_1805])
df_18_2_accumulated = pd.DataFrame()
df_accumulated_cov = pd.DataFrame()
index_18_2 = ['18_2B05','18_2B1','18_2B2','18_2B4','18_2B18','18_2B20']
for prot in protein_id_set:
    protein_seq = protein_dict[prot]

    for idx, val in enumerate(index_18_2):
        index_list = index_18_2[:idx+1]
        np_array = np.zeros(len(protein_seq))
        # print (index_list)
        pep_list = []
        for each_file in index_list:
            if prot in file_id_pep_dict[each_file]:
                pep_list+=list(file_id_pep_dict[each_file][prot])

            else:
                pep_list+=[]

        peptide_set = set(pep_list)

        for each in peptide_set:
            pep_loc = protein_seq.find(each)
            pep_end_loc = pep_loc+len(each)
            np_array[pep_loc:pep_end_loc]+=1
        coverage = np.count_nonzero(np_array)/len(np_array)*100

        df_18_2_accumulated.loc[val,prot] = coverage
        df_accumulated_cov.loc[val,protein_info_dict[prot][0]] = coverage
print (df_18_2_accumulated)
df_accumulated_cov.to_csv('D:/data/Naba_deep_matrisome/02152021_1/163_3_accumulated_cov_gene.csv')
# line plot
average_182A_coverage = df[df['file_name']=='18_2A']['coverage'].mean()

# line plot color map
ecm_class_color_dict = {"Collagens": '#F23F51', 'ECM-affiliated Proteins':'#23AECA',
                        'ECM Regulators':"#23CA66","Secreted Factors":"#E3EB09",
                        "ECM Glycoproteins":"#EBA709", "Proteoglycans":"#EB09DC"}
color_map = {prot:ecm_class_color_dict[ecm_class] for prot,ecm_class in zip(df['protein_id'],df['ecm_class'])}

colors = [v for v in ecm_class_color_dict.values()]
colors.append('black') # average flat line
labels = [k for k in ecm_class_color_dict.keys()]
labels.append('Average ECM protein coverage in 18_2A')
lines = [Line2D([0], [0], color=c, linewidth=3, linestyle='-') for c in colors[:-1]]
lines.append(Line2D([0], [0], color='black', linewidth=3, linestyle='--'))

# plotting

# ax = df_18_2_accumulated.plot.line(color=color_map, figsize=(10,15))  # dataframe plot


fig,ax = plt.subplots(1,1, figsize=(10,15))

#plot each line at a time
# x = range(len(df_18_2_accumulated.index))
# for column in df_18_2_accumulated.columns:
#     y = df_18_2_accumulated[column]
#     if y[-1]-y[0] == 0: # if line is flat, make it transparent
#         ax.plot(x,y, color=color_map[column], linestyle='-',alpha=0.1)
#     else:
#         ax.plot(x,y, color=color_map[column], linestyle='-')
#
# # ax.get_legend().remove()
# ax.legend(lines, labels, framealpha=0.5,loc='upper right',fontsize=15)
# ax.set_xticks(np.arange(len(df_18_2_accumulated.index)))
# ax.set_xticklabels(['0.5h','1h','2h','4h','18h', '20h'], fontsize=15)
# ax.set_title('time-series accumulated ECM protein coverage in 163_3', fontsize=22)
# ax.set_xlabel('time point', fontsize=20)
# ax.set_ylabel('%coverage',fontsize=20)
# ax.tick_params(axis='both', which='major', labelsize=15)
# plt.axhline(y=average_182A_coverage,xmin=0.04, xmax=0.96,linestyle='--',color='k',linewidth=3)
# plt.tight_layout()
# plt.show()


# box-plot plotting
from scipy.stats import ttest_ind
coverage_18_2A = df[df['file_name']=='18_2A']['coverage'].tolist()
coverage_18_2B_20 = [each for each in df_18_2_accumulated.iloc[-1,:].tolist() if each!=0]
# t, p = ttest_ind(coverage_18_2A, coverage_18_2B_20, equal_var=False)
# print (t,p)
# ax.boxplot([coverage_18_2A,coverage_18_2B_20], labels=['20 hour', 'time-series accumulated 20 hour'],showfliers=False)
# # plot points onto boxplot
# for i,cov_list in zip(range(1,3),[coverage_18_2A,coverage_18_2B_20]):
#     x = np.random.normal(i, 0.04, size=len(cov_list))  # jitter x-axis
#     ax.plot(x,cov_list,'r.',markersize=10 ,alpha=0.2)
# ax.set_title('ECM protein coverage distribution of 163_3 in two digestions', fontsize=22)
# ax.tick_params(axis='both', which='major', labelsize=15)
# ax.set_ylabel('%coverage',fontsize=20)
# plt.show()

# violin plot
ax.violinplot([coverage_18_2A,coverage_18_2B_20])
ax.get_xaxis().set_tick_params(direction='out')
ax.xaxis.set_ticks_position('bottom')
ax.set_xticks(np.arange(1, len(['20 hour', 'time-lapse accumulated 20 hour']) + 1))
ax.set_xticklabels(['20 hour','timelapse accumulated 20 hour'])
ax.set_ylabel('%coverage',fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_title('ECM protein coverage distribution of 18_2 in two digestions', fontsize=22)
# ax.set_xlim(0.25, len(labels) + 0.75)
# ax.set_xlabel('Sample name')
plt.show()

