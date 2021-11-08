from tsv_reader import combined_proteintsv_map, protein_info_from_combined,protein_reader,protein_info_from_fasta,peptide_counting, psm_reader
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
base_path = 'D:/data/native_protein_digestion/11052021/search_result/'
folders = [base_path+folder for folder in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, folder))]

psm_path_list = [each+'/psm.tsv' for each in folders]
pep_path_list = [each+'/peptide.tsv' for each in folders]
print (pep_path_list)

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
"""
df_ecm = pd.read_excel('D:/data/Naba_deep_matrisome/matrisome coverage_norepeat.xlsx')
df_native = pd.read_excel('D:/data/native_protein_digestion/dialysis_casset_aggre_cov.xlsx',index_col=0)
protein_list = df_native.index.values
ecm_protein_list = df_ecm['protein_id'].values

for i in protein_list:
    if i in ecm_protein_list:
        print (i)
"""

def protein_cov_reader_from_tsv(protein_tsv):
    """

    :param protein_tsv:
    :return: a dictionary that stores protein and coverage from protein.tsv
    """
    protein_cov_dict = {}
    with open(protein_tsv,'r') as f:
        next(f)
        for line in f:
            line_split = line.split('\t')
            protein_cov_dict[line_split[3]] = float(line_split[7])

    return protein_cov_dict


### missed cleavage analysis
"""
from calculations_and_plot import miss_cleavage_identify,num_miss_identify
from glob import glob
from tsv_reader import get_unique_peptide
base_path = 'D:/data/native_protein_digestion/11052021/search_result/'
folder_path = [base_path + folder+'/peptide.tsv' for folder in os.listdir(base_path)
               if os.path.isdir(os.path.join(base_path, folder))]
print (folder_path)
unique_pep_dict = get_unique_peptide(folder_path)
for each in unique_pep_dict:
    # peptide_list = peptide_counting(each)
    peptide_list = unique_pep_dict[each]
    miss_cleav = miss_cleavage_identify(peptide_list)
    bool_array = [v for v in miss_cleav.values()]
    num_miss_dict = num_miss_identify(peptide_list)
    print (each,len(peptide_list))
    for each in num_miss_dict:

        print ('%i miss cleavage: ratio %f' %(each,len(num_miss_dict[each])/len(peptide_list)))
"""

### compare coverage of different proteases
"""
tryp_cov_dict = protein_cov_reader_from_tsv('D:/data/deep_proteome/different_protease/Tryp/protein1.tsv')
thermo_cov_dict = protein_cov_reader_from_tsv('D:/data/deep_proteome/different_protease/thermolysin_30min/protein.tsv')
tryp_thermo_cov_dict = protein_cov_reader_from_tsv('D:/data/deep_proteome/different_protease/tryp_30_thermo_30/protein.tsv')

for each in tryp_thermo_cov_dict:
    if each in tryp_cov_dict and each in thermo_cov_dict:
        if tryp_thermo_cov_dict[each]>tryp_cov_dict[each] and tryp_thermo_cov_dict[each]>thermo_cov_dict[each]:
            print (each, tryp_cov_dict[each],thermo_cov_dict[each],tryp_thermo_cov_dict[each])
"""


df_native = pd.read_excel('D:/data/native_protein_digestion/native_digest_cassette_aggrecov.xlsx',index_col=0)
protein_list = df_native.index
df_raw = pd.read_excel('D:/data/native_protein_digestion/raw_result.xlsx',index_col=0)
columns = [x[0] for x in os.walk('D:/data/native_protein_digestion/11052021/search_result/')][1:]
print (columns)
# proteid_ids_dict = {time:df_raw.loc[df_raw[time+'_1_native_total_spec']!=0].index.tolist() for time in columns}
peptide_ids_dict = {time:peptide_counting(time+'/peptide.tsv') for time in columns}
combined_prot_tsv = 'D:/data/native_protein_digestion/11052021/search_result/combined_protein.tsv'
combined_prot_dict = combined_proteintsv_map(combined_prot_tsv)


### heatmap
"""
fig,ax = plt.subplots(1,1, figsize=(8,15))

df_heatmap = [[df_raw.at[each,'1h_1_native_coverage'],
               df_raw.at[each,'2h_1_native_coverage'],
               df_raw.at[each,'4h_1_native_coverage'],
               df_raw.at[each,'18h_1_native_coverage']] for each in protein_list]
# df_heatmap = df_native.iloc[:,-4:]
g = sns.heatmap(data=df_heatmap,ax=ax,
                cbar_kws={'label': 'absolute coverage','shrink': 0.5},cmap="coolwarm",yticklabels=False)
ax.set_xticklabels(labels=columns, rotation=30, fontsize=10, ha='right')
plt.savefig('D:/data/native_protein_digestion/dialysis_cassette_absolute_heatmap.png', dpi=300)
plt.show()
"""

### venn diagram
"""
from tsv_reader import venn_diagram_gen2
venn_diagram_gen2(peptide_ids_dict,png_output='D:/data/native_protein_digestion/figures/cassette_peptide_id_compare.png')
"""

### coverage calculation
"""
from multiprocessing_naive_algorithym import extract_UNID_and_seq, creat_total_seq_line, zero_line_for_seq
from aho_corasick import automaton_trie, automaton_matching
from calculations_and_plot import whole_proteome_cov

# total_protein_list = set([protein_id for time in proteid_ids_dict for protein_id in proteid_ids_dict[time]])
# proteine_sub_dict = {each:protein_dict[each] for each in total_protein_list}

proteine_sub_dict = {each:protein_dict[each] for each in protein_reader(combined_prot_tsv)}
print (len(proteine_sub_dict))
id_list, seq_list = extract_UNID_and_seq(proteine_sub_dict)
seq_line = creat_total_seq_line(seq_list)

# for time in columns:
#     pep_list = peptide_ids_dict[time]
#     aho_result = automaton_matching(automaton_trie(pep_list),seq_line)
#     proteome_cov = whole_proteome_cov(aho_result,seq_line)
#     print (time, proteome_cov)

for i in range(len(columns)):
    time_list = columns[:i+1]
    agg_peplist = [pep for each in time_list for pep in peptide_ids_dict[each]]
    aho_result = automaton_matching(automaton_trie(agg_peplist),seq_line)
    proteome_cov = whole_proteome_cov(aho_result,seq_line)
    print (time_list[-1], proteome_cov)
"""

### mass distribution

from tsv_reader import pep_mass_dist
from calculations_and_plot import protein_mass
import pickle
insilico_pep_dict = pickle.load(open('D:/data/native_protein_digestion/inslico_digest_human_fasta.p','rb'))
peptide_set = {pep for v in insilico_pep_dict.values() for pep in v}
# peptide_set = peptide_counting('D:/data/deep_proteome/20200915_tryp_37C_240min/peptide.tsv')

print (len(peptide_set))
mass_list = [protein_mass(each) for each in peptide_set]
sns.kdeplot(mass_list,linewidth=2, alpha=.5)
plt.axvline(x=3500, color='k', linestyle='--')
plt.xlabel('Mass in Da')
plt.ylabel('frequency')
plt.show()
# pep_mass_dist(pep_path_list,plot='kde')
"""
### covered distance analysis

fig,axs = plt.subplots(1,1, figsize=(10,8))

df = pd.read_excel('D:/data/native_protein_digestion/11052021/cov_distance_fillw_previous.xlsx',index_col=0)
# df = df.T.ffill().bfill()
df = df.dropna()
columns = list(df.columns)
idx = list(df.index)
new_columns = [str(int(columns[i][:4]))+'_'+str(int(columns[i+1][:4]))+'min' for i in range(len(columns)-1)]

print (new_columns)
#
df_new = pd.DataFrame(index=df.index, columns=new_columns)
for i in range(df.shape[0]):
    for j in range(df.shape[1]-1):
        df_new.at[idx[i],new_columns[j]] = (df.iloc[i,j+1]-df.iloc[i,j])/(int(columns[j+1][:4])-int(columns[j][:4]))
# df_new.to_excel('D:/data/native_protein_digestion/11052021/distance_change_rate.xlsx')
df_new = df_new.astype(float)

num_prot = df.shape[0]
for each in df:
    print (each, df[each].mean())
x = range(1,13)
#
#line plot
# for tp in df_new.itertuples(index=False):
#     axs.plot(x,[i for i in tp],linestyle='-',alpha=0.8)
# axs.set_xticks(x)
# axs.set_xticklabels(list(df_new.columns), fontsize=12,ha="center", rotation=45)

#heatmap
# g = sns.clustermap(data=df_new,cbar_kws={'label': 'distance to protein centroid','shrink': 0.5},
#                cmap="YlGnBu",yticklabels=False)
# g.set_xticklabels(labels=list(df_new.columns), ha='right',fontsize=10,rotation = 45)


### violin plot
# time_point_list = [each for each in df.columns for i in range(num_prot)]
# dist_list = [v for each in df.columns for v in df[each]]
# df_violin = pd.DataFrame(dict(timepoints=time_point_list,
#                               dist=dist_list))
# sns.violinplot(data=df_violin,x='timepoints',y='dist',ax =axs)
# plt.show()
"""