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
from glob import glob


fasta_path = 'D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_only.fasta'
protein_dict = fasta_reader(fasta_path)
### process tsv files
# base_path = 'D:/data/native_protein_digestion/11052021/search_result/'
# folders = [base_path+folder for folder in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, folder))]
#
# psm_path_list = [each+'/psm.tsv' for each in folders]
# pep_path_list = [each+'/peptide.tsv' for each in folders]
# print (pep_path_list)

### get raw summary report
"""
protein_info_dict = protein_info_from_fasta(fasta_path)

total_protein_set = protein_reader('D:/data/native_protein_digestion/11052021/search_result/combined_protein.tsv')

base_path = 'D:/data/native_protein_digestion/11052021/search_result/'
folders = [base_path+folder for folder in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, folder))]

print (folders)
psm_path_list = [each+'/psm.tsv' for each in folders]
pep_path_list = [each+'/peptide.tsv' for each in folders]

file_list = [each.split('\\')[-2] for each in glob(base_path+'*/')]

file_protein_cov_dict, file_id_pep_dict,file_unique_id_pep_dict, file_prot_spec_dict,file_unique_id_pep_count_dict,file_id_pep_count_dict  = dash_dataframe(pep_path_list,psm_path_list,protein_dict,psm_path_list,psm_path_list,psm_path_list)
file_protein_info_dict = {file:protein_info_dict for file in file_list}

protein_mass_dict = protein_mass_calculator(total_protein_set,protein_dict)

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
df_info.to_excel('D:/data/native_protein_digestion/11052021/raw_result.xlsx')
"""
### aggregate coverage from dialysis cassette digestion

# df = pd.read_excel('D:/data/native_protein_digestion/11052021/raw_result.xlsx',index_col=0)  # manually delete
# base_path = 'D:/data/native_protein_digestion/11052021/search_result/'
# time_points = [each.split('\\')[-2] for each in glob(base_path+'*/')]
"""
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
            total_peptides = df.at[each, each_time + '_total peptides identified']
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
        #     if unique_peptides !=  'set()' and unique_peptides != 0:
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

        df_aggregated.at[each, time + '_aggre_coverage'] = aggregated_cov
        # df_aggregated.at[each, sample + '_' + time + '_aggre_unique_pep_count'] = aggregated_unique_pep_count

df_aggregated.to_excel('D:/data/native_protein_digestion/11182021/search_result_XS/aggre_cov.xlsx')
"""

### get unique coverage at each time point (calculated by mapping peptides only showed at current time point, filter
### peptides from previous time points)
"""
df_unique = pd.DataFrame()
protein_aggre_peptide_dict = {}
for each in df.index:
    print (each)
    protein_seq = protein_dict[each]
    gene = df.at[each,'gene']
    MW = df.at[each,'MW']
    df_unique.at[each, 'gene'] = gene
    df_unique.at[each, 'MW_kDa'] = MW

    all_peptides = set()
    for ind, time in enumerate(time_points):
        np_array = np.zeros(len(protein_seq))
        peptide_list = df.at[each, time+'_total peptides identified']
        # print (peptide_list)
        if peptide_list!=0:
            if ', ' in peptide_list:
                peptide_list = peptide_list.split(', ')
            else:
                peptide_list = [peptide_list]

            unique_peptides = [each for each in peptide_list if each not in all_peptides]
            # print (time,unique_peptides)
            for i in unique_peptides:
                all_peptides.add(i)
            if unique_peptides != []:
                for pep in unique_peptides:
                    pep_loc = protein_seq.find(pep)
                    pep_end_loc = pep_loc + len(pep)
                    np_array[pep_loc:pep_end_loc] += 1
                unique_cov = np.count_nonzero(np_array) / len(np_array) * 100
                df_unique.at[each,time+'_unique_cov'] = unique_cov

            else:
                df_unique.at[each, time + '_unique_cov'] = 0
        # when no peptides identified
        else:
            df_unique.at[each,time+'_unique_cov'] = 0
df_unique.to_excel('D:/data/native_protein_digestion/11052021//unique_cov.xlsx')
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


# df_native = pd.read_excel('D:/data/native_protein_digestion/native_digest_cassette_aggrecov.xlsx',index_col=0)
# protein_list = df_native.index
# df_raw = pd.read_excel('D:/data/native_protein_digestion/raw_result.xlsx',index_col=0)
# columns = [x[0] for x in os.walk('D:/data/native_protein_digestion/11052021/search_result/')][1:]
# print (columns)
# # proteid_ids_dict = {time:df_raw.loc[df_raw[time+'_1_native_total_spec']!=0].index.tolist() for time in columns}
# peptide_ids_dict = {time:peptide_counting(time+'/peptide.tsv') for time in columns}
# combined_prot_tsv = 'D:/data/native_protein_digestion/11052021/search_result/combined_protein.tsv'
# combined_prot_dict = combined_proteintsv_map(combined_prot_tsv)


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
"""
from tsv_reader import pep_mass_dist
from calculations_and_plot import protein_mass
import pickle
# insilico_pep_dict = pickle.load(open('D:/data/native_protein_digestion/inslico_digest_human_fasta.p','rb'))
# peptide_set = {pep for v in insilico_pep_dict.values() for pep in v}
peptide_set = peptide_counting('D:/data/deep_proteome/20200915_tryp_37C_0240min/peptide.tsv')
protein_set = protein_reader('D:/data/deep_proteome/20200915_tryp_37C_0240min/protein.tsv')
total_prt_dict = fasta_reader('D:/data/proteome_fasta/uniprot-proteome_UP000005640.fasta')

print (len(peptide_set))
mass_list = [protein_mass(each) for each in peptide_set]
protein_mass_list = [protein_mass(total_prt_dict[each]) for each in protein_set if each in total_prt_dict]
fig,ax = plt.subplots(1,1)
sns.kdeplot(ax=ax,x=mass_list,linewidth=2, alpha=.5,log_scale=True)
sns.kdeplot(ax=ax,x=protein_mass_list,linewidth=2,alpha=.5,log_scale=True)
plt.axvline(x=3500, color='k', linestyle='--')
plt.xlabel('Mass in Da')
plt.ylabel('frequency')
plt.text(350, 2.0,'Peptide',color='k',fontsize='15')
plt.text(70000,1.2,'Protein',color='k',fontsize='15')
plt.text(3300,-0.15,'3.5kDa',color='k')
plt.show()
# pep_mass_dist(pep_path_list,plot='kde')
"""

### covered distance analysis
"""
import pymannkendall as mk
fig,axs = plt.subplots(1,1, figsize=(10,8))

df = pd.read_excel('D:/data/native_protein_digestion/12072021/control/cov_dist_unique.xlsx',index_col=0)
# df = df.T.ffill().bfill()
# print (df)
# df = df.dropna()
# print(df.mean())
mk_result = mk.original_test(df.median().tolist())
# df = df.fillna(df.median())
print(mk_result)
# df = df.T


# columns = list(df.columns)
#
# idx = list(df.index)
# new_columns = [str(int(columns[i][:4]))+'_'+str(int(columns[i+1][:4]))+'min' for i in range(len(columns)-1)]
#
# print (new_columns)
# ### calculate first derivative of cov distance
#
# df_new = pd.DataFrame(index=df.index, columns=new_columns)
# for i in range(df.shape[0]):
#     for j in range(df.shape[1]-1):
#         df_new.at[idx[i],new_columns[j]] = (df.iloc[i,j+1]-df.iloc[i,j])/(int(columns[j+1][:4])-int(columns[j][:4]))
# # df_new.to_excel('D:/data/native_protein_digestion/11052021/distance_change_rate.xlsx')
# df_new = df_new.astype(float)
#
#
# num_prot = df.shape[0]
# for each in df:
#     print (each, df[each].mean())
# x = range(1,len(new_columns)+1)

# df_plot = pd.DataFrame(dict(time=list(range(1,len(df.columns)+1))*df.shape[0], cov_KR_dens=df.to_numpy().flatten()))

# sns.regplot(x='time',y='cov_KR_dens',data=df_plot,color='k')
# sns.boxplot(x='time',y='cov_KR_dens',data=df_plot,linewidth=2.5)
# sns.kdeplot(data=df_plot, x="cov_plddt", hue="time",legend=False)
# add_stat_annotation(axs,data=df_plot, x='time',y='cov_KR_dens',box_pairs=[(1,2)],test='Wilcoxon',
#                     text_format='star',loc='inside', verbose=2)
# plt.xlim([0,8])
# axs.set_xticks(range(1,8))
# axs.set_xticklabels(list(df.columns), fontsize=12,ha="center", rotation=45)
# plt.legend(title='Time', loc='upper right', labels=list(df.columns))
# plt.show()

#
## line plot
x = range(1,len(df.columns)+1)
for tp in df.itertuples(index=False):
    axs.plot(x,[i for i in tp],linestyle='-',alpha=0.8)
axs.set_xticks(x)
axs.set_xticklabels(list(df.columns), fontsize=12,ha="center", rotation=45)
plt.show()

"""
#heatmap
"""
df = pd.read_excel('D:/data/native_protein_digestion/10282021/h20_cov_dist_centroid_mean_nadrop.xlsx',index_col=0)
df_new = df.iloc[:,:4]
g = sns.clustermap(data=df_new,cbar_kws={'label': 'distance to protein centroid','shrink': 0.5},
               cmap="viridis",yticklabels=False)
# g.set_xticklabels(labels=list(df_new.columns), ha='right',fontsize=10,rotation = 45)
plt.show()
"""
### violin plot
"""
# time_point_list = [each for each in df.columns for i in range(num_prot)]
# dist_list = [v for each in df.columns for v in df[each]]
# df_violin = pd.DataFrame(dict(timepoints=time_point_list,
#                               dist=dist_list))
# sns.violinplot(data=df_violin,x='timepoints',y='dist',ax =axs)
# plt.show()
"""

### swarmplot
"""
df = pd.read_excel('D:/data/native_protein_digestion/10282021/h20_cov_dist_centroid_mean_nadrop.xlsx',index_col=0)
df_swarmplot = pd.DataFrame(dict(distance=df['01h_h2o'].tolist()+df['02h_h2o'].tolist()+
                                df['04h_h2o'].tolist()+df['20h_h2o'].tolist(),
                       time_point=[1]*df.shape[0]+[2]*df.shape[0]+[4]*df.shape[0]+[20]*df.shape[0],
                       gene=df['gene'].tolist()*4,
                       description=df['description'].tolist()*4,len=df['length'].tolist()*4,
                       mass=df['mass'].tolist()*4),index=df.index.tolist()*4)

fig,ax = plt.subplots(1,1, figsize=(8,5))
colormap = {1: '#FF3F33', 2: '#36FF33', 4: '#F033FF', 20: '#33B8FF'}

g = sns.swarmplot(ax=ax,data=df_swarmplot, x='time_point', y='distance', hue='time_point',size=2, palette=colormap,
                   dodge=False)
add_stat_annotation(ax=ax, data=df_swarmplot, x='time_point', y='distance',
                        box_pairs=[(1,2),
                                   (2,4),
                                   (4,20)
                                   ],
                        test='Wilcoxon', text_format='star',loc='inside', verbose=2,comparisons_correction='bonferroni')
# plt.xticks(rotation=30)
# plt.show()
plt.savefig('D:/data/native_protein_digestion/10282021/search_result_4miss/h20/dist_swarmplot.png',dpi=300)
"""

### generate coverage excel between two replicates
"""
df_XS = pd.read_excel('D:/data/native_protein_digestion/11182021/search_result_XS/unique_cov_XS.xlsx', index_col=0)
df_RN = pd.read_excel('D:/data/native_protein_digestion/11052021/unique_cov_11052021.xlsx',index_col=0)

protein_list_XS, protein_list_RN = df_XS.index.tolist(), df_RN.index.tolist()
protein_to_check = [each for each in protein_list_XS if each in protein_list_RN]  # get overlapped proteins

xs_columns = df_XS.columns.tolist()[3:]
rn_columns = [each+'_unique_cov' for each in ['0030min','0060min','0120min','0240min','0480min','1320min','1800min']]
print (xs_columns,rn_columns)
df_plot = pd.DataFrame(columns=xs_columns+rn_columns)
for prot in protein_to_check:
    cov_xs = df_XS.loc[prot,xs_columns].tolist()
    cov_rn = df_RN.loc[prot,rn_columns].tolist()
    df_plot.loc[prot] = cov_xs+cov_rn

df_plot.to_excel('D:/data/native_protein_digestion/11182021/unique_cov_XS_1105_1108.xlsx')
"""

### generate residue-to-centroid distance between replicates
"""
df_XS = pd.read_excel('D:/data/native_protein_digestion/11182021/search_result_XS/cov_distance_each_unique_XS.xlsx',index_col=0)
df_RN = pd.read_excel('D:/data/native_protein_digestion/11052021/cov_distance_each_unique.xlsx',index_col=0)
protein_list_XS, protein_list_RN = df_XS.index.tolist(), df_RN.index.tolist()
protein_to_check = [each for each in protein_list_XS if each in protein_list_RN]

xs_columns = df_XS.columns.tolist()[1:]
rn_columns = [each for each in ['0030min','0060min','0120min','0240min','0480min','1320min','1800min']]
df_plot = pd.DataFrame(columns=xs_columns+rn_columns)
for prot in protein_to_check:
    cov_xs = df_XS.loc[prot,xs_columns].tolist()
    cov_rn = df_RN.loc[prot,rn_columns].tolist()
    df_plot.loc[prot] = cov_xs+cov_rn
df_plot.to_excel('D:/data/native_protein_digestion/11182021/residue_dist_XS_1105_1118.xlsx')
"""

### scatter plot of plddt and similarity between replicates
"""
from scipy import stats
df = pd.read_excel('D:/data/native_protein_digestion/11182021/cos_sim_plddt_xs_rn.xlsx',index_col=0)

# slope, intercept, r_value, p_value, std_err = stats.linregress(df.dropna()['Euclidean'],df.dropna()['ave_plddt'],)
# ax = sns.regplot(data=df.dropna(), x='Euclidean',y='ave_plddt',color='k',scatter_kws={'s':5},
#                  line_kws={'label':"y={0:.1f}x+{1:.1f}".format(slope,intercept)})
# ax.legend()
# sns.kdeplot(data=df.dropna(),x='Euclidean')
sns.kdeplot(data=df.dropna(),x='distance_cossim',color='k')
plt.show()
"""

### correlation plot heatmap
"""
d_corr = pd.read_excel('D:/data/native_protein_digestion/11182021/residue_dist_XS_1105_1118.xlsx',index_col=0).corr()
print (d_corr)
# mask = np.triu(np.ones_like(d_corr, dtype=bool))
mask =  np.triu(d_corr)
# Set up the matplotlib figure
f, ax = plt.subplots(figsize=(11, 9))

# Generate a custom diverging colormap
# cmap = sns.diverging_palette(230, 20, as_cmap=True)

# Draw the heatmap with the mask and correct aspect ratio
sns.heatmap(d_corr, cmap='viridis', mask=mask, annot=True,
            square=True, linewidths=.5, cbar_kws={"shrink": .5})

plt.show()
"""
### mapped cleavage density analysis

def cosine_sim_calculating(v1, v2):
    """
    calculate the cosine similarity beweeen two b/y ions binned vectors
    :param v1:
    :param v2:
    :return:
    """
    from scipy import spatial
    return 1-spatial.distance.cosine(v1,v2)
"""
df = pd.read_excel('D:/data/native_protein_digestion/11182021/search_result_XS/cov_KR_density.xlsx',index_col=0)
df_fill = df.fillna(df.mean())
mean_array = df_fill.mean().tolist()

for tp in df_fill.itertuples(index=True):
    protein = tp[0]
    kr_density = tp[1:]
    cos_sim = cosine_sim_calculating(kr_density,mean_array)
    df_fill.at[protein,'cos_sim_with_mean'] = cos_sim

df_fill.to_excel('D:/data/native_protein_digestion/11182021/search_result_XS/cov_KR_density_cosine.xlsx')
"""

### spearman correlation analysis of cleaved K/R densities between control and heat shock
"""
df_control = pd.read_excel('D:/data/native_protein_digestion/12072021/control/cov_dist_unique.xlsx',index_col=0)
# df_heatshock = pd.read_excel('D:/data/native_protein_digestion/12072021/heat_shock/cov_KR_density_heatshock.xlsx',index_col=0)
from scipy.stats import spearmanr

df_control_median = df_control.median().tolist()
# df_control_fill = df_control.fillna(df_control.median())

df_spearman = pd.DataFrame(index=df_control.index, columns=['spearman correlation', 'p value'])
for tp in df_control.itertuples():
    prot, kr_densities = tp[0], tp[1:]
    try:
        corr, p_val = spearmanr(kr_densities,df_control_median,nan_policy='omit')
        df_spearman.at[prot,'spearman correlation'] = corr
        df_spearman.at[prot,'p value'] = p_val
    except:
        df_spearman.at[prot, 'spearman correlation'] = 0
        df_spearman.at[prot,'p value'] = 0
df_spearman.to_excel('D:/data/native_protein_digestion/12072021/control/spearman_cov_dist_nofill.xlsx')

"""

df_spearman = pd.read_excel('D:/data/native_protein_digestion/12072021/control/spearman_corr_pval_nofill.xlsx',index_col=0)
df_spearman_cov_dist = pd.read_excel('D:/data/native_protein_digestion/12072021/control/spearman_cov_dist_nofill.xlsx',index_col=0)
df_spearman = df_spearman.dropna()
df_spearman_cov_dist = df_spearman_cov_dist.dropna()
# df_spearman['-log10p'] = -np.log10(df_spearman['p value'])
df_spearman_cov_dist['-log10p'] = [-np.log10(each+np.random.uniform(-0.05,0.05)) for each in df_spearman_cov_dist['p value']]
# df_spearman_cov_dist['spearman correlation'] = [each+np.random.uniform(-0.05,0.05) for each in df_spearman_cov_dist['spearman correlation']]

## scatter plot showing spearman correlation and -log10 pval
# df_spearman_cov_dist.plot.scatter(x='spearman correlation',y='-log10p',c='-log10p', colormap='viridis', s=8)
# plt.axhline(y=-np.log10(0.05), color="black", linestyle="--")
# plt.axvline(x=0, color="black", linestyle="--")
# # plt.xlim([-1,0])
# plt.show()

proteins_kr = df_spearman.loc[(df_spearman['spearman correlation']<0)&(df_spearman['p value']<0.05)].index
# proteins_distance = df_spearman_cov_dist.loc[(df_spearman_cov_dist['spearman correlation']<-0.99)&(df_spearman_cov_dist['p value']<0.05)].index

# print (len(proteins_distance))
# print ([each for each in proteins_kr if each in proteins_distance])
# print (df_spearman_cov_dist.loc['P60660',:])

### compare significant proteins with humap2.0 db
"""
humap_df = pd.read_excel('D:/data/native_protein_digestion/12072021/control/humap2.0_mapped.xlsx',index_col=0)
# humap_df = humap_df[humap_df['confidence']==1]
protein_ave_p_dict = {}
for prot in proteins_kr:
    try:
        min_pval = np.min([tp[-1]*tp[1] for tp in humap_df.itertuples() if prot in tp[4].split(" ")])
        if min_pval != np.nan:
            protein_ave_p_dict[prot] = -np.log10(min_pval)
    except ValueError:
        continue

protein_kr_positive = df_spearman.loc[(df_spearman['spearman correlation']>0)&(df_spearman['p value']<0.05)].index
# protein_distance_positive = df_spearman_cov_dist.loc[(df_spearman_cov_dist['spearman correlation']>0.99)&(df_spearman_cov_dist['p value']<0.05)].index
protein_ave_p_dict_positive = {}
for prot in protein_kr_positive:
    try:
        min_pval = np.min([tp[-1]*tp[1] for tp in humap_df.itertuples() if prot in tp[4].split(" ")])
        if min_pval != np.nan:
            protein_ave_p_dict_positive[prot] = -np.log10(min_pval)
    except ValueError:
        continue

neg_pval, pos_pval = [v for v in protein_ave_p_dict.values()], [v for v in protein_ave_p_dict_positive.values()]

# filtered_neg, filtered_pos = np.where(neg_pval>=np.quantile(neg_pval,0.75)), np.where(pos_pval>=np.quantile(pos_pval,0.75))

df_plot = pd.DataFrame(dict(direct=['neg']*len(neg_pval)+['pos']*len(pos_pval),
                            log10p=neg_pval+pos_pval))

import seaborn as sns
from statannot import add_stat_annotation
fig,ax = plt.subplots(1,1)
sns.kdeplot(data=df_plot,x='log10p',hue='direct',
            fill=True, common_norm=False, palette="Set2",
            alpha=.5, linewidth=0,
            )
# sns.boxplot(data=df_plot,x='direct',y='log10p',ax=ax)
# add_stat_annotation(ax, data=df_plot,x='direct',y='log10p',box_pairs=[("neg", "pos")],
#                     test='Mann-Whitney-gt', text_format='star',loc='inside', verbose=2)
plt.show()

# print (len(proteins_kr),len(protein_kr_positive))
# print (len(protein_ave_p_dict),len(protein_ave_p_dict_positive))
"""
