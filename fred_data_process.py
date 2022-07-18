from tsv_reader import combined_proteintsv_map, protein_info_from_combined,protein_reader,protein_info_from_fasta
from pandas import ExcelWriter
# import os
# os.environ["MODIN_ENGINE"] = "ray"
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

# import ray
# ray.init()
font = {'family' : 'Arial',
        'weight' : 'medium',
        'size'   : 8}

matplotlib.rc('font', **font)
fasta_path = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
protein_dict = fasta_reader(fasta_path)

### combine B and C sample and output a single spreadsheet
"""
protein_info_dict = protein_info_from_fasta(fasta_path)
# folder_path = 'D:/data/Naba_deep_matrisome/05142021_secondsearch/'
# combined_protein_file_list = [folder_path+each+'/combined_protein.tsv' for each in ['KOB','KOC', 'SNEDB', 'SNEDC']]
#
# total_protein_list = [prot for each in combined_protein_file_list for prot in protein_reader(each)]
# total_protein_set = set(total_protein_list)
total_protein_set = protein_reader('D:/data/Naba_deep_matrisome/06272022/combined_protein.tsv')

protein_mass_dict = protein_mass_calculator(total_protein_set,protein_dict)

# base_path = 'D:/data/Naba_deep_matrisome/05142021_secondsearch/'
# folders = [base_path+each+folder for each in ['KOB/','KOC/', 'SNEDB/', 'SNEDC/'] for folder in os.listdir(base_path+each) if '.' not in folder]
base_path = 'D:/data/Naba_deep_matrisome/06272022/'
folders = [base_path+folder for folder in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, folder))]

print (folders)

psm_path_list = [each+'/psm.tsv' for each in folders]
pep_path_list = [each+'/peptide.tsv' for each in folders]


# file_list = ['KO_B_0o5', 'KO_B_2','KO_B_4', 'KO_B_18',
#             'SNED_B_0o5', 'SNED_B_2','SNED_B_4', 'SNED_B_18',
#              'KO_C_0o5', 'KO_C_2', 'KO_C_4', 'KO_C_18',
#              'SNED_C_0o5', 'SNED_C_2', 'SNED_C_4', 'SNED_C_18']

# file_list = ['GFP_seq_30D','GFP_seq_120D','GFP_seq_240D','GFP_seq_1080D','GFP_120D','GFP_1080D',
#              'SNED1_seq_30D','SNED1_seq_120D','SNED1_seq_240D','SNED1_seq_1080D','SNED1_120D','SNED1_1080D',
#              'GFP_seq_30F','GFP_seq_120F','GFP_seq_240F','GFP_seq_1080F','GFP_120F','GFP_1080F',
#              'SNED1_seq_30F','SNED1_seq_120F','SNED1_seq_240F','SNED1_seq_1080F','SNED1_120F','SNED1_1080F']

file_list = [each.split('/')[-1] for each in folders]

file_protein_cov_dict, file_id_pep_dict,file_unique_id_pep_dict, file_prot_spec_dict,file_unique_id_pep_count_dict,file_id_pep_count_dict  = dash_dataframe(pep_path_list,psm_path_list,protein_dict,psm_path_list,psm_path_list,psm_path_list)
file_protein_info_dict = {file:protein_info_dict for file in file_list}


# file_path = 'D:/data/Naba_deep_matrisome/matrisome coverage.xlsx'
# df = pd.read_excel(file_path)
# df = df.drop_duplicates()
# ecm_prot_list = df['protein_id'].tolist()

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
df_info.to_excel('D:/data/Naba_deep_matrisome/06272022/summary_GFP_20220627.xlsx')
"""
# with ExcelWriter('D:/data/Naba_deep_matrisome/05142021_secondsearch/SNEDC/mat_protein_spec_SNEDC.xlsx') as writer:
#     for file in file_protein_spec_dict:
#         if file != 'Summarized':
#             info_list = [[prot,file_protein_spec_dict[file][prot], file_protein_cov_dict[file][prot], protein_info_dict[prot][1], protein_info_dict[prot][2]]
#                          for prot in file_protein_spec_dict[file] if prot in ecm_prot_list and file_protein_cov_dict[file][prot]!=0]
#             df = pd.DataFrame(info_list, columns=['protein','total_spec_count','coverage', 'gene_name', 'description'])
#             df.to_excel(writer, '%s' % file)



### get time-series aggregated data
"""
df = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_summary_D_F.xlsx',index_col=0)  # manually delete
# samples = ['KO_B','SNED_B','KO_C', 'SNED_C']
# time_points = ['0o5', '2', '4', '18']
samples = ['GFP_seq_D', 'SNED1_seq_D', 'GFP_seq_F', 'SNED1_seq_F']
time_points = ['30','120','240','1080']

df_aggregated = pd.DataFrame()
protein_aggre_peptide_dict = {}
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
                total_peptides = df.at[each,sample[:-1]+each_time+sample[-1]+'_total peptides identified']
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
                unique_peptides = df.at[each,sample[:-1]+each_time+sample[-1]+'_total peptides identified']
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

# df_aggregated.to_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_summary_aggregated_D_F.xlsx')
"""


### get ECM category info with gene and filter ECM gene in total genes
"""
df_ecm = pd.read_excel('D:/data/Naba_deep_matrisome/matrisome coverage_norepeat.xlsx')
df_ecm_aggre = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/8_1_matrisome_average_aggre.xlsx',index_col=0)
ecm_gene_category_dict = {gene:category for gene,category in zip(df_ecm['gene_id'], df_ecm['category'])}

df_summary = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_summary_aggregated_D_F.xlsx',index_col=0)
df_summary_slice = df_summary[df_summary['gene'].isin(df_ecm_aggre.gene)]
# df_summary_slice.to_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_ecm_aggregated_D_F.xlsx')
"""

### average coverage value from biological replicates data
"""
df_summary_slice=pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_ecm_aggregated_D_F.xlsx', index_col=0)
samples, times = ['GFP_seq_', 'SNED1_seq_'],['30','120','240','1080']

df_aggre_coverage = pd.DataFrame()
for prot in df_summary_slice.index:
    gene = df_summary_slice.at[prot,'gene']
    df_aggre_coverage.at[prot,'gene'] = gene
    df_aggre_coverage.at[prot, 'category'] = ecm_gene_category_dict[gene]
    for sample in samples:
        for time in times:
            aver_cov = np.mean([df_summary_slice.at[prot,sample+'D_'+time+'_aggre_coverage'],df_summary_slice.at[prot,sample+'F_'+time+'_aggre_coverage']])
            df_aggre_coverage.at[prot,sample+time+'_ave_aggre_cov']=aver_cov

df_aggre_coverage.to_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_ecm_aggregated_D_F_average.xlsx')
"""

"""
sp_tr_dict = defaultdict(set)
with open(fasta_path,'r') as f:
    f_split = f.read().split('\n>')
    for each in f_split[1:]:
        sp_tr_dict[each.split('|')[0]].add(each.split('|')[1])
    sp_tr_dict[f_split[0].split('|')[0].split('>')[1]].add(f_split[0].split('|')[1])
df_ecm_aggre = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/8_1_matrisome_average_aggre.xlsx',index_col=0)
df_ecm_aggre = df_ecm_aggre.copy()
df_ecm_aggre = df_ecm_aggre[df_ecm_aggre.index.isin(sp_tr_dict['sp'])]
category_list = df_ecm_aggre['category']
df_cov_derivative_delta = pd.DataFrame()
gfp_1080_agg = df_ecm_aggre.loc[df_ecm_aggre['category']=='ECM-affiliated Proteins']['GFP_seq_1080_ave_aggre_cov'].tolist()
sned1_1080_agg = df_ecm_aggre.loc[df_ecm_aggre['category']=='ECM-affiliated Proteins']['SNED1_seq_1080_ave_aggre_cov'].tolist()
"""

### calculate first derivatives
"""
for each in df_ecm_aggre.itertuples():
    prot = each[0]
    category = each[2]
    df_cov_derivative_delta.at[prot,'category'] = category
    for i,j, time_diff, time_interv in zip(range(3,6),range(7,10),[90,120,840],["0.5-2h","2-4h","4-18h"]):
        first_deriv_diffence =(each[j+1]-each[j])/time_diff - (each[i+1]-each[i])/time_diff

        df_cov_derivative_delta.at[prot,time_interv]=first_deriv_diffence

df_cov_derivative_delta_sub = df_cov_derivative_delta[df_cov_derivative_delta['category']=='Secreted Factors']

print (ttest_1samp(df_cov_derivative_delta_sub['0.5-2h'],popmean=0))
print (df_cov_derivative_delta_sub.loc[df_cov_derivative_delta_sub['0.5-2h']<0])
# sns.violinplot(data=df_cov_derivative_delta, x='category',y='0.5-2h')
# plt.title('Delta of coverage first derivative between GFP-SNED1 and GFP')
# plt.xticks(rotation=30)
# plt.show()
"""

### plotting aggregated coverage
"""
# fig,ax = plt.subplots(1,1, figsize=(10,15))

# line plot color map

# ecm_class_color_dict = {"Collagens": '#F23F51', 'ECM-affiliated Proteins':'#23AECA',
#                         'ECM Regulators':"#23CA66","Secreted Factors":"#E3EB09",
#                         "ECM Glycoproteins":"#EBA709", "Proteoglycans":"#EB09DC"}
ecm_class_color_dict = {"Collagens": '#0584B7', 'ECM-affiliated Proteins':'#F4511E',
                        'ECM Regulators':"#F9A287","Secreted Factors":"#FFE188",
                        "ECM Glycoproteins":"#133463", "Proteoglycans":"#59D8E6"}
color_map = {prot:ecm_class_color_dict[ecm_class] for prot,ecm_class in zip(df_ecm_aggre.index,df_ecm_aggre['category'])}

colors = [v for v in ecm_class_color_dict.values()]
colors.append('black') # average flat line
labels = [k for k in ecm_class_color_dict.keys()]
labels.append('Average ECM protein coverage in aggre_18h')
lines = [Line2D([0], [0], color=c, linewidth=3, linestyle='-') for c in colors[:-1]]
lines.append(Line2D([0], [0], color='black', linewidth=3, linestyle='--'))
"""

### line plot of time_series coverage/first derivates
"""
#plot each line at a time
fig,ax = plt.subplots(1,1, figsize=(10,15))
x = range(4)
for each in df_ecm_aggre.itertuples():
    # print (each)
    column = each[0]
    y = each[-4:]
    if y[-1]-y[0] == 0: # if line is flat, make it transparent
        ax.plot(x,y, color=color_map[column], linestyle='-',alpha=0.1)
    else:
        ax.plot(x,y, color=color_map[column], linestyle='-')
    # ax.plot(x,y, color=color_map[column], linestyle='-') # first derivates

# ax.get_legend().remove()
ax.legend(lines, labels, framealpha=0.5,loc='upper right',fontsize=15)
ax.set_xticks(x)
# ax.set_xticklabels(['0.5-2h','2-4h','4-18h'], fontsize=15) # first derivates
ax.set_xticklabels(['0.5h','2h','4h','18h'])
ax.set_title('ECM coverage in sequential GFP-SNED1', fontsize=22)
ax.set_xlabel('Time', fontsize=20)
#ax.set_ylabel('coverage_change/min',fontsize=20) # first derivates
ax.set_ylabel('Coverage%', fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=15)
plt.axhline(y=df_ecm_aggre['SNED1_seq_1080_ave_aggre_cov'].mean(),xmin=0.04, xmax=0.96,linestyle='--',color='k',linewidth=3)
plt.tight_layout()
plt.show()
"""

### break down ECM into categories and plot line
"""
from math import log10
fig,axs = plt.subplots(6,1, figsize=(10,20))
x = [log10(0.5), log10(2), log10(4), log10(18)]
for each_cat,ax in zip(ecm_class_color_dict,axs):
    sub_df = df_ecm_aggre[df_ecm_aggre['category']==each_cat]
    for each in sub_df.itertuples():
        y = np.array(each[4:8])-np.array(each[8:12])
        y_log10 = np.log10(y)
        ax.plot(x, y, color='lightcoral', linestyle='-')
        # if y[-1]-y[0] == 0: # if line is flat, make it transparent
        #     ax.plot(x,y_log10, color='lightcoral', linestyle='-',alpha=0.1)
        # else:
        #     ax.plot(x,y_log10, color='lightcoral', linestyle='-')
    ax.set_xticks([log10(0.5), log10(2), log10(4), log10(18)])
    ax.set_xticklabels(['0.5h','2h','4h','18h'], fontsize=10)
    ax.set_title(each_cat, fontsize=12)
    ax.set_xlabel('log time point', fontsize=12)
    ax.set_ylabel('log %coverage',fontsize=12)
    ax.set(yticklabels=[])
plt.show()
"""

"""
df_summary = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_summary_D_F_squential_standard.xlsx',index_col=0)
normal18GFP_cov = [np.mean([df_summary.at[prot,'GFP_1080D_coverage'], df_summary.at[prot,'GFP_1080F_coverage']]) for prot in df_ecm_aggre.index]
normal2GFP_cov = [np.mean([df_summary.at[prot,'GFP_120D_coverage'], df_summary.at[prot,'GFP_120F_coverage']]) for prot in df_ecm_aggre.index]
normal18SNED_cov = [np.mean([df_summary.at[prot,'SNED1_1080D_coverage'], df_summary.at[prot,'SNED1_1080F_coverage']]) for prot in df_ecm_aggre.index]
normal2SNED_cov = [np.mean([df_summary.at[prot,'SNED1_120D_coverage'], df_summary.at[prot,'SNED1_120F_coverage']]) for prot in df_ecm_aggre.index]

aggre_18_cov = df_ecm_aggre['SNED1_seq_1080_ave_aggre_cov'].tolist()
aggre_ko18_cov = df_ecm_aggre['GFP_seq_1080_ave_aggre_cov'].tolist()
aggre_2_cov = df_ecm_aggre['SNED1_seq_120_ave_aggre_cov'].tolist()
aggre_ko2_cov = df_ecm_aggre['GFP_seq_120_ave_aggre_cov'].tolist()

print (np.mean(aggre_ko18_cov), np.mean(aggre_18_cov))
print (ttest_rel(normal2GFP_cov,normal2SNED_cov,alternative='greater'))

df_new = pd.DataFrame(dict(gene=df_ecm_aggre.gene,category=df_ecm_aggre.category,GFP18_normal=normal18GFP_cov,
                           SNED18_normal=normal18SNED_cov,diff=np.array(normal18GFP_cov)-np.array(normal18SNED_cov)),index=df_ecm_aggre.index)
print (df_new.sort_values(by='diff',ascending=False))
# category_list = df_ecm_aggre['category'].tolist()
"""

### violin plot
"""
df_violin = pd.DataFrame(dict(category=category_list, cov=normal18GFP_cov,sample=['GFP 18h']*len(category_list))).append(pd.DataFrame(dict(category=category_list, cov=normal18SNED_cov,sample=['GFP-SNED1 18h']*len(category_list))))
fig,ax = plt.subplots(1,1, figsize=(8,5))
ax = sns.violinplot(data=df_violin, x='category', y='cov', hue='sample', palette='Set2', split=False)
ax = sns.swarmplot(data=df_violin, x='category', y='cov', hue='sample',size=4, palette='rocket', dodge=True)
add_stat_annotation(ax,data=df_violin,x='category',y='cov',hue='sample',
                    box_pairs=[(("ECM Glycoproteins","GFP 18h"),("ECM Glycoproteins","GFP-SNED1 18h")),
                               (("ECM-affiliated Proteins","GFP 18h"),("ECM-affiliated Proteins","GFP-SNED1 18h")),
                               (("Collagens", "GFP 18h"),("Collagens","GFP-SNED1 18h")),
                               (("Proteoglycans", "GFP 18h"),("Proteoglycans","GFP-SNED1 18h")),
                               (("ECM Regulators", "GFP 18h"),("ECM Regulators","GFP-SNED1 18h")),
                               (("Secreted Factors", "GFP 18h"),("Secreted Factors", "GFP-SNED1 18h"))],test='t-test_paired',text_format='star',loc='outside',verbose=2, comparisons_correction=None)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:2], labels[:2])
plt.xticks(rotation=30)
plt.show()
"""

### heatmap/cluster map
"""
# columns = ['standard_2h_GFP','standard_2h_SNED1','aggre_2h_GFP','aggre_2h_SNED1','standard_18h_GFP','standard_18h_SNED1','aggre_18h_GFP','aggre_18h_SNED1']
columns = ['30min', '2h','4h','18h']
# df_heatmap = pd.DataFrame(dict(standard_2h_GFP=normal2GFP_cov, standard_2h_SNED1=normal2SNED_cov, aggre_2h_GFP=aggre_ko2_cov,aggre_2h_SNED1=aggre_2_cov,
#                                standard_18h_GFP=normal18GFP_cov, standard_18h_SNED1=normal18SNED_cov, aggre_18h_GFP=aggre_ko18_cov, aggre_18h_SNED1=aggre_18_cov),
#                            index=df_ecm_aggre.index, columns=columns)
df_sub = df_ecm_aggre.sort_values(by=['category']).iloc[:, 3:11]
df_heatmap = pd.DataFrame([np.array(row[:4])-np.array(row[-4:]) for row in df_sub.itertuples(index=False)], columns=columns, index=df_sub.index)
# print (df_heatmap.head)

fig,ax = plt.subplots(1,1, figsize=(8,15))
g = sns.heatmap(data=df_heatmap, ax=ax,cbar_kws={'label': 'coverage difference of GFP and GFP-SNED1','shrink': 0.5},cmap="RdYlGn",yticklabels=True)
ax.set_xticklabels(labels=columns, rotation=30, fontsize=10, ha='right')
# ax.axes.yaxis.set_visible(False)
# g = sns.clustermap(data=df_heatmap,cbar_kws={'label': 'coverage','ticks': range(0,120,20),'shrink': 0.3},vmin=0,vmax=100, cmap="YlGnBu",yticklabels=False)
# g.ax_heatmap.set_ylabel('ECM protein')
# plt.setp(g.ax_heatmap.get_xticklabels(), rotation=30, ha='right')
plt.setp(ax.get_yticklabels(), fontsize=5)
# plt.savefig('D:/data/Naba_deep_matrisome/07232021_secondsearch/test.png', dpi=300)
# plt.show()
"""

### boxplot and dots connecting

"""
ax.legend(lines, labels, framealpha=0.5,loc='upper right',fontsize=15)

x_normal = np.random.normal(1, 0.04, size=len(aggre_ko18_cov))
x_aggre = np.random.normal(2, 0.04, size=len(aggre_18_cov))
color_list = [color_map[prot] for prot in df_ecm_aggre.index]
ax.boxplot([normal18GFP_cov,normal18SNED_cov], labels=['GFP 18h', 'GFP-SNED1 18h'],showfliers=False)
# plot points onto boxplot
for x,cov_list in zip([x_normal,x_aggre],[normal18GFP_cov,normal18SNED_cov]):
    ax.plot(x,cov_list,'ko',markersize=10,alpha=0.5)

#connect dots with different color
for x1,y1,x2,y2,color in zip(x_normal,normal18GFP_cov,x_aggre,normal18SNED_cov,color_list):
    ax.plot([x1,x2],[y1,y2],color=color, linestyle='--', alpha=0.8)

# for x1,y1,x2,y2 in zip(sub_x_163,sub_y_163,sub_x_182,sub_y_182):
#     ax.plot([x1, x2], [y1, y2], color=ecm_class_color_dict['ECM Glycoproteins'], linestyle='--', alpha=1)
# ax.set_title('Coverage comparison between 18h and aggregated 18h in GFP', fontsize=22)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_ylabel('%coverage',fontsize=20)
plt.show()
"""

# ptm analysis

# from tsv_reader import psm_ptm_reader
# base_path = 'D:/data/Naba_deep_matrisome/05142021_secondsearch/'
# folders = [base_path+each+folder for each in ['KOC/'] for folder in os.listdir(base_path+each) if '.' not in folder]
# print (folders)
# psm_path_list = [each+'/psm.tsv' for each in folders if 'Re4' not in each]
#
# ptm_psm_dict,total_psm = psm_ptm_reader(psm_path_list,gene_set=ecm_gene_category_dict,mod=0.9840)
# print ({each:len(ptm_psm_dict[each]) for each in ptm_psm_dict},total_psm)

### compare paraelle and sequential coverage
"""
# sequential_cov = 'D:/data/Naba_deep_matrisome/02152021_1/18_2_accumulated_cov.csv'
# seq_df = pd.read_csv(sequential_cov, index_col=0)
# 
# 
# paraelle_cov = 'D:/data/Naba_deep_matrisome/05142021_secondsearch/6_14_summary_B_C_Xinhao.xlsx'
# paraelle_df = pd.read_excel(paraelle_cov,index_col=0)
# paraelle_4h = []
# total_prot_set = set(seq_df.columns.tolist()+[prot for prot, gene in zip(paraelle_df.index.tolist(),paraelle_df['gene']) if gene in ecm_gene_category_dict])
# print (len(total_prot_set))
# 
# sequential_4h_cov = []
# paraelle_4h_cov = []
# 
# for each in total_prot_set:
#     if each in seq_df.columns and each in paraelle_df.index:
#         sequential_4h_cov.append(seq_df.at['18_2B4',each])
#         paraelle_4h_cov.append(np.mean([paraelle_df.at[each,'KO_B_4_coverage'],paraelle_df.at[each,'KO_C_4_coverage']],dtype=np.float64))
#     elif each in seq_df.columns and each not in paraelle_df.index:
#         sequential_4h_cov.append(seq_df.at['18_2B4',each])
#         paraelle_4h_cov.append(0)
#     elif each not in seq_df.columns and each in paraelle_df.index:
#         sequential_4h_cov.append(0)
#         paraelle_4h_cov.append(np.mean([paraelle_df.at[each,'KO_B_4_coverage'],paraelle_df.at[each,'KO_C_4_coverage']],dtype=np.float64))
# 
# plt.boxplot([sequential_4h_cov,paraelle_4h_cov], labels=['sequential_4h_cov', 'paraelle_4h_cov'],showfliers=False)
# plt.show()
"""

### coverage line_plot
"""
from scipy.interpolate import interp1d
base_path = 'D:/data/Naba_deep_matrisome/05142021_secondsearch/'
folders = [base_path+each+folder for each in ['KOB/','KOC/', 'SNEDB/', 'SNEDC/'] for folder in os.listdir(base_path+each) if '.' not in folder]
print (folders)
psm_path_list = [each+'/psm.tsv' for each in folders]
pep_path_list = [each+'/peptide.tsv' for each in folders]

# file_protein_peptide_spec_dict = dash_dataframe(pep_path_list,psm_path_list,protein_dict,psm_path_list,psm_path_list,psm_path_list)
file_protein_pep_spec_dict = ppp.load(open('B_C_file_protein_pep_spec_dict.p','rb'))

prot = 'Q60716'
protein_seq = protein_dict[prot]
ko_file = 'KO_C_0o5'
sned_file = 'SNED_C_0o5'
ko_protein_pep_dict,sned_protein_pep_dict = file_protein_pep_spec_dict[ko_file][prot],file_protein_pep_spec_dict[sned_file][prot]


aa_array_ko,aa_array_sned = np.zeros(len(protein_seq)), np.zeros(len(protein_seq))
for pep in ko_protein_pep_dict:
    pep_start_index = protein_seq.find(pep)
    pep_end_index = pep_start_index+len(pep)
    aa_array_ko[pep_start_index:pep_end_index]+=ko_protein_pep_dict[pep]

for pep in sned_protein_pep_dict:
    pep_start_index = protein_seq.find(pep)
    pep_end_index = pep_start_index+len(pep)
    aa_array_sned[pep_start_index:pep_end_index]+=sned_protein_pep_dict[pep]

x = range(len(protein_seq))
y = aa_array_ko
# f_cubic = interp1d(x,y,kind='cubic')
# x_new = np.linspace(0,len(protein_seq)-1,10000)
plt.plot(y) 
plt.show()

"""
### swarmplot on top of boxplot
"""

# df_new = pd.DataFrame(dict(gene=df_aggre_coverage.gene,category=df_aggre_coverage.category,GFP18_normal=normal18GFP_cov,
#                            SNED18_normal=normal18SNED_cov,diff=np.array(normal18GFP_cov)-np.array(normal18SNED_cov)),index=df_aggre_coverage.index)
# df_new.to_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/matrisome_18_standard_diff.xlsx')
# ecm_protein_list = df_aggre_coverage.index.tolist()
df_ecm_aggre['euclidean_dis/cosine_sim'] = df_ecm_aggre['euclidean_distance with +-']/df_ecm_aggre['cosine_sim']


fig, ax = plt.subplots(1,1,figsize=(6,6))
ax = sns.boxplot(x='category', y='euclidean_dis/cosine_sim', hue='category',palette=ecm_class_color_dict, data=df_ecm_aggre,linewidth=2.5,
                 order=["ECM Glycoproteins", "Collagens", "Proteoglycans", "ECM-affiliated Proteins",
                                "ECM Regulators","Secreted Factors"],
                 hue_order=["ECM Glycoproteins", "Collagens", "Proteoglycans", "ECM-affiliated Proteins",
                                "ECM Regulators","Secreted Factors"], dodge=False)
ax = sns.swarmplot(x='category', y='euclidean_dis/cosine_sim', data=df_ecm_aggre,color=".2",size=5,
                   order=["ECM Glycoproteins", "Collagens", "Proteoglycans", "ECM-affiliated Proteins",
                                "ECM Regulators","Secreted Factors"])
ax.get_legend().remove()
ax.set_ylabel('Dissimilarity index', fontsize=8)
# ax.set_ylabel('Coverage difference')
plt.xticks(rotation=30,fontsize=8)
plt.tight_layout()
plt.savefig('D:/data/Naba_deep_matrisome/07232021_secondsearch/figure_update/test.png', dpi=300)
plt.show()
"""

### nasf from combined protein tsv file
combined_df = pd.read_csv('F:/sned1_biotinalytion/07072022/search/combined_protein.tsv',delimiter='\t',index_col=0)
nsaf_df = pd.DataFrame(columns=['Protein ID','Gene', 'Protein Length','Description','ColI_FullReaction nsaf',
                                'ColI_NoBiotin nsaf','ColI_NoPrimary nsaf','FN_FullReaction nsaf',
                                'FN_NoBiotin nsaf','FN_NoPrimary nsaf'])

protein_list, protein_length_list = combined_df['Protein ID'].tolist(), combined_df['Protein Length'].tolist()
columns = [each for each in combined_df.columns if 'Total Spectral Count' in each]
for clm in columns:
        if clm.split(' ')[0]+' nsaf' in nsaf_df.columns:
                print(clm.split(' ')[0])
                spec_count_list = combined_df[clm].tolist()
                total_nsaf = sum([val/protein_length_list[idx] for idx,val in enumerate(spec_count_list)])
                nsaf_df[clm.split(' ')[0]+' nsaf'] = [val/protein_length_list[idx]/total_nsaf
                                                      for idx,val in enumerate(spec_count_list)]

nsaf_df['Protein ID'] = combined_df['Protein ID'].tolist()
nsaf_df['Gene'] = combined_df['Gene'].tolist()
nsaf_df['Protein Length'] = combined_df['Protein Length'].tolist()
nsaf_df['Description'] = combined_df['Description'].tolist()

nsaf_df.to_excel('F:/sned1_biotinalytion/07072022/search/nsaf.xlsx')