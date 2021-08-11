import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from protein_coverage import fasta_reader
from scipy.stats import ttest_1samp

# xls = pd.ExcelFile('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_summary_aggregated_D_F_Received_210727_AN.xlsx')
# df_matrisome_AN = pd.read_excel(xls,'Matrisome_Coverage',index_col=0)

def cosine_sim_calculating(v1, v2):
    """
    calculate the cosine similarity beweeen two b/y ions binned vectors
    :param v1:
    :param v2:
    :return:
    """
    from scipy import spatial
    return 1-spatial.distance.cosine(v1,v2)

fasta_path = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
protein_dict = fasta_reader(fasta_path)

# df_aggre_coverage = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/8_1_matrisome_average_aggre.xlsx',index_col=0)
# df_summary = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_summary_D_F_squential_standard.xlsx',index_col=0)
# normal18GFP_cov = [np.mean([df_summary.at[prot,'GFP_1080D_coverage'], df_summary.at[prot,'GFP_1080F_coverage']]) for prot in df_aggre_coverage.index]
# normal2GFP_cov = [np.mean([df_summary.at[prot,'GFP_120D_coverage'], df_summary.at[prot,'GFP_120F_coverage']]) for prot in df_aggre_coverage.index]
# normal18SNED_cov = [np.mean([df_summary.at[prot,'SNED1_1080D_coverage'], df_summary.at[prot,'SNED1_1080F_coverage']]) for prot in df_aggre_coverage.index]
# normal2SNED_cov = [np.mean([df_summary.at[prot,'SNED1_120D_coverage'], df_summary.at[prot,'SNED1_120F_coverage']]) for prot in df_aggre_coverage.index]

### combine BC DF
df_combined = pd.DataFrame()
df_BC = pd.read_excel('D:/data/Naba_deep_matrisome/05142021_secondsearch/6_14_summary_B_C_Xinhao.xlsx',index_col=0)
df_DF = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_summary_D_F_squential_standard.xlsx',index_col=0)
total_protein_list = set(df_BC.index.tolist()+df_DF.index.tolist())
for each in total_protein_list:
    df_combined.at[each,'gene'] = df_BC.at[each,'gene']
    df_combined.at[each,'MW'] = df_BC.at[each,'MW']


# df_aggre_coverage = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/8_1_matrisome_average_aggre.xlsx', index_col=0)
#
# for each_tp in df_aggre_coverage.itertuples():
#     print (each_tp[4:8],each_tp[-4:], cosine_sim_calculating(each_tp[4:8],each_tp[-4:]), np.linalg.norm(np.array(each_tp[4:8])-np.array(each_tp[-4:])))
#     df_aggre_coverage.at[each_tp[0],'cosine_sim'] = cosine_sim_calculating(each_tp[4:8],each_tp[-4:])
#     df_aggre_coverage.at[each_tp[0],'euclidean_distance with +-'] = np.linalg.norm(np.array(each_tp[4:8])-np.array(each_tp[8:12])) \
#         if each_tp[7]-each_tp[11]>=0 else -np.linalg.norm(np.array(each_tp[4:8])-np.array(each_tp[8:12]))
# df_aggre_coverage.to_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/8_1_matrisome_average_aggre.xlsx')

### swarmplot on top of boxplot
"""
df_new = pd.DataFrame(dict(gene=df_aggre_coverage.gene,category=df_aggre_coverage.category,GFP18_normal=normal18GFP_cov,
                           SNED18_normal=normal18SNED_cov,diff=np.array(normal18GFP_cov)-np.array(normal18SNED_cov)),index=df_aggre_coverage.index)
df_new.to_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/matrisome_18_standard_diff.xlsx')
ecm_protein_list = df_aggre_coverage.index.tolist()
# df_aggre_coverage['euclidean_dis/cosine_sim'] = df_aggre_coverage['euclidean_distance with +-']/df_aggre_coverage['cosine_sim']
fig, ax = plt.subplots(1,1,figsize=(6,9))
ax = sns.boxplot(x='category', y='diff', data=df_new)
ax = sns.swarmplot(x='category', y='diff', data=df_new,color=".2",size=4)
# ax.set_ylabel('Dissimilarity index')
ax.set_ylabel('Coverage difference')
plt.xticks(rotation=30)
plt.savefig('D:/data/Naba_deep_matrisome/07232021_secondsearch/test.png', dpi=300)
plt.show()
"""
# for each in df_aggre_coverage['category'].unique():
#     sub_df = df_aggre_coverage[df_aggre_coverage['category']==each]
#     print (each, ttest_1samp(sub_df['euclidean_distance with +-'],0.0,alternative='greater'))

# summary_df = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_summary_D_F_squential_standard.xlsx', index_col=0)
# summary_df_slice = summary_df[summary_df.index.isin(ecm_protein_list)]

# print (summary_df_slice['GFP_1080D_total peptides identified'].unique())
# peptide_set = set()
# for column in ['SNED1_1080D_total peptides identified', 'SNED1_1080F_total peptides identified']:
#     for each in summary_df_slice[column].unique():
#         print (each)
#         if each != 0:
#             if ',' in each:
#                 for pep in each.split(', '):
#                     peptide_set.add(pep)
#             else:
#                 peptide_set.add(each)
#         else:
#             continue
# print (len(peptide_set))
