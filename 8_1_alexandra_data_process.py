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

matrisomedb_excel = 'D:/data/Naba_deep_matrisome/matrisome coverage_norepeat.xlsx'
df_matrisome = pd.read_excel(matrisomedb_excel)
ecm_category_dict = {i:j for i,j in zip(df_matrisome['protein_id'].tolist(),df_matrisome['category'].tolist())}
print (len(ecm_category_dict))

df_aggre_coverage = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/8_1_matrisome_average_aggre.xlsx',index_col=0)
ecm_category_dict_2  = {i:j for i,j in zip(df_aggre_coverage.index.tolist(),df_aggre_coverage['category'].tolist())}
print (len(ecm_category_dict_2))
ecm_category_dict.update(ecm_category_dict_2)

# df_summary = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_summary_D_F_squential_standard.xlsx',index_col=0)
# normal18GFP_cov = [np.mean([df_summary.at[prot,'GFP_1080D_coverage'], df_summary.at[prot,'GFP_1080F_coverage']]) for prot in df_aggre_coverage.index]
# normal2GFP_cov = [np.mean([df_summary.at[prot,'GFP_120D_coverage'], df_summary.at[prot,'GFP_120F_coverage']]) for prot in df_aggre_coverage.index]
# normal18SNED_cov = [np.mean([df_summary.at[prot,'SNED1_1080D_coverage'], df_summary.at[prot,'SNED1_1080F_coverage']]) for prot in df_aggre_coverage.index]
# normal2SNED_cov = [np.mean([df_summary.at[prot,'SNED1_120D_coverage'], df_summary.at[prot,'SNED1_120F_coverage']]) for prot in df_aggre_coverage.index]


# df_aggre_coverage = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/8_1_matrisome_average_aggre.xlsx', index_col=0)
#
# for each_tp in df_aggre_coverage.itertuples():
#     print (each_tp[4:8],each_tp[-4:], cosine_sim_calculating(each_tp[4:8],each_tp[-4:]), np.linalg.norm(np.array(each_tp[4:8])-np.array(each_tp[-4:])))
#     df_aggre_coverage.at[each_tp[0],'cosine_sim'] = cosine_sim_calculating(each_tp[4:8],each_tp[-4:])
#     df_aggre_coverage.at[each_tp[0],'euclidean_distance with +-'] = np.linalg.norm(np.array(each_tp[4:8])-np.array(each_tp[8:12])) \
#         if each_tp[7]-each_tp[11]>=0 else -np.linalg.norm(np.array(each_tp[4:8])-np.array(each_tp[8:12]))
# df_aggre_coverage.to_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/8_1_matrisome_average_aggre.xlsx')


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
"""
bcdf_combined_agg = pd.read_excel('D:/data/Naba_deep_matrisome/BCDF_combined/BCDF_ecm_aggregated.xlsx',index_col=0)
samples, times = ['GFP_seq_', 'SNED1_seq_'],['30','120','240','1080']

df_aggre_coverage = pd.DataFrame()
for prot in bcdf_combined_agg.index:
    gene = bcdf_combined_agg.at[prot,'gene']
    category = bcdf_combined_agg.at[prot,'ECM_category']
    df_aggre_coverage.at[prot,'gene'] = gene
    df_aggre_coverage.at[prot, 'ECM category'] = category
    for sample in samples:
        for time in times:
            # if time != '30':
            #     aver_cov = np.mean([bcdf_combined_agg.at[prot,sample+'D_'+time+'_aggre_coverage'],bcdf_combined_agg.at[prot,sample+'F_'+time+'_aggre_coverage']])
            #     df_aggre_coverage.at[prot,sample+time+'_ave_aggre_cov']=aver_cov
            # else:
            #     # add 30min data BC into average
            #     aver_cov = np.mean([bcdf_combined_agg.at[prot,sample+'D_'+time+'_aggre_coverage'],bcdf_combined_agg.at[prot,sample+'F_'+time+'_aggre_coverage'],
            #                         bcdf_combined_agg.at[prot, sample.split('_')[0]+'_parallel_B_30_aggre_coverage'], bcdf_combined_agg.at[prot, sample.split('_')[0]+'_parallel_C_30_aggre_coverage']])
            #     df_aggre_coverage.at[prot,sample+time+'_ave_aggre_cov']=aver_cov

            aver_cov = np.mean([bcdf_combined_agg.at[prot, sample + 'D_' + time + '_aggre_coverage'],
                                bcdf_combined_agg.at[prot, sample + 'F_' + time + '_aggre_coverage'],
                                bcdf_combined_agg.at[prot, sample.split('_')[0] + '_parallel_B_'+time+'_aggre_coverage'],
                                bcdf_combined_agg.at[prot, sample.split('_')[0] + '_parallel_C_'+time+'_aggre_coverage']])
            df_aggre_coverage.at[prot, sample + time + '_ave_aggre_cov'] = aver_cov
"""
df_aggre_coverage=pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/8_1_matrisome_average_aggre.xlsx',index_col=0)
columns = list(df_aggre_coverage)[3:7]

df_group = df_aggre_coverage.groupby(['category'])[columns].mean().transpose()
print (df_group)

ecm_class_color_dict = {"Collagens": '#3DA8FB', 'ECM-affiliated Proteins':'#23AECA',
                        'ECM Regulators':"#23CA66","Secreted Factors":"#E3EB09",
                        "ECM Glycoproteins":"#EBA709", "Proteoglycans":"#EB09DC"}
fig, ax = plt.subplots()

# Initialize the bottom at zero for the first set of bars.
bottom = np.zeros(len(df_group))

# Plot each layer of the bar, adding each bar to the "bottom" so
# the next bar starts higher.
for i, col in enumerate(df_group.columns):
  ax.bar(df_group.index, df_group[col], bottom=bottom, label=col, )
  bottom += np.array(df_group[col])
ax.set_xticklabels(['30min', '2h', '4h', '18h'], fontsize=15)
ax.set_title('Tips by Day and Gender')
ax.legend()
plt.show()