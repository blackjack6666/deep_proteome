import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from protein_coverage import fasta_reader

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


# df_aggre_coverage = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/8_1_matrisome_average_aggre.xlsx', index_col=0)
#
# for each_tp in df_aggre_coverage.itertuples():
#     print (each_tp[4:8],each_tp[-4:], cosine_sim_calculating(each_tp[4:8],each_tp[-4:]), np.linalg.norm(np.array(each_tp[4:8])-np.array(each_tp[-4:])))
#     df_aggre_coverage.at[each_tp[0],'cosine_sim'] = cosine_sim_calculating(each_tp[4:8],each_tp[-4:])
#     df_aggre_coverage.at[each_tp[0],'euclidean_distance'] = np.linalg.norm(np.array(each_tp[4:8])-np.array(each_tp[-4:]))

df_aggre_coverage = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/8_1_matrisome_average_aggre.xlsx',index_col=0)
ecm_protein_list = df_aggre_coverage.index.tolist()

# ax = sns.boxplot(x='category', y='euclidean_distance', data=df_aggre_coverage)
# plt.xticks(rotation=30)
# plt.show()

summary_df = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_summary_D_F_squential_standard.xlsx', index_col=0)
summary_df_slice = summary_df[summary_df.index.isin(ecm_protein_list)]

print (summary_df_slice['GFP_1080D_total peptides identified'].unique())
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

print (summary_df_slice[(summary_df_slice['GFP_1080D_total_spec']!=0)|(summary_df_slice['GFP_1080F_total_spec']!=0)].shape)