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
"""
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

df_summary = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_summary_D_F_squential_standard.xlsx',index_col=0)
"""
# df_aggre_coverage = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/8_1_matrisome_average_aggre.xlsx', index_col=0)
# normal18GFP_cov = [np.mean([df_summary.at[prot,'GFP_1080D_coverage'], df_summary.at[prot,'GFP_1080F_coverage']]) for prot in df_aggre_coverage.index]
# normal2GFP_cov = [np.mean([df_summary.at[prot,'GFP_120D_coverage'], df_summary.at[prot,'GFP_120F_coverage']]) for prot in df_aggre_coverage.index]
# normal18SNED_cov = [np.mean([df_summary.at[prot,'SNED1_1080D_coverage'], df_summary.at[prot,'SNED1_1080F_coverage']]) for prot in df_aggre_coverage.index]
# normal2SNED_cov = [np.mean([df_summary.at[prot,'SNED1_120D_coverage'], df_summary.at[prot,'SNED1_120F_coverage']]) for prot in df_aggre_coverage.index]



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
### get matrisome non-aggregated coverage data
"""
columns = ['gene','MW','category','GFP_30_ave','GFP_120_ave','GFP_240_ave','GFP_1080_ave',
           'SNED1_30_ave','SNED1_120_ave','SNED1_240_ave','SNED1_1080_ave']
df_ecm_original = pd.DataFrame(columns=columns,index=df_aggre_coverage.index)
for each_prot in df_aggre_coverage.index:
    df_ecm_original.at[each_prot,'gene'] = df_aggre_coverage.at[each_prot,'gene']
    df_ecm_original.at[each_prot,'MW'] = df_aggre_coverage.at[each_prot,'MW_kDa']
    df_ecm_original.at[each_prot,'category'] = df_aggre_coverage.at[each_prot,'category']
    for sample in ['GFP','SNED1']:
        for time in ['_30','_120','_240','_1080']:
            df_ecm_original.at[each_prot,sample+time+'_ave'] = np.mean([df_summary.at[each_prot,sample+'_seq'+time+'D_coverage'], df_summary.at[each_prot,sample+'_seq'+time+'F_coverage']])


df_ecm_original.to_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/ecm_original_cov_ave.xlsx')
"""


### NSAF analysis

from collections import defaultdict
fasta_path = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
protein_dict = fasta_reader(fasta_path)

sp_tr_dict = defaultdict(set)
with open(fasta_path,'r') as f:
    f_split = f.read().split('\n>')
    for each in f_split[1:]:
        sp_tr_dict[each.split('|')[0]].add(each.split('|')[1])
    sp_tr_dict[f_split[0].split('|')[0].split('>')[1]].add(f_split[0].split('|')[1])

ecm_df = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/8_1_matrisome_average_aggre.xlsx',index_col=0)
ecm_list = ecm_df.index.tolist()
ecm_list_sp = [each for each in ecm_list if each in sp_tr_dict['sp']]

df_DF = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_summary_D_F_squential_standard.xlsx',index_col=0)
selected_columns = [each for each in df_DF.columns if '_total_spec' in each]
new_columns = [each.replace('_total_spec','nsaf') for each in selected_columns]
# get the sum of nsaf for each condition
column_total_nsaf_dict = {each:sum([df_DF.at[prot,each]/len(prot) for prot in ecm_list_sp]) for each in selected_columns}
df = pd.DataFrame(index=ecm_list_sp, columns=['gene','ECM category']+new_columns)

for protein in ecm_list_sp:
    df.at[protein,'gene'] = ecm_df.at[protein,'gene']
    df.at[protein,'ECM category'] = ecm_df.at[protein,'category']
    for clm,old_clm in zip(new_columns,selected_columns):
        df.at[protein,clm] = df_DF.at[protein,old_clm]/len(protein_dict[protein])/column_total_nsaf_dict[old_clm]

df.to_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/nsaf_raw.xlsx')