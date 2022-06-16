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
ecm_class_color_dict = {"Collagens": '#0584B7', 'ECM-affiliated Proteins':'#F4511E',
                        'ECM Regulators':"#F9A287","Secreted Factors":"#FFE188",
                        "ECM Glycoproteins":"#133463", "Proteoglycans":"#59D8E6"}

### NSAF analysis
"""
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

df_DF = pd.read_excel('D:/data/Naba_deep_matrisome/BCDF_combined/BCDF_combined_summary_sp_only.xlsx',index_col=0)
total_protein_list = df_DF.index
# selected_columns = [each for each in df_DF.columns if '_total_spec' in each]
# new_columns = [each.replace('_total_spec','nsaf') for each in selected_columns]
# print (selected_columns)
# # get the sum of nsaf for each condition
#
# column_total_nsaf_dict = {each:sum([df_DF.at[prot,each]/len(protein_dict[prot]) for prot in total_protein_list]) for each in selected_columns}
# df = pd.DataFrame(index=ecm_list_sp, columns=['gene','ECM category']+new_columns)
#
# for protein in ecm_list_sp:
#     df.at[protein,'gene'] = ecm_df.at[protein,'gene']
#     df.at[protein,'ECM category'] = ecm_df.at[protein,'category']
#     for clm,old_clm in zip(new_columns,selected_columns):
#         df.at[protein,clm] = (df_DF.at[protein,old_clm]/len(protein_dict[protein]))/column_total_nsaf_dict[old_clm]
#
# df.to_excel('D:/data/Naba_deep_matrisome/BCDF_combined/NSAF_bg_all_proteins.xlsx')
"""

### aggregated NSAF analysis
"""
selected_columns = [each for each in df_DF.columns if 'seq' in each and 'total_spec' in each]
time_points = ['30','120','240','1080']
print (selected_columns)
column_total_nsaf_dict ={}
for each in ['GFP','SNED1']:
    total_nsaf = 0
    for prot in total_protein_list:
        total_spec = sum([df_DF.at[prot,each+'_seq_'+j+i+'_total_spec'] for i in ['D','F'] for j in time_points])
        prot_nsaf = total_spec/len(protein_dict[prot])
        total_nsaf+=prot_nsaf
    column_total_nsaf_dict[each]=total_nsaf

nsaf_df = pd.DataFrame(index=ecm_list_sp,columns=['gene','MW_kDa','category','GFP_NSAF','SNED1_NSAF', 'GFP_SNED1_diff'])
for prot in ecm_list_sp:
    nsaf_df.at[prot,'gene'] = ecm_df.at[prot,'gene']
    nsaf_df.at[prot,'MW_kDa'] = ecm_df.at[prot,'MW_kDa']
    nsaf_df.at[prot,'category'] = ecm_df.at[prot,'category']
    gfp_nsaf = sum([df_DF.at[prot,'GFP_seq_'+j+i+'_total_spec'] for i in ['D','F'] for j in time_points])/len(protein_dict[prot])/column_total_nsaf_dict['GFP']
    sned_nsaf = sum([df_DF.at[prot,'SNED1_seq_'+j+i+'_total_spec'] for i in ['D','F'] for j in time_points])/len(protein_dict[prot])/column_total_nsaf_dict['SNED1']
    nsaf_df.at[prot,'GFP_NSAF'] = gfp_nsaf
    nsaf_df.at[prot,'SNED1_NSAF'] = sned_nsaf
    nsaf_df.at[prot,'GFP_SNED1_diff'] = gfp_nsaf-sned_nsaf

nsaf_df.to_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/DF_average_NSAF.xlsx')
"""

### plot NSAF difference versus MW
"""
df_average_nsaf = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/DF_average_NSAF.xlsx',index_col=0)
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
q95,q5 = np.percentile(df_average_nsaf['GFP_SNED1_diff'],[95,5])
print (q95,q5)
df_sub = df_average_nsaf[(df_average_nsaf['GFP_SNED1_diff']>q95)|(df_average_nsaf['GFP_SNED1_diff']<q5)]
print (df_sub)
fig,ax = plt.subplots(1,1,figsize=(10,6))
g = sns.scatterplot(data=df_average_nsaf,ax=ax,x='MW_kDa',y='GFP_SNED1_diff',hue='category',palette=ecm_class_color_dict,alpha=0.9,s=45)
ax.get_legend().remove()
for each in df_sub.index:
    gene = df_sub.at[each,'gene']
    category = df_sub.at[each,'category']
    x = df_sub.at[each,'MW_kDa']
    y = df_sub.at[each,'GFP_SNED1_diff']
    plt.text(x=x+0.3,y=y+0.00001,s=gene,fontdict=dict(color=ecm_class_color_dict[category],size=12))

plt.savefig('D:/data/Naba_deep_matrisome/07232021_secondsearch/figure_update/GFP_SNED1_aggre_NSAF_diff_MW.png', dpi=300)
plt.show()
"""
### peptide index mapping from tsv
"""
from multiprocessing_naive_algorithym import map_peptide_index
from glob import glob
import os
fasta_path = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
protein_dict = fasta_reader(fasta_path)

# base_path = 'D:/data/Naba_deep_matrisome/05142021_secondsearch/'
# folders = [base_path+each+folder for each in ['KOB/','KOC/', 'SNEDB/', 'SNEDC/'] for folder in os.listdir(base_path+each) if '.' not in folder]
# files = [each+'/peptide.tsv' for each in folders if 'Re4' not in each]


base_path = 'D:/data/Naba_deep_matrisome/07232021_secondsearch/'
files = [base_path+folder+'/peptide.tsv' for folder in os.listdir(base_path) if 'GFP' in folder or 'SNED1' in folder]


with pd.ExcelWriter('D:/data/Naba_deep_matrisome/PSM_info_DF_XS.xlsx') as writer:
    for each_file in files:
        with open(each_file,'r') as f_o:
            next(f_o)
            peptide_gene_spec_list = [(line.split('\t')[0],line.split('\t')[11], line.split('\t')[4]) for line in f_o]
            pep_list, gene_list, spec_count_list = zip(*peptide_gene_spec_list)
            pep_start_end_dict = map_peptide_index(protein_dict,pep_list)
            df = pd.DataFrame(dict(gene=gene_list,
                                   protein_start_end=[pep_start_end_dict[each] for each in pep_list],
                                   spec_count = spec_count_list),index=pep_list)
            df.to_excel(writer,sheet_name=each_file.split('/')[-2])
            print(each_file)
"""
### get total spec count from summary excel
excel_path = 'D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_summary_D_F.xlsx'
df = pd.read_excel(excel_path,index_col=0)

ecm_excel_path = 'D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_ecm_aggregated_D_F.xlsx'
df_ecm = pd.read_excel(ecm_excel_path, index_col=0)

df_mask = df.loc[df_ecm.index,:]
# print (df_mask)
columns = [each for each in df_mask.columns if 'total_spec' in each and '_seq_' in each]
for col in columns:
    print (col, df_mask[col].sum())