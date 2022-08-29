import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
# excel_peptide_excel = 'F:/native_digestion/Uchicago_TMT/tmt_search_0824/tmt-report/abundance_peptide_GN.tsv'
# df = pd.read_csv(excel_peptide_excel,delimiter='\t',index_col=0)
# df = df.fillna(0)
# columns = df.columns
# data_columns = [col for col in df.columns if '_' in col]
# print (data_columns)
#
# new_df = pd.DataFrame(index=df.index,columns=df.columns)
# # convert log2 int back to original value and normalize as percentage to reference
# for each in df.itertuples():
#     index = each[0]
#     new_df.loc[index,columns[:4]]=each[1:5]
#     reference = pow(2,each[5])
#     new_df.loc[index,columns[4]]=reference
#     for ind, val in enumerate(data_columns):
#
#         new_df.loc[index,val] = pow(2,each[6+ind])/reference
# new_df.to_csv('F:/native_digestion/Uchicago_TMT/tmt_search_0824/tmt-report/abundance_peptide_GN_normalized_0824.tsv',sep='\t')


### combine peptides intensities (divided by pool) from 3 tmt peptide.tsv together
"""
df1 = pd.read_csv('F:/native_digestion/Uchicago_TMT/tmt_search_0826/TMT1/peptide.tsv',delimiter='\t',index_col=0)
df2 = pd.read_csv('F:/native_digestion/Uchicago_TMT/tmt_search_0826/TMT2/peptide.tsv',delimiter='\t',index_col=0)
df3 = pd.read_csv('F:/native_digestion/Uchicago_TMT/tmt_search_0826/TMT3/peptide.tsv',delimiter='\t',index_col=0)
df_list = [df1,df2,df3]

pep_list = df1.index.tolist()+df2.index.tolist()+df3.index.tolist()
pep_set = set(pep_list)
df = pd.DataFrame(index=pep_set)

for pep in pep_set:
    # print (pep)
    for data in df_list:
        peptide_ids = data.index.tolist()
        last_10_column = data.columns[-10:]
        bridge_column = last_10_column[0]
        if pep in peptide_ids:
            pool_int = data.at[pep, bridge_column]
            if pool_int!=0:
                for time in last_10_column[1:]:
                    df.at[pep,time] = data.at[pep,time]/pool_int
            else:
                df.at[pep,last_10_column[1:]] = [0]*9
        else:
            df.at[pep,last_10_column[1:]] = [0]*9
df.to_csv('F:/native_digestion/Uchicago_TMT/tmt_search_0826/pep_normalized_int.tsv',sep='\t')
"""
### normalize tmt ratio between 0-1
"""
df = pd.read_csv('F:/native_digestion/Uchicago_TMT/tmt_search_0826/pep_normalized_int.tsv',delimiter='\t',index_col=0)

array = df.to_numpy()
max_int = np.max(array,axis=1) # 1d array
max_int_array = np.reshape(np.repeat(max_int,df.shape[1]),(df.shape[0],df.shape[1]))

new_array = np.divide(array,max_int_array)
new_df = pd.DataFrame(new_array, index=df.index,columns=df.columns).fillna(0)
# new_df.to_csv('F:/native_digestion/Uchicago_TMT/tmt_search_0826/peptide_tmt_normalized_0_1.tsv',sep='\t')
"""
### plot
import pymannkendall as mk
df = pd.read_csv('F:/native_digestion/Uchicago_TMT/tmt_search_0826/distance_tmt_weighted_0826.tsv',delimiter='\t',index_col=0)
new_df = df.loc[:,[each for each in list(df.columns) if 'chymo' in each][:-1]]
mk_result = mk.original_test(new_df.median().tolist())
print(mk_result)
#
plot_array = new_df.to_numpy().flatten()
x = new_df.columns.tolist()*new_df.shape[0]
df_plot = pd.DataFrame(dict(sample=x,tmt_weight_distance=plot_array))

fig,axs = plt.subplots(figsize=(10,6))

sns.boxplot(data=df_plot,ax=axs, x='sample',y='tmt_weight_distance',linewidth=2.5)
axs.set_xticklabels(list(new_df.columns), fontsize=10,ha="center", rotation=45)
plt.show()


