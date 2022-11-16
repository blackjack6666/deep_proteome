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

### combine peptides intensities (divided by pool) from 3 tmt peptide.tsv together, * to do: pool value needs to be normalized?
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
### combine trypsin/chymotrypsin data from pep_normalized_int.tsv
"""
df = pd.read_csv('F:/native_digestion/Uchicago_TMT/tmt_search_0826/pep_normalized_int.tsv',delimiter='\t',index_col=0)
array = df.to_numpy()
array_trypsin, array_chymo = array[:,:12], array[:,14:26]
combined_array = np.add(array_trypsin,array_chymo)
new_df = pd.DataFrame(combined_array,index=df.index, columns=[str(i)+'h_combined' for i in range(1,13)])
new_df.to_csv('F:/native_digestion/Uchicago_TMT/tmt_search_0826/pep_normalized_int_combined.tsv',sep='\t')
"""
### normalize tmt ratio between 0-1, didn't work
"""
df = pd.read_csv('F:/native_digestion/Uchicago_TMT/tmt_search_0826/pep_normalized_int.tsv',delimiter='\t',index_col=0)
df_tryp, df_chy = df.iloc[:,:14], df.iloc[:,14:]

## normalized bewteen 0 and 1
array = df.to_numpy()
# max_int = np.max(array,axis=1) # 1d array
# max_int_array = np.reshape(np.repeat(max_int,df.shape[1]),(df.shape[0],df.shape[1]))
# new_array = np.divide(array,max_int_array)
# new_df = pd.DataFrame(new_array, index=df.index,columns=df.columns).fillna(0)
# new_df.to_csv('F:/native_digestion/Uchicago_TMT/tmt_search_0826/peptide_tmt_normalized_combined_0_1.tsv',sep='\t')

## normalize each column
max_int = np.max(array,axis=0)
print (max_int)
max_int_array = np.tile(max_int,(df.shape[0],1))
new_array = np.divide(array,max_int_array)
# normalize each row (peptide) after column (channel)
max_int_new = np.max(new_array,axis=1)
max_int_array_new = np.reshape(np.repeat(max_int_new,df.shape[1]),(df.shape[0],df.shape[1]))
array_final = np.divide(new_array,max_int_array_new)
new_df = pd.DataFrame(array_final, index=df.index,columns=df.columns)

new_df.to_csv('F:/native_digestion/Uchicago_TMT/tmt_search_0826/peptide_tmt_normalized_bycolumnandrow_01.tsv',sep='\t')

## normalized between -1 and 1, negative value is penalty

# for each_df,name in zip([df_tryp,df_chy],['tryp','chymo']):
#     new_df = pd.DataFrame(index=each_df.index,columns=each_df.columns)
#     array2d = each_df.to_numpy()
#     for each, peptide in zip(array2d,each_df.index):
#         array1d = each/(np.max(each)/2)-1
#         new_df.loc[peptide,:] = array1d
#     new_df = new_df.dropna()
#     new_df.to_csv('F:/native_digestion/Uchicago_TMT/tmt_search_0826/normalized_-1_1_'+name+'.tsv',sep='\t')
"""
### normalize tmt within 50% quantile, 0 those data below 25%, didn't work
"""
from scipy import stats
import pickle
df = pd.read_csv('F:/native_digestion/Uchicago_TMT/tmt_search_0826/pep_normalized_int.tsv',delimiter='\t',index_col=0)
pep_list = df.index.tolist()
columns = df.columns.tolist()
file_tmt_dict = {}
for col in columns:
    # df[col] = (df[col]-df[col].mean())/df[col].std()
    pep_tmt_dict = {pep:tmt for pep,tmt in zip(pep_list,df[col]) if tmt!=0}
    new_pep = [k for k in pep_tmt_dict]
    zscore = stats.zscore([v for v in pep_tmt_dict.values()])
    new_pep_tmt_dict = {pep:z for pep,z in zip(new_pep,zscore)}
    file_tmt_dict[col] = new_pep_tmt_dict
pickle.dump(file_tmt_dict,open('F:/native_digestion/Uchicago_TMT/tmt_search_0826/peptide.std.p','wb'))
# df.to_csv('F:/native_digestion/Uchicago_TMT/tmt_search_0826/peptide_std.tsv',sep='\t')
"""

### only keep peptide with intensity above 25 quantile from normalized data
df = pd.read_csv('F:/native_digestion/Uchicago_TMT/tmt_search_0826/peptide_tmt_normalized_bycolumn_01.tsv',delimiter='\t',index_col=0)
newdf = pd.DataFrame(index=df.index, columns=df.columns)
array = df.to_numpy()
for i,j in enumerate(array):
    new_1d_array = np.where(j>np.quantile(j,0.95),1,0)
    newdf.iloc[i,:] = new_1d_array
newdf.to_csv('F:/native_digestion/Uchicago_TMT/tmt_search_0826/tmt_convert_95quantile.tsv',sep='\t')

### plot

import pymannkendall as mk
"""
df = pd.read_csv('F:/native_digestion/Uchicago_TMT/tmt_search_0826/distance_11_13.tsv',delimiter='\t',index_col=0)
new_df = df.loc[:,[each for each in list(df.columns) if 'tryps' in each]]
# new_df = df
mk_result = mk.original_test(new_df.median().tolist())
print(mk_result)
#
plot_array = new_df.to_numpy().flatten()
x = new_df.columns.tolist()*new_df.shape[0]
df_plot = pd.DataFrame(dict(sample=x,dist=plot_array))

fig,axs = plt.subplots(figsize=(10,6))

sns.boxplot(data=df_plot,ax=axs, x='sample',y='dist',linewidth=2.5)
axs.set_xticklabels(list(new_df.columns), fontsize=10,ha="center", rotation=45)
plt.show()
"""