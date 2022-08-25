import pandas as pd

excel_peptide_excel = 'F:/native_digestion/Uchicago_TMT/tmt_search_0824/tmt-report/abundance_peptide_GN.tsv'
df = pd.read_csv(excel_peptide_excel,delimiter='\t',index_col=0)
df = df.fillna(0)
columns = df.columns
data_columns = [col for col in df.columns if '_' in col]
print (data_columns)

new_df = pd.DataFrame(index=df.index,columns=df.columns)
# convert log2 int back to original value and normalize as percentage to reference
for each in df.itertuples():
    index = each[0]
    new_df.loc[index,columns[:4]]=each[1:5]
    reference = pow(2,each[5])
    new_df.loc[index,columns[4]]=reference
    for ind, val in enumerate(data_columns):

        new_df.loc[index,val] = pow(2,each[6+ind])/reference
new_df.to_csv('F:/native_digestion/Uchicago_TMT/tmt_search_0824/tmt-report/abundance_peptide_GN_normalized_0824.tsv',sep='\t')

