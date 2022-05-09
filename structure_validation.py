import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import numpy as np
import pickle
from glob import glob

### combine full and paritial uniprot-PDB mapping csv file
"""
header_names = ['uniprot',
                'PDB_chain',
                'identity_against_uniprot',
                'identity_against_PDBchain',
                'seqLength_uniprot',
                'coverage_uniprot',
                'coverage_PDBchain']
df_full = pd.read_csv('C:/tools/seqmappdb/human/fully_covered_unique_PDB.csv')

df_parital = pd.read_csv('C:/tools/seqmappdb/human/partial_covered_nonoverlapping_PDB.csv')
df_combined = pd.concat([df_full,df_parital],ignore_index=True)
print (df_combined.head)
#
# df_combined.to_csv('C:/tools/seqmappdb/demo_data/human/full_parital_combined.csv')

# df = pd.read_csv('C:/tools/seqmappdb/demo_data/full_parital_combined.csv',index_col=0)
df_control = pd.read_excel('D:/data/native_protein_digestion/12072021/control/KR_atoms_linear_reg.xlsx',index_col=0).dropna()
# df_control['-log10p'] = -np.log10(df_control['p value'])
df_control['significant_-slope'] = [1 if slope<0 and p <0.05 else 0 for slope, p in zip(df_control['slope'],df_control['p_val'])]
df_control = df_control[df_control['slope']!=0]
#
uniprot_ids = df_combined['queryID'].tolist()
for each in df_control.index.tolist():
    if '>'+each not in uniprot_ids:
        df_control.at[each,'blast_identity'] = 0
    else:
        df_control.at[each,'blast_identity'] = df_combined.loc[df_combined['queryID']=='>'+each]['identityQuery'].max()

plddt_dict = pickle.load(open('D:/data/alphafold_pdb/pLDDT_human_dict.pkl','rb'))
df_control['plddt'] = [np.median(plddt_dict[each]) for each in df_control.index]
print (df_control.head)
fig,ax = plt.subplots()
sns.scatterplot(data=df_control,x='blast_identity',y='plddt', hue='significant_-slope',
                alpha=0.6,palette='Set2')

coef = np.polyfit(df_control['blast_identity'],df_control['plddt'],1)
poly1d_fn = np.poly1d(coef)
ax.plot(df_control['blast_identity'],poly1d_fn(df_control['blast_identity']),'--k')
# ax.set_ylim(-0.1,1.1)
# ax.get_legend().remove()
plt.show()
r, p = pearsonr(df_control.dropna()['blast_identity'],df_control.dropna()['plddt'])
print (f'r:{r},p:{p}')
# print (df_control.loc[(df_control['spearman correlation']>0)&(df_control['p value']<0.05)]['blast_identity'].mean())
"""
### compare cleavage density data with probe having different size
import pymannkendall as mk
file_list = glob('D:/data/native_protein_digestion/12072021/control/cov_KR_density_*.xlsx')
print (file_list)
color_list = ['black','brown','cyan','gold','green']

fig,axs = plt.subplots(figsize=(10,6))
for each, color in zip(file_list,color_list):
    df = pd.read_excel(each,index_col=0)
    median_list = np.log2(df.median().tolist())
    mk_result = mk.original_test(df.median().tolist())
    print (each.split('\\')[-1], mk_result)
    x = range(1, len(df.columns) + 1)
    axs.plot(x, median_list, linestyle='-', color=color, linewidth=3, label=each.split('\\')[-1].split('_')[-1].split('.')[0])
axs.set_xticks(x)
axs.set_xticklabels(list(df.columns), fontsize=12,ha="center", rotation=45)
axs.set_ylabel('log2(cleavage density)',fontsize=12)
axs.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.show()