import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, mannwhitneyu, kruskal
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
# print (df_combined.head)
#
# df_combined.to_csv('C:/tools/seqmappdb/demo_data/human/full_parital_combined.csv')

# df = pd.read_csv('C:/tools/seqmappdb/demo_data/full_parital_combined.csv',index_col=0)
df_control = pd.read_excel('F:/native_digestion/chymotrypsin_4_16/search/spearman_atom_dens_5_240min.xlsx',index_col=0).dropna()
# df_control['-log10p'] = -np.log10(df_control['p value'])
df_control['significant_-spearman'] = [1 if spear<0 and p <0.05 else 0 for spear, p in zip(df_control['spearman correlation'],df_control['p value'])]
# df_control = df_control[df_control['spearman correlation']!=0].replace({'p value':{0:0.001}})
df_control = df_control[(df_control['spearman correlation']!=0)&(df_control['p value']!=0)]
df_control['pos_or_neg'] = ['pos' if num >0 else 'neg' for num in df_control['spearman correlation']]
df_control['absolute_spearman'] = np.abs(df_control['spearman correlation'])
#

uniprot_ids = df_combined['queryID'].tolist()
for each in df_control.index.tolist():
    if '>'+each not in uniprot_ids:
        # df_control.at[each,'blast_identity'] = np.random.uniform(0,1) # should be 0, random used for jitter plotting
        df_control.at[each,'blast_identity'] = np.nan
    else:
        df_control.at[each,'blast_identity'] = df_combined.loc[df_combined['queryID']=='>'+each]['identityQuery'].max()
df_control = df_control.dropna()  # drop nan blast identity
plddt_dict = pickle.load(open('D:/data/alphafold_pdb/pLDDT_human_dict.pkl','rb'))
df_control['plddt'] = [np.median(plddt_dict[each]) for each in df_control.index]
df_control['-log_p'] = -np.log10(df_control['p value'].values)

# filter by top and bottom 5% blast identity n
top5, low5 = np.quantile(df_control['blast_identity'],0.75), np.quantile(df_control['blast_identity'],0.25)
df_top5 = df_control[df_control['blast_identity']>=top5]
print (df_top5)
# df_0 = df_control[df_control['blast_identity']==0]
print (df_top5['pos_or_neg'].value_counts())
"""
### plot plddt, blast_identity and spearman correlation
"""
fig,ax = plt.subplots()
sns.scatterplot(data=df_top5,x='blast_identity',y='plddt', size='absolute_spearman', hue='-log_p',
                style='pos_or_neg',
                alpha=0.5, sizes=(50,500))
# sns.scatterplot(data=df_top5,x='spearman correlation',y='-log_p', alpha=0.5,s=50)
# sns.boxplot(data=df_control,x='pos_or_neg',y='blast_identity', palette='Set2')

# sns.violinplot(data=df_control, x='significant_-spearman',y='blast_identity', palette='Set2')
# res = kruskal(df_control[df_control['significant_-spearman']==1]['blast_identity'].values,
#                    df_control[df_control['significant_-spearman']==0]['blast_identity'].values)
# print (res)
# ax.yaxis.grid(True)

# coef = np.polyfit(df_control['blast_identity'],df_control['plddt'],1)
# poly1d_fn = np.poly1d(coef)
# ax.plot(df_control['blast_identity'],poly1d_fn(df_control['blast_identity']),'--k')
ax.set_ylim(0,110)
# ax.get_legend().remove()
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.axhline(y=-np.log10(0.05), color='k', linestyle='--')
# plt.show()
# r, p = pearsonr(df_control.dropna()['blast_identity'],df_control.dropna()['plddt'])
# print (f'r:{r},p:{p}')
# print (df_control.loc[(df_control['spearman correlation']>0)&(df_control['p value']<0.05)]['blast_identity'].mean())
"""
### compare cleavage density data with probe having different size
"""
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
"""
### get overlap from sasa, atom density and distance analysis
"""
from tsv_reader import venn_diagram_gen
df_sasa_spearman = pd.read_excel('D:/data/native_protein_digestion/12072021/control/sasa_spearman_10_240min.xlsx', index_col=0).dropna()
sasa_candidate = df_sasa_spearman[(df_sasa_spearman['spearman correlation']<0)].index.values

df_atom_spearman = pd.read_excel('D:/data/native_protein_digestion/12072021/control/atom_spearman_10_240min.xlsx',index_col=0).dropna()
atoms_candidate = df_atom_spearman[(df_atom_spearman['spearman correlation']<0)].index.values

df_dist_spearman = pd.read_excel('D:/data/native_protein_digestion/12072021/control/dist_spearman_10_240min.xlsx',index_col=0).dropna()
dist_candidate = df_dist_spearman[(df_dist_spearman['spearman correlation']<0)].index.values

print (np.intersect1d(np.intersect1d(atoms_candidate,dist_candidate),sasa_candidate))

ven_dict = {'sasa':sasa_candidate,'atoms_density':atoms_candidate,'distance':dist_candidate}
venn_diagram_gen(ven_dict,'proteins with spearman < 0')
"""

### correlate sasa, density and distance
from scipy.stats import pearsonr
df_sasa = pd.read_excel('D:/data/native_protein_digestion/12072021/control/sasa.xlsx',index_col=0)
df_density = pd.read_excel('D:/data/native_protein_digestion/12072021/control/cov_KR_density_15A.xlsx', index_col=0)
df_distance = pd.read_excel('D:/data/native_protein_digestion/12072021/control/cov_dist_unique.xlsx', index_col=0)

df_corr = pd.DataFrame(index=df_sasa.index,columns=df_sasa.columns)
x, y = [], []
for ind in df_sasa.index:
    sasa_vals = df_sasa.loc[ind,:].values
    sasa_clean = sasa_vals[np.isfinite(sasa_vals)]
    density_vals = df_density.loc[ind,:].values
    density_clean = density_vals[np.isfinite(density_vals)]
    dist_vals = df_distance.loc[ind,:].values
    dist_clean = dist_vals[np.isfinite(dist_vals)]
    if all([len(sasa_clean)>=3, len(density_clean)>=3, len(dist_clean)>=3]):
        x.append(pearsonr(sasa_clean,density_clean))
        y.append(pearsonr(sasa_clean,dist_clean))

significant_x_y = [(val[0],y[ind][0]) for ind,val in enumerate(x) if val[1]<0.05 and y[ind][1]<0.05]
non_signi_x_y = [(val[0],y[ind][0]) for ind,val in enumerate(x) if val[1]>=0.05 or y[ind][1]>=0.05]
x_, y_ = zip(*significant_x_y)
nonsig_x, nonsig_y = zip(*non_signi_x_y)
plt.scatter(x=x_,y=y_,s=8,alpha=0.6, c='red', label='p val < 0.05')
plt.scatter(x=nonsig_x,y=nonsig_y,s=8,alpha=0.3, c='grey',label='p val >= 0.05' )
plt.xlabel('Pearson corr between sasa and atom density')
plt.ylabel('Pearson corr between sasa and center distance')
plt.axhline(0,-1,1,c='k',linestyle='--')
plt.axvline(0,-1,1,c='k',linestyle='--')
plt.xlim(-1.1,1.1)
plt.ylim(-1.1,1.1)
plt.legend(loc='upper right',framealpha=0.5)
plt.show()

# df_plot = pd.DataFrame(dict(pearson=[each[0] for each in x]+[each[0] for each in y],x_or_y=['x']*len(x)+['y']*len(y)))
# sns.violinplot(x='x_or_y',y='pearson',data=df_plot)
# plt.show()