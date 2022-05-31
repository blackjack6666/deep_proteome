import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, mannwhitneyu, kruskal
import numpy as np
import pickle
from glob import glob
from tsv_reader import protein_tsv_reader, get_unique_peptide, map_k_r
from protein_coverage import fasta_reader
import os
import umap
import hdbscan
import sklearn.cluster as cluster
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
from sklearn.decomposition import PCA

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
"""
from scipy.stats import pearsonr, wilcoxon

df_sasa = pd.read_excel('D:/data/native_protein_digestion/12072021/control/sasa.xlsx',index_col=0)
df_density = pd.read_excel('D:/data/native_protein_digestion/12072021/control/cov_KR_density_15A.xlsx', index_col=0)
df_distance = pd.read_excel('D:/data/native_protein_digestion/12072021/control/cov_dist_unique.xlsx', index_col=0)

df_corr = pd.DataFrame(index=df_sasa.index,columns=df_sasa.columns)
x, y = [], []
result_dict = {}
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

# significant_x_y = [(val[0],y[ind][0]) for ind,val in enumerate(x) if val[1]<0.05 and y[ind][1]<0.05]
# non_signi_x_y = [(val[0],y[ind][0]) for ind,val in enumerate(x) if val[1]>=0.05 or y[ind][1]>=0.05]
# x_, y_ = zip(*significant_x_y)
# nonsig_x, nonsig_y = zip(*non_signi_x_y)
# plt.scatter(x=x_,y=y_,s=8,alpha=0.6, c='red', label='p val < 0.05')
# plt.scatter(x=nonsig_x,y=nonsig_y,s=8,alpha=0.3, c='grey',label='p val >= 0.05' )
# plt.xlabel('Pearson corr between sasa and atom density')
# plt.ylabel('Pearson corr between sasa and center distance')
# plt.axhline(0,-1,1,c='k',linestyle='--')
# plt.axvline(0,-1,1,c='k',linestyle='--')
# plt.xlim(-1.1,1.1)
# plt.ylim(-1.1,1.1)
# plt.legend(loc='upper right',framealpha=0.5)
# plt.show()

# df_plot = pd.DataFrame(dict(pearson=[each[0] for each in x]+[each[0] for each in y],x_or_y=['x']*len(x)+['y']*len(y)))
# sns.violinplot(x='x_or_y',y='pearson',data=df_plot)
# plt.show()
"""

### implement amino acid exposure data from structuremap (full exposure, value is number of neiboring amino acids)
# https://github.com/MannLabs/structuremap_analysis/blob/master/data_analysis_structuremap.ipynb
"""
fasta_file = 'D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_tr.fasta'
protein_dict = fasta_reader(fasta_file)

# df_full_exp = pd.read_csv('D:/data/alphafold_pdb/full_sphere_exposure.csv') 

## split big csv into chunks and save individually for higher reading speed

# pd_dicts = {key: df_full_exp.loc[value] for key, value in df_full_exp.groupby("protein_id").groups.items()}
# for k in pd_dicts:
#     print (k)
#     pd_dicts[k].to_csv(os.path.join('D:/data/alphafold_pdb/full_sphere_expo_split/',k+'.csv'))
# print (df_full_exp.head)

protein_tsv = 'D:/data/native_protein_digestion/12072021/control/combined_protein.tsv'
protein_list = protein_tsv_reader(protein_tsv, protein_column=3)
sub_protein_dict = {prot:protein_dict[prot] for prot in protein_list}

base_path = 'D:/data/native_protein_digestion/12072021/control/'
folders = [base_path + folder for folder in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, folder))]
time_points = [each.split('/')[-1] for each in folders]
pep_path_list = [each + '/peptide.tsv' for each in folders]
psm_path_list = [each + '/psm.tsv' for each in folders]
unique_peptide_dict = get_unique_peptide(pep_path_list)

df_native_exposure = pd.DataFrame(index=protein_list,columns=time_points)

for pep_tsv in pep_path_list:
    print (pep_tsv)
    peptide_list = unique_peptide_dict[pep_tsv.split('/')[-2]]
    freq_array_dict, freq_array_index_dict = map_k_r(peptide_list,sub_protein_dict)

    for prot in protein_list:
        print (prot)
        try:
            df_expo = pd.read_csv('D:/data/alphafold_pdb/full_sphere_expo_split/'+prot+'.csv')

            freq_index = freq_array_index_dict[prot]
            if len(freq_index)!=0: # has at least one peptide matched
                total_expo = sum([df_expo.at[each,'nAA_24_180_pae'] for each in freq_index])
                ave_expo = total_expo/len(freq_index)
                df_native_exposure.at[prot,pep_tsv.split('/')[-2]] = ave_expo
            else:
                df_native_exposure.at[prot,pep_tsv.split('/')[-2]] = np.nan
        except FileNotFoundError:
            df_native_exposure.at[prot, pep_tsv.split('/')[-2]] = np.nan
df_native_exposure.to_excel('D:/data/native_protein_digestion/12072021/control/aa_exposure_structuremap.xlsx')

"""

### correlation
from matplotlib.collections import PatchCollection

df_sasa = pd.read_excel('D:/data/native_protein_digestion/12072021/control/sasa.xlsx',index_col=0)
df_density = pd.read_excel('D:/data/native_protein_digestion/12072021/control/cov_KR_density_15A.xlsx', index_col=0)
df_distance = pd.read_excel('D:/data/native_protein_digestion/12072021/control/cov_dist_unique.xlsx', index_col=0)
df_atom_exposure = pd.read_excel('D:/data/native_protein_digestion/12072021/control/aa_exposure_structuremap.xlsx',index_col=0)
"""
df_atom_exposure_filter = pd.DataFrame(columns=df_distance.columns)
for tp in df_distance.itertuples():
    prot_id, values = tp[0], tp[1:]
    if sum(~np.isnan(values)) >=3: # filter out those with less than 3 data points
        df_atom_exposure_filter.loc[prot_id,:] = values


# df_corr = pd.DataFrame(columns=['pearson_corr', 'p_val'])
# for ind in df_sasa.index:
#     for ind2 in df_sasa.index:
#         if ind!=ind2:
#             sasa_vals = df_sasa.loc[ind,:].values
#             sasa_clean = sasa_vals[np.isfinite(sasa_vals)]
#             density_vals = df_density.loc[ind,:].values
#             density_clean = density_vals[np.isfinite(density_vals)]
#         # dist_vals = df_distance.loc[ind,:].values
#         # dist_clean = dist_vals[np.isfinite(dist_vals)]
#         # atom_exposure = df_atom_exposure.loc[ind,:].values
#         # atom_exposure_clean = atom_exposure[np.isfinite(atom_exposure)]
#             if all([len(sasa_clean)>=3, len(density_clean)>=3, len(dist_clean)>=3, len(atom_exposure_clean)>=3]):
#                 df_corr.loc[ind,:] = pearsonr(sasa_clean,density_clean)
#
# df_corr['-log10p'] = -np.log10(df_corr.p_val.astype(np.float64))
fig,ax = plt.subplots(figsize=(12,10))
# df_atom_exposure_filter = df_atom_exposure_filter.astype(float)
df_corr = df_atom_exposure_filter.T.corr()
print (df_corr)
# sns.heatmap(df_corr,ax=ax,xticklabels=[], yticklabels=[],cmap='viridis',cbar_kws={"shrink": 0.5})
# plt.show()
"""
### umap clustering
"""
sasa_spearman_df = pd.read_excel('D:/data/native_protein_digestion/12072021/control/sasa_spearman_10_240min.xlsx')
denstiy_spearman_df = pd.read_excel('D:/data/native_protein_digestion/12072021/control/atom_spearman_10_240min.xlsx')
distance_spearman_df = pd.read_excel('D:/data/native_protein_digestion/12072021/control/dist_spearman_10_240min.xlsx')
sasa_spearman,denstiy_spearman,distance_spearman = sasa_spearman_df['spearman correlation'].values, \
                                                   denstiy_spearman_df['spearman correlation'].values, \
                                                   distance_spearman_df['spearman correlation'].values
sasa_pval, density_pval, distance_pval = sasa_spearman_df['p value'].values, \
                                         denstiy_spearman_df['p value'].values, \
                                         distance_spearman_df['p value'].values


aver_spearman = np.sum([sasa_spearman,denstiy_spearman,distance_spearman],axis=0)/3
min_pval = np.min([sasa_pval,density_pval,distance_pval],axis=0)
bool_spearman_negative = np.where((aver_spearman<0)&(min_pval<0.05),1,0)  # get average spearman <0 one of the pval<0.05
bool_spearman_positive = np.where((aver_spearman>0)&(min_pval<0.05),-1,0)
color_column = np.sum([bool_spearman_negative,bool_spearman_positive],axis=0)


df_plot = pd.concat([df_sasa,df_distance,df_density],axis=1)
df_plot_fill = df_plot.fillna(0)
# df_plot_fill = df_plot.dropna()
# print (df_plot_fill.head)

## PCA filter out noisy data
pca = PCA(df_plot_fill.shape[1])

df_pca = pca.fit_transform(df_plot_fill)  # fill na with extreme values
print (df_pca.shape)
print('Explained variation per principal component: {}'.format(pca.explained_variance_ratio_))
## variance explained by PC1 is too low, not a good idea to use PCA eigenvalues

## umap visualization/dimension reduction
# clusterable_embedding = umap.UMAP(
#     n_neighbors=30,
#     min_dist=0.0,
#     n_components=2,
#     random_state=42,
# ).fit_transform(df_plot_fill)

# differentiate colors by bool_spearman condition
# plt.scatter(clusterable_embedding[:, 0], clusterable_embedding[:, 1],c=color_column, cmap='plasma',
#              s=0.3,alpha=0.8)
# plt.show()

## HDBSCAN clustering
# labels = hdbscan.HDBSCAN(
#     min_samples=10,
#     min_cluster_size=60,
# ).fit_predict(clusterable_embedding)
#
# clustered = (labels >= 0)
#
# plt.scatter(clusterable_embedding[~clustered, 0],
#             clusterable_embedding[~clustered, 1],
#             color=(0.5, 0.5, 0.5),
#             s=0.3,
#             alpha=0.5)
# plt.scatter(clusterable_embedding[clustered, 0],
#             clusterable_embedding[clustered, 1],
#             c=labels[clustered],
#             s=0.3,
#             cmap='Spectral')
#
# plt.show()
"""
### protein fragments length analysis in native digestion
"""
df_cleav_index = pd.read_excel('D:/data/native_protein_digestion/12072021/control/cleavage_index_4_24.xlsx',index_col=0)
selected_cols = df_cleav_index.columns[1:-2]
df_new = pd.DataFrame(index=df_cleav_index.index,columns=selected_cols)

# print (selected_cols)
for ind in df_cleav_index.index:
    seq_len = len(df_cleav_index.at[ind,'sequence'])
    for each_col in selected_cols:
        cells = df_cleav_index.loc[ind, '0010min':each_col].tolist()  # multiple columns

        cleav_array = ''.join([each for each in cells]).replace('[',' ').replace(']',' ').replace('\n','')
        np_array = np.append(np.fromstring(cleav_array,dtype=int, sep=' '),seq_len)
        unique_sort_array = np.sort(np.unique(np_array))
        substract_array = np.insert(unique_sort_array, 0, 0)[:-1]
        # print (unique_sort_array,substract_array)

        frag_ratio = np.subtract(unique_sort_array,substract_array)/seq_len # normalized by length of protein
        max_ind = np.argmax(frag_ratio)
        max_idx_first, max_idx_latter = substract_array[max_ind],unique_sort_array[max_ind]
        print (max_idx_first,max_idx_latter)
        # df_new.at[ind,each_col] = np.max(frag_ratio)
        df_new.at[ind,each_col] = (max_idx_first,max_idx_latter)
# df_new.to_excel('D:/data/native_protein_digestion/12072021/control/digestion_max_peptide_relative_length.xlsx')
df_new.to_excel('D:/data/native_protein_digestion/12072021/control/digestion_max_peptide_index.xlsx')
"""

### filter by blast identity
# df_blast_identity = pd.read_csv('C:/tools/seqmappdb/demo_data/full_parital_combined.csv')
# quantile95 = np.quantile(df_blast_identity['identity_against_uniprot'], 0.95)
#
# filtered_df = df_blast_identity[df_blast_identity['identity_against_uniprot']>quantile95]


### text cloud
from collections import Counter
with open('D:/data/ASMS/New Text Document.txt','r') as f_o:
    text_list = f_o.readline().split(' ')

counted = Counter(text_list)
df = pd.DataFrame({'Weight':[v for v in counted.values()], 'Word':[k for k in counted.keys()]})
df.to_csv('D:/data/ASMS/background.csv')