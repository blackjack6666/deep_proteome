import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from protein_coverage import fasta_reader
from scipy.stats import pearsonr
import numpy as np
from tsv_reader import combined_proteintsv_map,venn_diagram_gen


combined_tsv = 'F:/deion_biotin/05152022/search/combined_protein.tsv'
# df = pd.read_csv(combined_tsv,delimiter='\t', header=0, index_col=False,usecols=['Protein ID','Gene',
#                                                     'Fluc_TMV_1 Spectral Count',
#                                                     'Fluc_TMV_2 Spectral Count',
#                                                     'Fluc_ctrl_1 Spectral Count',
#                                                     'Fluc_ctrl_2 Spectral Count',
#                                                     'Hoil_TMV_1 Spectral Count',
#                                                     'Hoil_TMV_2 Spectral Count',
#                                                     'Hoil_ctrl_1 Spectral Count',
#                                                     'Hoil_ctrl_2 Spectral Count'])

### NSAF calculation
"""
protein_dict = fasta_reader('F:/deion_biotin/uniprot-proteome_UP000000589.fasta')
# df.index = df['Protein ID']
df_filter = df.iloc[:,2:]
df_nsaf = pd.DataFrame(columns=df_filter.columns, index=df['Gene'])
proteinid_list = df['Protein ID'].tolist()
for column in df_filter.columns:
    spectra_count_list = df_filter[column].tolist()
    nsaf_sum = sum([val/len(protein_dict[proteinid_list[ind]]) for ind,val in enumerate(spectra_count_list)])
    df_nsaf[column] = [val/len(protein_dict[proteinid_list[ind]])/nsaf_sum for ind, val in enumerate(spectra_count_list)]
"""

### heatmap
# df.index = df.Gene
# df = df.iloc[:,1:]
# print (df.T)
# fig,ax = plt.subplots(figsize=(20,120))
# sns.set(font_scale = 1)
# sns.heatmap(data=df_nsaf,ax=ax,xticklabels=1,yticklabels=1, annot=True,cmap='viridis')
# ax.set_xticklabels([column.replace('Spectral Count','NSAF') for column in df_nsaf.columns])
# plt.xticks(rotation=45)
# plt.show()


# df_sample1, df_sample2 = df_nsaf.iloc[:,[0,2,4,6]], df_nsaf.iloc[:,[1,3,5,7]]
# result = [pearsonr(sample1[::],sample2[::]) for sample1, sample2 in zip(df_sample1.itertuples(index=False),df_sample2.itertuples(index=False))]
# plt.boxplot([each[0] for each in result if type(each[0])==np.float64])
# plt.tick_params(
#     axis='x',          # changes apply to the x-axis
#     which='both',      # both major and minor ticks are affected
#     bottom=False,      # ticks along the bottom edge are off
#     top=False,         # ticks along the top edge are off
#     labelbottom=False) # labels along the bottom edge are off
# plt.ylabel('Pearson correlation coefficient r')
# plt.title('Pearson correlation of NSAF bewteen sample 1/2')
# plt.show()
# print ([each[0] for each in result if type(each[0])==np.float64])

### venn_diagram
info_dict = combined_proteintsv_map(combined_tsv)
print ([k for k in info_dict.keys()])
venn_dict = {'Hoil_ctrl_1':[k for k in info_dict['Hoil_ctrl_1'].keys()],
             'Hoil_ctrl_2':[k for k in info_dict['Hoil_ctrl_2'].keys()]}
venn_diagram_gen(venn_dict,'protein IDs from Hoil_ctrl')