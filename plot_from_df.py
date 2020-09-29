import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from protein_coverage import fasta_reader2
from calculations_and_plot import length_KR_ratio, specific_aa_length_ratio

fasta_path = 'D:/data/proteome_fasta/uniprot-proteome_UP000005640.fasta'
protein_dict = fasta_reader2(fasta_path)
kr_len_ratio = specific_aa_length_ratio(protein_dict,spec_aa='FWY')


protein_len_dict = {each:len(protein_dict[each]) for each in protein_dict}

pd.set_option('display.max_columns', None)

df = pd.read_excel("D:/data/deep_proteome/9_20_identified_cov_total.xlsx")


df['protein_num'] = range(len(df))
df = df.set_index('Unnamed: 0')
protein_list = df.index.tolist()

kr_len_ratio_list = [kr_len_ratio[each] for each in protein_list]
df['kr_len_ratio'] = kr_len_ratio_list
# kr length ratio plot

# plt.hist(kr_len_ratio_list,bins=100,color='black')
# plt.xlabel('FWY freq/protein length')
# plt.ylabel('frequency')
# plt.title('The ratio of FWY frequency to protein length for all proteins identified')
# plt.show()

# protein_len_list = [protein_len_dict[each] for each in protein_list]
# df['protein_len'] = protein_len_list

df = df.sort_values('kr_len_ratio')
print (df[-100:])
# # df_split = np.split(df,[200],axis=0)[0]
# # df_split = df_split[df_split.columns[:-1]]
# # print (df_split.head())

fig, ax = plt.subplots(1,1, figsize=(40,20))
# for each_cond in df_split.columns[1:-1]:
#     print (each_cond)
#     ax.plot(df_split['protein_num'],df_split[each_cond], '-', c=np.random.rand(3,), label=each_cond)

g = sns.heatmap(df[df.columns[:-2]],ax=ax,yticklabels=False)
#
plt.xticks(rotation=60)
plt.show()

# # charge plot
# df_charge = pd.read_excel('D:/data/deep_proteome/9_20_peptide_charge.xlsx')
# fig, ax = plt.subplots(1,1, figsize=(40,20))
# g = sns.heatmap(df_charge,ax=ax,yticklabels=False)
# plt.xticks(rotation=60)
# plt.show()