import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

pd.set_option('display.max_columns', None)

df = pd.read_excel("D:/data/deep_proteome/9_20_identified_cov_total.xlsx")
df['protein_num'] = range(len(df))
df = df.set_index('Unnamed: 0')

# df_split = np.split(df,[200],axis=0)[0]
# df_split = df_split[df_split.columns[:-1]]
# print (df_split.head())
#
# fig, ax = plt.subplots(1,1, figsize=(40,20))
# # for each_cond in df_split.columns[1:-1]:
# #     print (each_cond)
# #     ax.plot(df_split['protein_num'],df_split[each_cond], '-', c=np.random.rand(3,), label=each_cond)
#
# g = sns.heatmap(df[df.columns[:-1]],ax=ax,yticklabels=False)
# #
# plt.xticks(rotation=45)
# plt.show()

#