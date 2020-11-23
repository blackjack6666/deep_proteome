import pickle as ppp
import matplotlib.pyplot as plt
import statistics
import pandas as pd
import numpy as np
# bioplex_cov_dict = ppp.load(open('bioplex_protein_coverage_dict.p','rb'))
# coverage_list = [each for each in bioplex_cov_dict.values()]
# print (statistics.mean(coverage_list),statistics.median(coverage_list))
# print (len(coverage_list))
# plt.hist(coverage_list,bins=100,color='black')
# plt.xlabel('sequence coverage %')
# plt.ylabel('frequency')
# plt.show()

# matrisome db data

matrisome_cov_csv = 'D:/uic/lab/data/naba/matrisome coverage.xlsx'

df = pd.read_excel(matrisome_cov_csv)
ecm_class = set(df['class'].tolist())

df = df.drop_duplicates()
print (df.shape)

# subplots
rows = [0,0,0,1,1,1]
colums = [0,1,2,0,1,2]
fig,axs = plt.subplots(2,3,figsize=(10,5))
for each,row,col in zip(ecm_class,rows,colums):
    sub_df = df[df['class']==each]
    cov_list = sub_df['sequence_cov'].tolist()
    print (statistics.mean(cov_list),statistics.median(cov_list))
    axs[row,col].hist(cov_list,bins=50,color='black')
    axs[row,col].set_xticks(np.arange(0,1.2,0.2))
    #axs[row,col].set_ylim(0,3500)
    axs[row,col].set_xlabel('sequence coverage')
    axs[row,col].set_ylabel('frequency')
    axs[row,col].set_title(each)

plt.show()
# cov_list = df['cov'].tolist()
# print (len(cov_list))
# print (statistics.mean(cov_list),statistics.median(cov_list))
# plt.hist(cov_list,bins=100,color='black')
# plt.xlabel('sequence coverage %')
# plt.ylabel('frequency')
# plt.show()