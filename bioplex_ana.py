import pickle as ppp
import matplotlib.pyplot as plt
import statistics
import pandas as pd

# bioplex_cov_dict = ppp.load(open('bioplex_protein_coverage_dict.p','rb'))
# coverage_list = [each for each in bioplex_cov_dict.values()]
# print (statistics.mean(coverage_list),statistics.median(coverage_list))
# print (len(coverage_list))
# plt.hist(coverage_list,bins=100,color='black')
# plt.xlabel('sequence coverage %')
# plt.ylabel('frequency')
# plt.show()

# matrisome db data

matrisome_cov_csv = 'D:/data/Naba_deep_matrisome/seq_cov_matrisomedb.csv'

df = pd.read_csv(matrisome_cov_csv)
cov_list = df['cov'].tolist()
print (len(cov_list))
print (statistics.mean(cov_list),statistics.median(cov_list))
plt.hist(cov_list,bins=100,color='black')
plt.xlabel('sequence coverage %')
plt.ylabel('frequency')
plt.show()