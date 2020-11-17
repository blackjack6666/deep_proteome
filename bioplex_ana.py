import pickle as ppp
import matplotlib.pyplot as plt
bioplex_cov_dict = ppp.load(open('bioplex_protein_coverage_dict.p','rb'))
coverage_list = [each for each in bioplex_cov_dict.values()]
plt.hist(coverage_list,bins=100,color='black')
plt.xlabel('sequence coverage %')
plt.ylabel('frequency')
plt.show()