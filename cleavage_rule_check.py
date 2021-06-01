from tsv_reader import peptide_counting, id_pep_from_peptsv
from protein_coverage import fasta_reader2
import re
import pandas as pd
from collections import defaultdict
import numpy as np
from data_preprocess_3 import polymer_gen
import matplotlib.pyplot as plt

fasta_path = 'D:/data/proteome_fasta/uniprot-proteome_UP000005640.fasta'
pep_tsv = 'D:/data/deep_proteome/non_specfic_search/ct_4h/peptide.tsv'
psm_tsv = 'D:/data/deep_proteome/non_specfic_search/ct_4h/psm.tsv'

protein_dict = fasta_reader2(fasta_path)
id_peptide_dict = id_pep_from_peptsv(pep_tsv)
df_psm = pd.read_csv(psm_tsv,sep='\t')
psm_score_dict = defaultdict(list)  # psm and hyperscore
for each_row in df_psm.itertuples():
    psm_score_dict[each_row[3]].append(each_row[16])
pep_average_score_dict = {each:np.mean(psm_score_dict[each]) for each in psm_score_dict}

# check if pep fits customized regex cleavage rule
rule = r'(?<!D)[YFLWM](?!EE|ED|DE|DD|PE|PD)'

matched_peptide = 0
no_matched_peptide = 0
matched_peptide_score_dict, no_matched_peptide_score_dict = {}, {}
for prot in id_peptide_dict:
    prot_seq = protein_dict[prot]
    for pep in id_peptide_dict[prot]:
        pep_loc = prot_seq.find(pep) - 1
        if pep_loc!=-2:
            pep_end_loc = pep_loc+len(pep)
            start_31mer, end_31_mer = polymer_gen(prot_seq,pep_loc), polymer_gen(prot_seq,pep_end_loc)
            start_match, end_match = [m.end() for m in re.finditer(rule,start_31mer)], [m.end() for m in re.finditer(rule,end_31_mer)]
            if 16 in start_match and 16 in end_match:  # start and end position of a peptide matches with cleavage rule
                matched_peptide+=1
                matched_peptide_score_dict[pep] = pep_average_score_dict[pep]
            elif (16 not in start_match) and (16 not in end_match):
            #
                print (prot,pep, start_31mer, start_match,end_31_mer, end_match)
                no_matched_peptide+=1
                no_matched_peptide_score_dict[pep] = pep_average_score_dict[pep]


matched_hyper_values = [v for v in matched_peptide_score_dict.values()]
no_matched_hyper_values = [v for v in no_matched_peptide_score_dict.values()]

print (np.mean(matched_hyper_values), np.mean(no_matched_hyper_values))

# plt.hist(matched_hyper_values, bins=50, alpha=0.5, label='regex-matched psm, average score=20.6',log=True)
# plt.hist(no_matched_hyper_values,bins=50, alpha=0.5, color='pink', label='regex-filtered psm, average score=16.2',log=True)
# plt.xlabel('Hyperscore')
# plt.ylabel('Frequency')
# plt.legend()
# plt.show()


from scipy import stats
print (stats.ttest_ind(matched_hyper_values, no_matched_hyper_values))
