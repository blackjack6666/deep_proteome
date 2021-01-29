import pandas as pd
import pickle as p
from protein_coverage import fasta_reader
import numpy as np
from tsv_reader import venn_diagram_gen


p_path = 'D:/uic/lab/data/naba/18_2_id_pep_dict.p'
file_id_pep_dict = p.load(open(p_path,'rb'))

id_pep_dict_1805 = file_id_pep_dict['18_2B05']


fasta_path = 'D:/uic/lab/data/naba/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
protein_dict = fasta_reader(fasta_path)

# get time-series accumulated covearge for ECM proteins
df = pd.read_csv('D:/uic/lab/data/naba/dash_info_new_1_20.csv')

df = df.drop(columns=['Unnamed: 0', 'Unnamed: 0.1'])
print (df.head())
protein_id_set = df['protein_id'].unique()
print (protein_id_set)
file_set = df['file_name'].unique()

print (len(protein_id_set))

df_18_2_accumulated = pd.DataFrame()
index_18_2 = ['18_2B05','18_2B1','18_2B2','18_2B4','18_2B18','18_2B20']
for prot in protein_id_set:
    protein_seq = protein_dict[prot]
    np_array = np.zeros(len(protein_seq))
    for idx, val in enumerate(index_18_2):
        index_list = index_18_2[:idx+1]
        # print (index_list)
        pep_list = []
        for each_file in index_list:
            if prot in file_id_pep_dict[each_file]:
                pep_list+=list(file_id_pep_dict[each_file][prot])

            else:
                pep_list+=[]

        peptide_set = set(pep_list)

        for each in peptide_set:
            pep_loc = protein_seq.find(each)
            pep_end_loc = pep_loc+len(each)
            np_array[pep_loc:pep_end_loc]+=1
        coverage = np.count_nonzero(np_array)/len(np_array)*100

        df_18_2_accumulated.loc[val,prot] = coverage

for column in df_18_2_accumulated.columns:
    column_list = df_18_2_accumulated[column].tolist()
    if column_list[0]!=column_list[1]:
        print (column)

