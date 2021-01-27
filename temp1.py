import pandas as pd
import pickle as p


p_path = 'D:/uic/lab/data/naba/18_2_id_pep_dict.p'
file_id_pep_dict = p.load(open(p_path,'rb'))


# get time-series accumulated covearge for ECM proteins
df = pd.read_csv('D:/uic/lab/data/naba/dash_info_new_1_20.csv')

df = df.drop(columns=['Unnamed: 0', 'Unnamed: 0.1'])
print (df.head())
protein_id_set = df['protein_id'].unique()
print (protein_id_set)
file_set = df['file_name'].unique()
print (file_set)



df_18_2_accumulated = pd.DataFrame()
index_18_2 = ['18_2B05','18_2B1','18_2B2','18_2B4','18_2B18','18_2B20']
for prot in protein_id_set[:2]:

    for idx, val in enumerate(index_18_2):
        index_list = index_18_2[:idx+1]
        print (index_list)
        pep_list = []
        for each_file in index_list:
            try:

                pep_list+=list(file_id_pep_dict[each_file][prot])

            except KeyError:
                continue


        df_18_2_accumulated.loc[val,prot] = coverage

print (df_18_2_accumulated)
