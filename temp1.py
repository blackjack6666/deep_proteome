import pandas as pd

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
        coverage = 0
        for each_file in index_list:
            try:

                coverage+=df[(df['protein_id']==prot)&(df['file_name']==each_file)]['coverage'].values[0]
                print (coverage)
            except IndexError:
                continue
        df_18_2_accumulated.loc[val,prot] = coverage

print (df_18_2_accumulated)
