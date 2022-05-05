import pandas as pd


### combine full and paritial uniprot-PDB mapping csv file
# header_names = ['uniprot',
#                 'PDB_chain',
#                 'identity_against_uniprot',
#                 'identity_against_PDBchain',
#                 'seqLength_uniprot',
#                 'coverage_uniprot',
#                 'coverage_PDBchain']
# df_full = pd.read_csv('C:/tools/seqmappdb/demo_data/fully_covered_unique_PDB.csv',names=header_names)
# df_parital = pd.read_csv('C:/tools/seqmappdb/demo_data/partial_covered_nonoverlapping_PDB.csv',names=header_names)
# df_combined = pd.concat([df_full,df_parital],ignore_index=True)
#
# df_combined.to_csv('C:/tools/seqmappdb/demo_data/full_parital_combined.csv')

df = pd.read_csv('C:/tools/seqmappdb/demo_data/full_parital_combined.csv',index_col=0)
df_density_linear_reg = pd.read_excel('D:/data/native_protein_digestion/12072021/control/KR_atoms_linear_reg.xlsx',index_col=0)

uniprot_ids = df['uniprot'].tolist()
for each in df_density_linear_reg.index.tolist():
    if '>'+each in uniprot_ids:
        print (each)