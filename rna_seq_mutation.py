import pandas as pd
from protein_coverage import fasta_reader,read_prot_id_gene, read_description
import re


fasta_path = 'D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_only.fasta'
tsv_path = 'D:/data/pats/results/PXD021865/protein.tsv'

protein_dict = fasta_reader(fasta_path)
description_dict = read_description((fasta_path))

df = pd.read_csv('D:/data/pats/Site_mutsTue Aug 24 20_47_01 2021.csv')
df_protein = pd.read_csv(tsv_path,'r',delimiter='\t')


# sub_df = df[df['Gene Name']=='EIF3CL']
# print (sub_df['AA Mutation'].value_counts())


for protein_id,gene in zip(df_protein['Protein ID'], df_protein['Gene']):
    protein_seq = protein_dict[protein_id]
    sub_df = df[df['Gene Name']==gene]
    mutation_value_sort = sub_df['AA Mutation'].value_counts()
    for each in mutation_value_sort.index.tolist():
        each = each.split('.')[1]
        if each[-3:] == 'del': # site deletion
            original_aa,index = each[0], int(re.sub('[a-zA-Z]','', each))-1
            print (each, original_aa, index)

        elif 'delins' in each:

        else:
            rules = ['?' not in each, '=' not in each, '*' not in each]
            if all(rules): # substitution
                original_aa, mut_aa, index = each[0], each[-1], int(each[1:-1])-1
                print (original_aa,mut_aa,index)


