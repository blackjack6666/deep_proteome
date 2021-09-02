import pandas as pd
from protein_coverage import fasta_reader,read_prot_id_gene, read_description
import re
import pickle as ppp



fasta_path = 'D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_only.fasta'
tsv_path = 'D:/data/pats/results/PXD022915_standard_sp_only/protein.tsv'

protein_dict = fasta_reader(fasta_path)
description_dict = read_description((fasta_path))
### get mutant protein seq

df = pd.read_csv('D:/data/pats/Site_muts_stomach_cancer.csv')
df_protein = pd.read_csv(tsv_path,'r',delimiter='\t')


# sub_df = df[df['Gene Name']=='EIF3CL']
# print (sub_df['AA Mutation'].value_counts())
mut_protein_dict = {}

counter = 0
for protein_id,gene in zip(df_protein['Protein ID'], df_protein['Gene']):
    protein_seq = protein_dict[protein_id]
    sub_df = df[df['Gene Name']==gene]
    mutation_value_sort = sub_df['AA Mutation'].value_counts()
    del_aa_index = []
    subs_aa_index = []
    for each in mutation_value_sort.index.tolist():
        each = each.split('.')[1]

        if each[-3:] == 'del' and '_' not in each: # site deletion

            original_aa,index = each[0], int(re.sub('[a-zA-Z]','', each))-1
            if index<len(protein_seq):
                if protein_seq[index] == original_aa:
                    del_aa_index.append((original_aa,index))
                else:
                    print (protein_id, index, 'deletion not match')
            # print (each, original_aa, index)


        # elif 'dup' in each:


        else:
            rules = ['?' not in each, '=' not in each, '*' not in each]
            if all(rules): # substitution

                try:
                    original_aa, mut_aa, index = each[0], each[-1], int(each[1:-1])-1
                    # print (each,protein_id,index)
                    if index < len(protein_seq):
                        if protein_seq[index] != original_aa:
                            print (protein_id, index, 'not matching')
                    # assert protein_seq[index] == original_aa
                        else:
                            subs_aa_index.append((mut_aa, index))
                    # print (original_aa,mut_aa,index)
                except ValueError:
                    continue
    mod_prot_seq = protein_seq
    if subs_aa_index:
        for tp in subs_aa_index:
            mod_prot_seq = mod_prot_seq[:tp[1]]+tp[0]+mod_prot_seq[(tp[1]+1):]
    if del_aa_index:
        for tp in del_aa_index:
            mod_prot_seq = mod_prot_seq[:tp[1]]+mod_prot_seq[(tp[1]+1):]

    mut_protein_dict[protein_id]=mod_prot_seq
    counter +=1
    print(counter)

ppp.dump(mut_protein_dict,open('D:/data/pats/breast_cancer_raw/stomach_cancer_mutant_protein_dict_0826.p','wb'))

"""
mutant_protein_dict = ppp.load(open('D:/data/pats/breast_cancer_raw/breast_cancer_mutant_protein_dict_0826.p','rb'))
with open('D:/data/pats/breast_cancer_raw/breast_cancer_mutant_human.fasta','w',newline='\n') as f:
    for prot in mutant_protein_dict:
        prefix,descpt = description_dict[prot]
        f.write('>'+prefix+'|'+prot+'_mut'+'|'+descpt+'\n')
        seq = mutant_protein_dict[prot]
        block = range(0,len(seq)+60,60)
        for i in range(len(block)-1):
            f.write(seq[block[i]:block[i+1]]+'\n')
"""
