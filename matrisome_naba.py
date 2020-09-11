import pandas as pd
from protein_coverage import fasta_reader, read_description

file_path = 'D:/data/Naba_deep_matrisome/matrisome coverage.xlsx'

df = pd.read_excel(file_path)

ecm_protein_id_set = set(df['protein_id'].to_list())

fasta_path = 'D:/data/proteome_fasta/uniprot-proteome_UP000000589_mouse.fasta'
protein_dict = fasta_reader(fasta_path)
protein_descript_dict = read_description(fasta_path)

count=0
with open('D:/data/Naba_deep_matrisome/mouse_ecm_costom_proteome_db.fasta', 'w', newline='\n') as file_w:
    for each in ecm_protein_id_set:
        if each in protein_dict:
           count+=1
           seq = protein_dict[each]
           prefix,description = protein_descript_dict[each]
           file_w.write('>'+str(prefix)+'|'+str(each)+'|'+str(description)+'\n')
           bins = range(0,len(seq)+60,60)
           for i in range(len(bins)-1):
               file_w.write(seq[bins[i]:bins[i+1]]+'\n')

print (count,'fasta done')


