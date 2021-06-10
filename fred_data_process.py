from tsv_reader import combined_proteintsv_map, protein_info_from_combined
from pandas import ExcelWriter
import pandas as pd
from dash_dataframe_naba import dash_dataframe
import os
from protein_coverage import fasta_reader

fasta_path = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
protein_dict = fasta_reader(fasta_path)


base_path = 'D:/data/Naba_deep_matrisome/05142021_secondsearch/SNEDC/'
folders = [f for f in os.listdir(base_path) if '.' not in f]
psm_path_list = [base_path + each + '/psm.tsv' for each in folders]
pep_path_list = [base_path + each + '/peptide.tsv' for each in folders]
file_protein_cov_dict, id_pep_dict = dash_dataframe(pep_path_list,psm_path_list,protein_dict,psm_path_list,psm_path_list,psm_path_list)

file_path = 'D:/data/Naba_deep_matrisome/matrisome coverage.xlsx'
df = pd.read_excel(file_path)
df = df.drop_duplicates()
ecm_prot_list = df['protein_id'].tolist()

combined_protein_path = 'D:/data/Naba_deep_matrisome/05142021_secondsearch/SNEDC/combined_protein.tsv'
file_protein_spec_dict = combined_proteintsv_map(combined_protein_path)
protein_info_dict = protein_info_from_combined(combined_protein_path)

with ExcelWriter('D:/data/Naba_deep_matrisome/05142021_secondsearch/SNEDC/nonmat_protein_spec_SNEDC.xlsx') as writer:
    for file in file_protein_spec_dict:
        if file != 'Summarized':
            info_list = [[prot,file_protein_spec_dict[file][prot], file_protein_cov_dict[file][prot], protein_info_dict[prot][1], protein_info_dict[prot][2]]
                         for prot in file_protein_spec_dict[file] if prot not in ecm_prot_list and file_protein_cov_dict[file][prot]!=0]
            df = pd.DataFrame(info_list, columns=['protein','total_spec_count','coverage', 'gene_name', 'description'])
            df.to_excel(writer, '%s' % file)
