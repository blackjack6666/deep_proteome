
from tsv_reader import combined_proteintsv_map, protein_info_from_combined, venn_diagram_gen,protein_info_from_fasta
import pandas as pd
from pandas import ExcelWriter
from multiprocessing_naive_algorithym import extract_UNID_and_seq,creat_total_seq_line,creat_ID_pep_dict,\
    read_position_ID_into_dict, creat_pep_ID_dict
import aho_corasick
from protein_coverage import fasta_reader
import os
import pickle as p

# protein_info_dict = protein_info_from_fasta('D:/data/proteome_fasta/uniprot-proteome_UP000005640.fasta')

# combined_prot_tsv = 'G:/XS/Burdette/4_8_21_search_result/combined_protein.tsv'
# info_dict = combined_proteintsv_map(combined_prot_tsv)
# protein_info_dict = protein_info_from_combined(combined_prot_tsv)
#
# with ExcelWriter('G:/XS/Burdette/4_8_21_search_result/protein_spec_count.xlsx') as writer:
#         for each in info_dict:
#             info_list = [[prot,info_dict[each][prot],protein_info_dict[prot][0],protein_info_dict[prot][1],protein_info_dict[prot][2]]
#                          for prot in info_dict[each]]
#             df = pd.DataFrame(info_list, columns=['protein id', 'spectra count', 'entry name', 'gene name', 'description'])
#             df.to_excel(writer,'%s' % each)


combined_prot_tsv = 'G:/XS/Burdette/4_8_21_search_result/combined_protein.tsv'
info_dict = combined_proteintsv_map(combined_prot_tsv)
# protein_info_dict = protein_info_from_combined(combined_prot_tsv)
print ([each for each in info_dict])
DMSO_47D_01, DMSO_47D_02 = [k for k in info_dict['47D_DMSO_01']], [k for k in info_dict['47D_DMSO_02']]
dmso_bai_overlap = set(DMSO_47D_01)&set(DMSO_47D_02)

bai_47d_01, bai_47d_02 = [k for k in info_dict['47D_P4_Bai_01']], [k for k in info_dict['47D_P4_Bai_02']]
bai_47d_overlap = set(bai_47d_01)&set(bai_47d_02)
bai_47d_unique = [each for each in bai_47d_overlap if each not in dmso_bai_overlap]
dmso_unique = [each for each in dmso_bai_overlap if each not in bai_47d_overlap]
for each in bai_47d_unique:
    if not each.startswith('CT'):
        print (each)
# print ([(each,int(info_dict['A12_P4_Bai_01'][each]+info_dict['A12_P4_Bai_02'][each]/2)) for each in bai_47d_unique if not each.startswith('CT')])

# venn_diagram_gen({'A12_DMSO':dmso_bai_overlap,'A12_P4_Bai':bai_47d_overlap})
# file_path = 'D:/data/Naba_deep_matrisome/matrisome coverage.xlsx'
# df = pd.read_excel(file_path)
# df = df.drop_duplicates()
# ecm_prot_list = df['protein_id'].tolist()
# df = df.set_index('protein_id')
# ecm_info_dict = {each:(df.loc[each,'gene_id'], df.loc[each, 'loc'], df.loc[each,'category']) for each in ecm_prot_list}
#
#
# base_path = 'D:/data/Naba_deep_matrisome/03122021/KOC/'
# folders = [f for f in os.listdir(base_path) if '.' not in f]
#
# id_pep_path_list = [base_path + each + '/id_peptide_dict.p' for each in folders]
# print (id_pep_path_list)
#
# with ExcelWriter('D:/data/Naba_deep_matrisome/03122021/KOC/KOC_id_pep.xlsx') as writer:
#
#         for each in id_pep_path_list:
#             file = '_'.join(each.split('/')[-3:-1])
#             id_pep_dict = p.load(open(each,'rb'))
#
#             info_list = [[id, ecm_info_dict[id][0], ecm_info_dict[id][-1], id_pep_dict[id]]for id in id_pep_dict]
#             df = pd.DataFrame(info_list, columns=['protein id', 'gene', 'category', 'peptides'])
#             df.to_excel(writer,'%s' % file)