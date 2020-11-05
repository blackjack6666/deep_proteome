# from tsv_reader import psm_reader,venn_diagram_gen, peptide_counting, venn_diagram_gen2
#
# trypsin_37C_15min = 'D:/data/deep_proteome/20200806_tryp_RT_25min/psm.tsv'
# trypsin_37C_2h = 'D:/data/deep_proteome/20200806_tryp_RT_40min/psm.tsv'
# trypsin_37C_19h = 'D:/data/deep_proteome/20200806_tryp_37C_19H/psm.tsv'
# trypsin_RT_19h = 'D:/data/deep_proteome/20200806_tryp_RT_19H/psm.tsv'
#
# trypsin_37C_15min_peptsv = 'D:/data/deep_proteome/20200806_tryp_RT_25min/peptide.tsv'
# trypsin_37C_2h_peptsv = 'D:/data/deep_proteome/20200806_tryp_RT_40min/peptide.tsv'
# trypsin_37C_19h_peptsv = 'D:/data/deep_proteome/20200806_tryp_RT_19H/peptide.tsv'
#
# trypsin_37C_15min_psm_dict = psm_reader(trypsin_37C_15min)[0]
# trypsin_37C_2h_psm_dict = psm_reader(trypsin_37C_2h)[0]
# trypsin_37C_19h_psm_dict = psm_reader(trypsin_37C_19h)[0]
#
# trypsin_37C_15min_peplist = peptide_counting(trypsin_37C_15min_peptsv)
# trypsin_37C_2h_peplist = peptide_counting(trypsin_37C_2h_peptsv)
# trypsin_37C_19h_peplist = peptide_counting(trypsin_37C_19h_peptsv)
#
# trypsin_37C_15min_psm_list = [pep+'_'+str(i) for pep in trypsin_37C_15min_psm_dict for i in range(trypsin_37C_15min_psm_dict[pep])]
# trypsin_37C_2h_psm_list = [pep+'_'+str(i) for pep in trypsin_37C_2h_psm_dict for i in range(trypsin_37C_2h_psm_dict[pep])]
# trypsin_37C_19h_psm_list = [pep+'_'+str(i) for pep in trypsin_37C_19h_psm_dict for i in range(trypsin_37C_19h_psm_dict[pep])]
#
# venn_dict = {'tryp_RT_25min':trypsin_37C_15min_peplist,
#              'tryp_RT_40min': trypsin_37C_2h_peplist,
#              'tryp_RT_19h':trypsin_37C_19h_peplist}
#
# venn_diagram_gen(venn_dict, title='pep comparison')

from data_preprocess2 import *
import pickle as ppp
from tsv_reader import peptide_counting,psm_reader
import numpy as np

fasta_path = 'D:/data/proteome_fasta/uniprot-proteome_UP000000589_mouse.fasta'
proteome_dict = fasta_reader2(fasta_path)
protein_list = ['P62918']
protein_miss_clea_loc_dict = cleavage_site_identify(protein_list, proteome_dict, 'trypsin')

# therotical
polymers_dict = polymer_miss_cleav_charc(protein_miss_clea_loc_dict, proteome_dict)
polymers_list = [polymers_dict[each][each_i] for each in polymers_dict for each_i in polymers_dict[each]]

matrix = []
for each_polymer in polymers_list:
    one_d_arry = []
    for each_aa in each_polymer:
        one_d_arry.append(each_aa)
    matrix.append(one_d_arry)
matrix = np.array(matrix)
print (matrix.shape)
ppp.dump(matrix, open('P62918_matrix_2d_array.p', 'wb'))


# peptide_tsv_path = "D:/data/Naba_deep_matrisome/10_30/search_result/18_2_dec/peptide.tsv"
# psm_tsv_path = "D:/data/Naba_deep_matrisome/10_30/search_result/18_2_dec/psm.tsv"
# pep_list = peptide_counting(peptide_tsv_path)
# id_pep_dict,protein_list = protein_id_peplist_dict_getter(proteome_dict, pep_list)
# psm_dict = psm_reader(psm_tsv_path)[0]
# protein_polymer_sc_dict = map_to_cleavage(id_pep_dict,protein_miss_clea_loc_dict,polymers_dict,proteome_dict,psm_dict)
# cleavage_site_label_dict, protein_poly_dict, uncertain_polymer_no = cleavage_site_label(protein_polymer_sc_dict, protein_polymer_convert(polymers_dict))
# print (cleavage_site_label_dict)
# ppp.dump(cleavage_site_label_dict, open('P62918_polymer_dict.p', 'wb'))