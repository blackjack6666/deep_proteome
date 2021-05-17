from tsv_reader import peptide_counting

pep_path1 = 'D:/data/pats/results/hek_zz_sp_isoform_no_razor/peptide.tsv'
pep_path2 = 'D:/data/pats/results/hek_zz_sp_isoform_default_closed_param/peptide.tsv'

pep_list1, pep_list2 = peptide_counting(pep_path1), peptide_counting(pep_path2)

print ([pep for pep in pep_list2 if pep not in pep_list1])
print (len([pep for pep in pep_list2 if pep not in pep_list1]))