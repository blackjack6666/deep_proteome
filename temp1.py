from tsv_reader import psm_reader,venn_diagram_gen, peptide_counting, venn_diagram_gen2

trypsin_37C_15min = 'D:/data/deep_proteome/20200806_tryp_RT_25min/psm.tsv'
trypsin_37C_2h = 'D:/data/deep_proteome/20200806_tryp_RT_40min/psm.tsv'
trypsin_37C_19h = 'D:/data/deep_proteome/20200806_tryp_37C_19H/psm.tsv'
trypsin_RT_19h = 'D:/data/deep_proteome/20200806_tryp_RT_19H/psm.tsv'

trypsin_37C_15min_peptsv = 'D:/data/deep_proteome/20200806_tryp_RT_25min/peptide.tsv'
trypsin_37C_2h_peptsv = 'D:/data/deep_proteome/20200806_tryp_RT_40min/peptide.tsv'
trypsin_37C_19h_peptsv = 'D:/data/deep_proteome/20200806_tryp_RT_19H/peptide.tsv'

trypsin_37C_15min_psm_dict = psm_reader(trypsin_37C_15min)[0]
trypsin_37C_2h_psm_dict = psm_reader(trypsin_37C_2h)[0]
trypsin_37C_19h_psm_dict = psm_reader(trypsin_37C_19h)[0]

trypsin_37C_15min_peplist = peptide_counting(trypsin_37C_15min_peptsv)
trypsin_37C_2h_peplist = peptide_counting(trypsin_37C_2h_peptsv)
trypsin_37C_19h_peplist = peptide_counting(trypsin_37C_19h_peptsv)

trypsin_37C_15min_psm_list = [pep+'_'+str(i) for pep in trypsin_37C_15min_psm_dict for i in range(trypsin_37C_15min_psm_dict[pep])]
trypsin_37C_2h_psm_list = [pep+'_'+str(i) for pep in trypsin_37C_2h_psm_dict for i in range(trypsin_37C_2h_psm_dict[pep])]
trypsin_37C_19h_psm_list = [pep+'_'+str(i) for pep in trypsin_37C_19h_psm_dict for i in range(trypsin_37C_19h_psm_dict[pep])]

venn_dict = {'tryp_RT_25min':trypsin_37C_15min_peplist,
             'tryp_RT_40min': trypsin_37C_2h_peplist,
             'tryp_RT_19h':trypsin_37C_19h_peplist}

venn_diagram_gen(venn_dict, title='pep comparison')