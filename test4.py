from tsv_reader import peptide_counting

psm_set = set()
counter = 0
file_path = 'C:/Users/gao lab computer/Downloads/163-3A_1_psm_fido.tsv'
with open(file_path,'r') as f:
    next(f)
    for line in f:
        line_split = line.split('\t')
        pep_seq = line_split[4]
        counter+=1
        print (pep_seq)
        psm_set.add(pep_seq)

print (counter, len(psm_set))