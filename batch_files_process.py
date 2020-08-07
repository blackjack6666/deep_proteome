"""
read psm and charge for each file and combine them together into a dictionary
"""
def psm_reader(psm_path):
    from collections import defaultdict
    pep_spec_count_dict = defaultdict(int)
    pep_charge_dict = defaultdict(list)
    with open(psm_path, 'r') as f:
        for i in range(1):
            next(f)
        for line in f:
            line_split = line.split('\t')
            pep_seq = line_split[1]
            charge = int(line_split[0].split('.')[-1])
            pep_spec_count_dict[pep_seq] += 1
            pep_charge_dict[pep_seq].append(charge)
    return pep_spec_count_dict, pep_charge_dict


if __name__=='__main__':
    import glob
    import pickle as ppp
    path = 'D:/data/deep_proteome/'
    psm_files = glob.glob(path+'*.tsv')
    print (psm_files)
    sum_dict = {}
    file_dict = {}
    folder_name = psm_files[0].split('\\')[-1].split('_')[0]
    for each_file in psm_files:

        psm_file_name = each_file.split('\\')[-1]
        psm_file_name_split = psm_file_name.split('_')

        enzyme_time = '_'.join(psm_file_name_split[1:3])
        pep_spec_count_dict,pep_charge_dict = psm_reader(each_file)
        print (each_file)
        file_dict[enzyme_time]=(pep_spec_count_dict,pep_charge_dict)
    sum_dict[folder_name]=file_dict

    print(sum_dict['20200706']['tryp_5h'][1])