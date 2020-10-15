from collections import defaultdict


def peptide_info_reader(psm_path:str):
    """
    read peptide sequence as key, [(filename, spectram number,mass, charge,retention time),()] as key
    :param psm_path:
    :return:
    """
    info_dict = defaultdict(list)
    with open(psm_path,'r') as f:
        for i in range(1):
            next(f)
        for line in f:
            line_split = line.split('\t')
            pep_seq = line_split[1]
            file_name = line_split[0].split('.')[0]
            spec_no = int(line_split[0].split('.')[-2])
            mass = float(line_split[6])
            charge = int(line_split[4])
            ret_time = float(line_split[5])
            info_dict[pep_seq].append([file_name,spec_no,mass,charge,ret_time])
    return info_dict


def peptide_charger_reader(pep_tsv):
    """
    get charge of the peptide
    :param pep_tsv:
    :return:
    """
    peptide_charge_dict = {}
    with open(pep_tsv,'r') as file_open:
        next(file_open)
        for line in file_open:
            line_split = line.split('\t')
            peptide_seq = line_split[0]
            charge = [line_split[2]] if ', ' not in line_split[2] else line_split[2].split(', ')
            peptide_charge_dict[peptide_seq] = charge
        return peptide_charge_dict


def protein_tsv_reader(protein_tsv_file):

    with open(protein_tsv_file, 'r') as file_open:
        next(file_open)
        return [line.split("\t")[3] for line in file_open]


def protein_tsv_reader_no_contam(protein_tsv_file):
    with open(protein_tsv_file, 'r') as file_open:
        next(file_open)
        return [line.split("\t")[3] for line in file_open if 'contaminant' not in line.split("\t")[3]]

def peptide_counting(peptide_tsv_file):

    with open(peptide_tsv_file, 'r') as file_open:
        next(file_open)

        peptide_list = [line.split("\t")[0] for line in file_open]
    return peptide_list

def psm_reader(psm_path,fragpipe_ver=13.0):
    pep_spec_count_dict = defaultdict(int)
    ret_pep_dict = {}
    with open(psm_path, 'r') as f:
        for i in range(1):
            next(f)
        for line in f:
            line_split = line.split('\t')
            pep_seq = line_split[2] if fragpipe_ver==13.0 else line_split[1]
            retention_time = float(line_split[5])/60  if fragpipe_ver==13.0 else float(line_split[4])/60 # in minute
            pep_spec_count_dict[pep_seq]+=1
            ret_pep_dict[retention_time] = pep_seq
    return pep_spec_count_dict, ret_pep_dict

def venn_diagram_gen(dictionary, title=''): # parameter could be a dictionary of proteins or peptides from different samples {'sample1': [], 'sample2': []}
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2, venn3

    value_list_of_sets = [set(l) for l in dictionary.values()]
    sample_name_list = [n for n in dictionary.keys()]
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))

    if len(dictionary) == 2:  # two samples for venn diagram
        venn2(value_list_of_sets, set_labels=sample_name_list)

    elif len(dictionary) == 3:  # 3 samples for venn diagram
        venn3(value_list_of_sets, set_labels=sample_name_list)

    else:
        print ('Error: only 2 or 3 comparison for venn diagram are accepted in this script.')

    plt.title(title)
    plt.show()


def venn_diagram_gen2(dictionary, title=None): # parameter could be a dictionary of proteins or peptides from different samples {'sample1': [], 'sample2': []}
    import matplotlib.pyplot as plt
    import matplotlib
    import venn

    if len(dictionary)>=6:
        raise ValueError('maximum groups to compare is 5')

    else:
        matplotlib.rcParams.update({'font.size': 15})

        venn_dict = {each:set(dictionary[each]) for each in dictionary}
        ax = venn.venn(venn_dict)
        if title:
            ax.set_title(title)
    plt.show()

if __name__=="__main__":
    from glob import glob
    import numpy as np
    path = 'D:/data/deep_proteome/20200915_ct_50C*/protein.tsv'
    file_list = glob(path)
    print (file_list)

    venn_dict = {'_'.join(each_file.split('\\')[-2].split('_')[1:]):protein_tsv_reader(each_file)
                 for each_file in file_list[1:6]}
    venn_diagram_gen2(venn_dict)

    # from calculations_and_plot import miss_cleavage_identify
    # for each_file in file_list:
    #     print (each_file)
    #     pep_list = peptide_counting(each_file)
    #     miss_cleav_dict = miss_cleavage_identify(pep_list,regex_pattern=r'(?:F|W|Y)\w+')
    #     print (float(np.count_nonzero([each for each in miss_cleav_dict.values()]))/len(pep_list))

