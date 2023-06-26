def dta_charge_reader(file_location_path):
    """
    -----
    read peptide seq, PSM, protein ID from dta files
    -----
    :param file_location_path: could be str folder path or dta file path or list of dta files or folder path
    :return: peptide list, PSM dictionary showing number of PSM, protein ID list
    """
    from glob import glob
    from collections import defaultdict

    dta_files = []

    # single dta read
    if isinstance(file_location_path,str) and file_location_path.endswith('.dta'):
        dta_files = [file_location_path]
        print ('reading single dta file')

    # dta folder read
    elif isinstance(file_location_path,str) and not file_location_path.endswith('.dta'):
        dta_files = glob(file_location_path + '*.dta')
        print ('reading single folder')

    # multiple dta files read or multiple dta folder read
    elif isinstance(file_location_path,list):
        for each_dta in file_location_path:
            # when parameter is dta file

            if each_dta.endswith('.dta'):
               dta_files.append(each_dta)

            # when parameter is a folder path
            else:
                dta_files += glob(each_dta+'*.dta')
    else:
        raise ValueError('parameter should be string folder path or dta file path or list of dta files or folder paths')

    # exclude wash and hela files
    clean_dta_files = []
    for each in dta_files:
        wash_hela = 0
        for word in ['wash', 'Wash', 'WASH', 'Hela', 'hela', 'HELA']:
            if word in each:
                wash_hela += 1
                break
        if wash_hela == 0:
            clean_dta_files.append(each)

    print (clean_dta_files)

    # read info.

    seq_charge_dict = defaultdict(list)
    for dta_file in clean_dta_files:
        print ('reading: ', dta_file)
        with open(dta_file, 'r') as file_open:
            for i in range(29):
                next(file_open)
            Reverse_start = 0
            for line in file_open:
                line_split = line.split('\t')
                if line.startswith('Reverse_') or line.startswith('Rev_'):
                    Reverse_start = 1
                elif line.startswith('sp') or line.startswith('tr'):
                    Reverse_start = 0

                elif len(line_split) == 15 and Reverse_start == 0:
                    pep_seq = line_split[-1].split('.')[1]
                    charge = int(line_split[1].split('.')[-1])
                    seq_charge_dict[pep_seq].append(charge)

    seq_charge_dict = {each:list(set(seq_charge_dict[each])) for each in seq_charge_dict}
    print ('reading done.')
    return seq_charge_dict


if __name__=='__main__':
    from calculations_and_plot import identified_proteome_cov
    import aho_corasick, multiprocessing_naive_algorithym
    from glob import glob
    from protein_coverage import fasta_reader2
    import matplotlib.pyplot as plt
    import pickle as ppp
    import seaborn as sns
    import numpy as np
    import pandas as pd

    # fasta_path = 'D:/data/proteome_fasta/uniprot-proteome_UP000005640.fasta'
    # protein_dict = fasta_reader2(fasta_path)
    #
    # ID_list, seq_list = multiprocessing_naive_algorithym.extract_UNID_and_seq(protein_dict)
    # seq_line = multiprocessing_naive_algorithym.creat_total_seq_line(seq_list)
    # path = 'D:/data/dta_files/dta_result/'
    # peptide_list = [pep for pep in dta_charge_reader(path)]
    # automaton = aho_corasick.automaton_trie(peptide_list)
    # aho_result = aho_corasick.automaton_matching(automaton, seq_line)
    # prot_cov_dict = identified_proteome_cov(aho_result,protein_dict)[1]
    # cov_list = [v for v in prot_cov_dict.values()]

    # ppp.dump(prot_cov_dict, open('bioplex_protein_coverage_dict.p','wb'))
    # prot_cov_dict = ppp.load(open('bioplex_protein_coverage_dict.p', 'rb'))
    # cov_array = np.array([v for v in prot_cov_dict.values()])
    # sort = -np.sort(-cov_list)
    # reshape = sort_array.reshape((250,334))

    # matrisome_cov_csv = 'D:/data/Naba_deep_matrisome/matrisome coverage_norepeat.xlsx'
    # df_mat = pd.read_excel(matrisome_cov_csv, index_col=0)
    # cov_array = df_mat['cov'].to_numpy()*100
    matrisome_cov_csv = 'F:/matrisomedb2.0/statistics/glob_seq_coverage_1.tsv'
    df_mat = pd.read_csv(matrisome_cov_csv,sep='\t',index_col=0)
    cov_array = df_mat['Sequence coverage'].to_numpy()

    sort = -np.sort(-cov_array)
    fig, ax = plt.subplots(figsize=(6.5,5))
    ax.plot(range(len(sort)),sort, '-', color='#FFFFFF', alpha=0.3)
    # ax.fill_between(range(len(sort)),sort, color='#fcba03',alpha=0.3)
    # ax.fill_between(range(len(sort)), sort, [100]*len(sort), color='#e6e5e3')
    ax.set_xlim(0,len(sort))
    ax.set_ylim(0,100)
    ax.set_xticks([])
    ax.tick_params(axis='both', which='major', labelsize=15)

    plt.show()

    # reshape = sort.reshape((49,34))
    # pd.DataFrame(reshape).to_excel('D:/data/bioplex/matrisomedb_coverage.xlsx')
    # sns.heatmap(data=reshape,cmap="viridis")
    # plt.show()
    # fig, ax = plt.subplots(1,1)
    # sns.distplot(cov_list, color='black', ax=ax)
    # plt.xlabel('sequence coverage')
    # plt.ylabel('density')
    # plt.show()

