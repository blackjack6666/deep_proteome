import pandas as pd
from multiprocessing_naive_algorithym import extract_UNID_and_seq,creat_total_seq_line,creat_ID_pep_dict,\
    read_position_ID_into_dict, creat_pep_ID_dict, create_unique_id_peptide_dict
from tsv_reader import psm_reader,peptide_counting, protein_tsv_reader
from calculations_and_plot import identified_proteome_cov
import aho_corasick
import os
from protein_coverage import fasta_reader
import matplotlib.pyplot as plt
import statistics
import pickle as p
import seaborn as sns
from math import log10,log2


def combined_coverage(pep_path_list,protein_dict,ecm_protein_list):
    id_list, seq_list = extract_UNID_and_seq(protein_dict)
    seq_line = creat_total_seq_line(seq_list, sep='|')
    pos_id_dict = read_position_ID_into_dict(id_list, seq_list, seq_line)
    peptide_list = [pep for pep_path in pep_path_list for pep in peptide_counting(pep_path)]
    automaton = aho_corasick.automaton_trie(peptide_list)
    aho_result = aho_corasick.automaton_matching(automaton, seq_line)
    coverage_dict = identified_proteome_cov(aho_result, protein_dict)[1]
    ecm_coverage_dict = {p:coverage_dict[p] for p in coverage_dict if p in ecm_protein_list}

    return ecm_coverage_dict


def dash_dataframe(pep_path_list, psm_path_list, protein_dict, ecm_prot_list, ecm_info_dict, output_path):
    file_name_number_dict = {pep_tsv.split('/')[-2]: i for pep_tsv, i in zip(pep_path_list, range(len(pep_path_list)))}
    id_list, seq_list = extract_UNID_and_seq(protein_dict)
    seq_line = creat_total_seq_line(seq_list,sep="|")
    pos_id_dict = read_position_ID_into_dict(id_list, seq_list, seq_line)

    info_list = []
    file_id_peptide_dict = {}
    file_id_pep_count_dict = {}
    file_protein_cov_dict = {}
    file_protein_pep_spec_dict = {}
    file_protein_spec_dict = {}
    file_unique_id_pep_dict = {}
    file_unique_id_pep_count_dict = {}
    for pep_tsv, psm_tsv in zip(pep_path_list, psm_path_list):
        file_name = pep_tsv.split('/')[-2]
        # protein_list=protein_tsv_reader(pep_tsv.replace("peptide","protein"))
        print(file_name)
        pep_list = peptide_counting(pep_tsv)
        automaton = aho_corasick.automaton_trie(pep_list)
        aho_result = aho_corasick.automaton_matching(automaton, seq_line)
        coverage_dict = identified_proteome_cov(aho_result, protein_dict)[-1]
        file_protein_cov_dict[file_name] = coverage_dict
        id_pep_dict = creat_ID_pep_dict(aho_result, pos_id_dict)
        id_pep_count_dict = {id:len(id_pep_dict[id]) for id in id_pep_dict}
        file_id_pep_count_dict[file_name] = id_pep_count_dict
        # id_pep_dict = {k:id_pep_dict[k] for k in id_pep_dict if k in ecm_prot_list}
        pep_id_dict = creat_pep_ID_dict(aho_result,pos_id_dict)
        unique_id_pep_dict, unique_id_pep_count_dict = create_unique_id_peptide_dict(pep_id_dict)

        file_unique_id_pep_count_dict[file_name] = unique_id_pep_count_dict
        file_unique_id_pep_dict[file_name] = unique_id_pep_dict
        # pep_id_dict = {k:pep_id_dict[k] for k in pep_id_dict if pep_id_dict[k] in ecm_prot_list}
        file_id_peptide_dict[file_name] = id_pep_dict
        # p.dump(id_pep_dict,open(pep_tsv.replace('peptide.tsv','id_peptide_dict.p'),'wb'))
        psm_dict = psm_reader(psm_tsv)[0]
        protein_pep_spec_count = {prot_id:{pep:psm_dict[pep] for pep in id_pep_dict[prot_id]} for prot_id in id_pep_dict}
        file_protein_pep_spec_dict[file_name] = protein_pep_spec_count
        prot_spec_dict = {}
        for id in id_pep_dict:
            spec = 0
            for pep in id_pep_dict[id]:
                spec += psm_dict[pep]
            prot_spec_dict[id] = spec
        file_protein_spec_dict[file_name]=prot_spec_dict
    #     protein_id_ls = [k for k in id_pep_dict if k in ecm_prot_list]
    #     # protein_id_ls = [k for k in id_pep_dict]
    #     file_list = [[prot,
    #                   len(protein_dict[prot]),
    #                   coverage_dict[prot],
    #                   ecm_info_dict[prot][0],
    #                   prot_spec_dict[prot],
    #                   ecm_info_dict[prot][2],
    #                   file_name,
    #                   file_name_number_dict[file_name]] for prot in protein_id_ls]
    #     info_list+=file_list
    #
    # info_df = pd.DataFrame(info_list,
    #                        columns=['protein_id', 'length', 'coverage', 'gene_name', 'spec_count', 'ecm_class','file_name',
    #                                 'file_number'])
    # info_df.to_csv(output_path,index=False)
    return file_protein_cov_dict, file_id_peptide_dict,file_unique_id_pep_dict,file_protein_spec_dict,file_unique_id_pep_count_dict,file_id_pep_count_dict
    # return file_protein_pep_spec_dict


if __name__=='__main__':
    import plotly.express as px
    import random
    import numpy as np
    file_path = 'D:/data/Naba_deep_matrisome/matrisome coverage.xlsx'
    df = pd.read_excel(file_path)
    df = df.drop_duplicates()
    ecm_prot_list = df['protein_id'].tolist()
    df = df.set_index('protein_id')
    ecm_info_dict = {each:(df.loc[each,'gene_id'], df.loc[each, 'loc'], df.loc[each,'category']) for each in ecm_prot_list}
    #
    base_path = 'D:/data/Naba_deep_matrisome/02152021_1/'
    folders = [f for f in os.listdir(base_path) if '.' not in f]
    psm_path_list = [base_path + each + '/psm.tsv' for each in folders]
    pep_path_list = [base_path + each + '/peptide.tsv' for each in folders]
    peptide_163_05_20 = list(set([pep for pep_file in pep_path_list[8:] for pep in peptide_counting(pep_file)]))
    peptide_163_3A = peptide_counting(pep_path_list[7])
    print (pep_path_list[1:7])

    fasta_path = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
    protein_dict = fasta_reader(fasta_path)

    #
    dash_dataframe(pep_path_list[:7],psm_path_list[:7],protein_dict,ecm_prot_list,ecm_info_dict)

    # calculate coverage ratio for each protein between time-series and non-time series


    from tsv_reader import venn_diagram_gen
    file_path = 'D:/data/Naba_deep_matrisome/matrisome coverage.xlsx'

    df = pd.read_excel(file_path)
    df = df.drop_duplicates()
    ecm_prot_list = df['protein_id'].tolist()
    new_protein_dict = {prot:protein_dict[prot] for prot in ecm_prot_list if prot in protein_dict}

    df = df.set_index('protein_id')

    ecm_info_dict = {each: (df.loc[each, 'gene_id'], df.loc[each, 'loc'], df.loc[each, 'category']) for each in
                     ecm_prot_list}
    id_list, seq_list = extract_UNID_and_seq(new_protein_dict)
    seq_line = creat_total_seq_line(seq_list,sep='|')
    automaton = aho_corasick.automaton_trie(peptide_163_05_20)
    aho_result = aho_corasick.automaton_matching(automaton,seq_line)
    coverage_dict = identified_proteome_cov(aho_result, new_protein_dict)[1]

    automaton2 = aho_corasick.automaton_trie(peptide_163_3A)
    aho_result2 = aho_corasick.automaton_matching(automaton2,seq_line)
    coverage_dict2 = identified_proteome_cov(aho_result2,new_protein_dict)[1]

    new_cov_dict = {}
    for prot in coverage_dict:
        new_cov_dict[prot] = coverage_dict[prot]/coverage_dict2[prot] if prot in coverage_dict2 else 100

    print (np.mean([new_cov_dict[each] for each in new_cov_dict if new_cov_dict[each]!=100]))

    print (new_cov_dict)
    df_dash = pd.read_csv('D:/data/Naba_deep_matrisome/02152021_1/dash_info.csv')
    ecm_class = df_dash.ecm_class.unique()
    length_of_color_scheme = len(ecm_class)
    color = ["#" + ''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(length_of_color_scheme)]

    color_map = {each: c for each, c in zip(ecm_class, color)}
    dash_info_dict = {row['protein_id']:(row['length'],row['gene_name'],row['uniprot_num'],row['ecm_class']) for index, row in df_dash.iterrows()}

    info_list = [[each, new_cov_dict[each],log10(dash_info_dict[each][0]),dash_info_dict[each][1],dash_info_dict[each][2],dash_info_dict[each][3]]
                 for each in new_cov_dict if new_cov_dict[each] != 100]


    cov_ratio_df = pd.DataFrame(info_list, columns=['protein_id','coverage_ratio', 'log10length', 'gene','uniprot_num','ecm_class'])
    pd.set_option('display.max_columns', None)
    print (cov_ratio_df[cov_ratio_df['coverage_ratio']==1])
    print (cov_ratio_df.shape)
    # sizes_dict = {cov:int(cov*20) for cov in cov_ratio_df['coverage_ratio']}
    ax = sns.scatterplot(data=cov_ratio_df,x="uniprot_num",y="log10length",size="coverage_ratio",sizes=(20,200),hue='ecm_class',palette='deep',alpha=0.5)
    ax.set_xlabel('Uniprot num', fontsize=15)
    ax.set_ylabel('log10(length)',fontsize=15)
    ax.legend(framealpha=0.5)
    # plt.show()

    # fig = px.scatter(cov_ratio_df, x="uniprot_num", y="length",
    #                  size="coverage_ratio", color='ecm_class', hover_name="gene",
    #                  log_x=False, size_max=55, color_discrete_map=color_map)
    # fig.show()

    # coverage_array = [v for v in coverage_dict.values()]
    # coverage_array2 = [v for v in coverage_dict2.values()]
    # print (statistics.mean(coverage_array))
    # plt.hist(coverage_array,bins=15, color='black',alpha=0.5,rwidth=1)
    # plt.hist(coverage_array2,bins=15,color='grey',alpha=0.2,rwidth=1.5)
    # plt.xlim(0,100)
    # plt.xlabel('sequence coverage%')
    # plt.ylabel('frequency')
    # # plt.title('ECM sequence coverage distribution in time lapsed and standard digestion in 18_2 sample')
    # plt.show()



    df = pd.read_csv('D:/data/Naba_deep_matrisome/02152021_1/dash_info.csv')
    protein_18_2A_list = df[df['file_number']==7]['protein_id'].tolist()
    protein_1805_1820_list = list(set([prot for i in range(8,14) for prot in df[df['file_number']==i]['protein_id'].tolist()]))

    venn_dict = {'18_2_time_lapsed_digestion': protein_1805_1820_list,'18_2_20hour_digestion': protein_18_2A_list}
    venn_diagram_gen(venn_dict)
    # print ([ecm_info_dict[prot] for prot in protein_1805_1820_list if prot not in protein_18_2A_list])
