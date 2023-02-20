"""data processing and analysis of native limited proteolysis, DIA-NN output, written by blackjack, 01/27/2023
"""
import pandas as pd
import numpy as np
import pickle as pk
from collections import defaultdict


def peptide_getter():
    """
    output peptide list for each time point in a dictionary
    :return:
    """
    file_peptide_dict = defaultdict(set)
    df = pd.read_csv('F:/native_digestion/01242023/time_points/report.tsv',sep='\t')
    file_list, psm_list = df['File.Name'], df['Stripped.Sequence']
    for f, psm in zip(file_list,psm_list):
        file_peptide_dict[f.split('\\')[-1].split('_rerun.mzML')[0]].add(psm)
    print ([k for k in file_peptide_dict])
    return pk.dump(file_peptide_dict, open('F:/native_digestion/01242023/time_points/f_peptides_dict.p','wb'))


def get_unique_peptides():
    """
    get unique peptides at each time point
    :return:
    """
    unique_peptides_dict = {}
    peptide_list = []
    f_list = ['tryps_0005min','tryps_0010min','tryps_0015min','tryps_0020min','tryps_0030min','tryps_0040min',
              'tryps_0050min','tryps_0060min','tryps_0120min','tryps_0180min','tryps_0240min','tryps_1440min',
              'tryps_leftover'] # order by time
    f_peptide_dict = pk.load(open('F:/native_digestion/01242023/time_points/f_peptides_dict.p','rb'))
    for f in f_list:
        unique_pep_list = [each for each in f_peptide_dict[f] if each not in peptide_list]
        peptide_list += unique_pep_list
        unique_peptides_dict[f] = unique_pep_list
    print ([(each, len(unique_peptides_dict[each])) for each in unique_peptides_dict])
    # pk.dump(unique_peptides_dict,open('F:/native_digestion/01242023/time_points/f_unique_peptides_dict.p','wb'))


def get_protein_set():
    """
    get total protein list
    :return:
    """
    protein_set = set()
    df = pd.read_csv('F:/native_digestion/01242023/time_points/report.pg_matrix.tsv',sep='\t')
    protein_list = df['Protein.Ids']
    for prot_str in protein_list:
        if ';' in prot_str:
            for prot in prot_str.split(';'):
                protein_set.add(prot)
        else:
            protein_set.add(prot_str)
    print (len(protein_set), list(protein_set)[:10])
    # pk.dump(protein_set,open('F:/native_digestion/01242023/time_points/proteinid_set.p','wb'))


def protein_intensity():
    """
    get peptide intensity vector from pg report
    :return:
    """
    protein_int_dict = {}
    df = pd.read_csv('F:/native_digestion/01242023/time_points/report.pg_matrix.tsv', sep='\t')
    for row in df.itertuples(index=False):
        protein_id_list = [row[1]] if ';' not in row[1] else row[1].split(';')
        int_array = np.array(row[-13:])
        for each in protein_id_list:
            protein_int_dict[each] = int_array
    pk.dump(protein_int_dict,open('F:/native_digestion/01242023/time_points/prot_intensity_dict.p','wb'))


def combine_distance_intensity():
    """
    combine cleavage to center distance and protein intensity
    :return:
    """
    protein_int_dict = pk.load(open('F:/native_digestion/01242023/time_points/prot_intensity_dict.p','rb'))
    distance_df = pd.read_excel('F:/native_digestion/01242023/analysis/distance_to_center_all.xlsx',index_col=0)
    distance_df = distance_df.fillna(0)
    new_df = pd.DataFrame(index=distance_df.index,columns=distance_df.columns)
    for row in distance_df.itertuples(index=True):
        prot = row[0]
        int_array = protein_int_dict[prot]
        int_array_fillna = np.where(np.isnan(int_array),0,int_array)
        normalized_int = (int_array_fillna-np.min(int_array_fillna))/(np.max(int_array_fillna)-np.min(int_array_fillna))
        multiply = np.multiply(np.array(row[1:]), normalized_int)
        new_df.loc[prot,:] = multiply
    new_df.to_csv('F:/native_digestion/01242023/analysis/distance_to_center_times_normIntensity.tsv',sep='\t')


def filter_df():
    # filter out all zeros in df
    import time
    time.sleep(3)
    df = pd.read_csv('F:/native_digestion/01242023/analysis/distance_to_center_times_normIntensity.tsv',sep='\t',index_col=0)
    df = df.copy()
    data = []
    index = []
    for row in df.itertuples():
        # filter rows with all 0s
        if np.count_nonzero([row[i]==0 for i in range(-13,0)])>=6:  # if there are more than 6 zeros
        # if all([row[i] == "0" for i in range(1, len(row))]):
            continue
        else:
            data.append([i for i in row][1:])
            index.append(row[0])
    new_df = pd.DataFrame(data, columns=df.columns,index=index)
    new_df.to_csv('F:/native_digestion/01242023/analysis/distance_to_center_times_normIntensity_filter.tsv',sep='\t')


def cluster_map():
    import matplotlib.pyplot as plt
    import seaborn as sns

    data = pd.read_csv('F:/native_digestion/01242023/analysis/distance_to_center_times_normIntensity_filter_norm01.tsv',sep='\t',index_col=0)
    # new_data = pd.DataFrame(index=data.index, columns=data.columns)
    # normalize each protein between 0 and 1
    # for prot, array in zip(data.index, data.to_numpy()):
    #     normalize_0_1 = (array-np.min(array))/(np.max(array)-np.min(array))
    #     new_data.loc[prot,:] = normalize_0_1
    # new_data.to_csv('F:/native_digestion/01242023/analysis/distance_to_center_times_normIntensity_filter_norm01.tsv', sep='\t')
    fig, ax = plt.subplots(1, 1, figsize=(8, 15))
    g = sns.clustermap(data=data,col_cluster=False,cmap="YlGnBu",yticklabels=False)
    plt.show()

if __name__ == '__main__':
    # get_unique_peptides()
    # get_protein_set()
    # print (protein_intensity())
    # combine_distance_intensity()
    # filter_df()
    cluster_map()