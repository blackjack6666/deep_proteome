"""data processing and analysis of native limited proteolysis, DIA-NN output, written by blackjack, 01/27/2023
"""
import pandas as pd
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


if __name__ == '__main__':
    # get_unique_peptides()
    get_protein_set()