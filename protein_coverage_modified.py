from protein_coverage import read_fasta_into_dict
from collections import OrderedDict
import glob
import numpy as np
import time
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns


time_start = time.clock()
fasta_name = 'C:/uic/lab/data/xinhao_data1/uniprot-proteome_UP000005640.fasta'
aa_frequency_dict = read_fasta_into_dict(fasta_name)[1]
fasta_protein_seq_dict = read_fasta_into_dict(fasta_name)[0]
aa_frequency_identified_dict = {}


def read_peptide_into_dict(dta_file):
    pep_dict = OrderedDict()
    ID = ''
    pep_seq_list = []
    with open(dta_file, 'r') as fileopen:
        for i in range(29):
            next(fileopen)
        Reverse_start = 0
        for line in fileopen:
            if line.startswith('Reverse_'):
                Reverse_start = 1
            elif line.startswith('sp'):
                Reverse_start = 0
                pep_dict[ID] = list(set(pep_seq_list))
                pep_seq_list = []
                ID = line.split('\t')[0].split('|')[1]
            elif len(line.split('\t')) == 15 and Reverse_start == 0:
                peptide = line.split('\t')[-1].split('.')[1]
                pep_seq_list.append(peptide)
    if ID not in pep_dict:
        pep_dict[ID] = list(set(pep_seq_list))
    del pep_dict['']
    for i in np.arange(len(pep_dict)):
        if pep_dict.values()[i] == []:
            for j in np.arange(len(pep_dict) - i):
                if pep_dict.values()[i + j] == []:
                    continue
                else:
                    pep_dict.values()[i] += pep_dict.values()[i + j]
                    break
    pep_dict = dict(pep_dict)
    return pep_dict


def creat_aa_frequency_dict(dta_file):
    pep_dict = read_peptide_into_dict(dta_file)
    for ID in pep_dict.keys(): #match peptide ID with protein ID in fasta file
        if ID in fasta_protein_seq_dict and ID not in aa_frequency_identified_dict:
            protein_seq1 = fasta_protein_seq_dict[ID]
            aa_frequency_identified_dict[ID] = [0] * len(protein_seq1)
            for pep in pep_dict[ID]:
                startpos = protein_seq1.find(pep)
                endpos = startpos + len(pep)
                aa_frequency_dict[ID][startpos:endpos] = [x+1 for x in aa_frequency_dict[ID][startpos:endpos]]
                aa_frequency_identified_dict[ID][startpos:endpos] = [x + 1 for x in
                                                                     aa_frequency_identified_dict[ID][startpos:endpos]]
        elif ID in fasta_protein_seq_dict and ID in aa_frequency_identified_dict:
            protein_seq2 = fasta_protein_seq_dict[ID] # protein_seq1 = protein_seq2
            for pep in pep_dict[ID]:
                startpos = protein_seq2.find(pep)
                endpos = startpos + len(pep)
                aa_frequency_dict[ID][startpos:endpos] = [x+1 for x in aa_frequency_dict[ID][startpos:endpos]]
                aa_frequency_identified_dict[ID][startpos:endpos] = [x + 1 for x in
                                                                     aa_frequency_identified_dict[ID][startpos:endpos]]
    return aa_frequency_identified_dict, aa_frequency_dict


def coverage():
    identified_proteins_ID_list = []
    average_precentage = []
    identified_zero_count = 0
    identified_total_len = 0
    total_all_protein_len = 0
    path = 'C:/uic/lab/data/xinhao_data1/'
    dtafiles = glob.glob(path+'*.dta')
    for dtafile in dtafiles:
        Total_identified_protein_freq = creat_aa_frequency_dict(dtafile)[0]
        Total_proteins_freq = creat_aa_frequency_dict(dtafile)[1]
    for ID in Total_identified_protein_freq:
        identified_proteins_ID_list.append(ID)
        zero_count = Total_identified_protein_freq[ID].count(0)
        identified_zero_count += zero_count
        total_len = len(Total_identified_protein_freq[ID])
        identified_total_len += total_len
        percentage = float(total_len-zero_count)/total_len*100
        average_precentage.append(percentage)

    for ID in Total_proteins_freq:
        total_all_protein_len += len((Total_proteins_freq[ID]))
    total_average = np.mean(average_precentage)
    overall_coverage_for_identified_protein = float(identified_total_len-identified_zero_count)/identified_total_len*100
    overall_coverage_for_all_proteins = float(identified_total_len-identified_zero_count)/total_all_protein_len*100
    return total_average, overall_coverage_for_identified_protein, overall_coverage_for_all_proteins, average_precentage, identified_proteins_ID_list

#print "The average identified protein coverage is %.2f%%\nThe overall coverage for all "
if __name__ == '__main__':
    percentage_list, ID_list = coverage()[3:5]
    df = pd.DataFrame(dict(identified_protein_ID=ID_list, coverage_percentage=percentage_list))
    df.to_excel('identified_proteins_coverage.xlsx', sheet_name='sheet1', index=False)


    num_list = [0]*101  #creat a list with 100 zeros, from 0 to 100%, 1% range
#protein_percentage_list = []

    for percentage in percentage_list:
        for i in range(0,101):
            if i<=percentage<i+1:
                num_list[i] += 1

    df1 = pd.DataFrame(dict(coverage_percentage=np.arange(0,101,1), number_of_proteins=num_list))

#for num in num_list:
 #   protein_percentage = float(num)/sum(num_list)*100
 #   protein_percentage_list.append(protein_percentage)
    sns.set()
    sns_plot = sns.lineplot(x="coverage_percentage", y="number_of_proteins", data=df1)
    fig = sns_plot.get_figure()

    plt.show()

#plt.xticks(np.arange(5,101,5))
#plt.yticks(np.arange(0,101,10))
#plt.show()