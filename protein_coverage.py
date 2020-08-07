from collections import defaultdict
import glob
import numpy as np
def read_fasta_into_dict(filename):
    dict_human_ref = dict()
    #key = ""
    value = ""
    aa_frequency_all_zero_dict = {}
    with open(filename, 'r') as file_open:
        Reverse = 0
        for line in file_open:
            if line.startswith('>sp') and value == '':
                key = line.split('|')[1]
                Reverse = 0
                #dict_human_ref[key] = value
                #value = ""
                #key = line.split('|')[1]
            elif line.startswith('>sp') and value != '':
                dict_human_ref[key] = value
                key = line.split('|')[1]
                value = ''
                Reverse = 0
            elif line.startswith('>Rev'):
                Reverse = 1
            elif Reverse == 0:
                value += line.rstrip('\n')

    if key not in dict_human_ref.keys():
        dict_human_ref[key] = value
    for ID in [ID for ID in dict_human_ref.keys()][1:]:
        sequence = dict_human_ref[ID]
        aa_frequency_all_zero_dict[ID] = [0]*len(sequence)


    return dict_human_ref, aa_frequency_all_zero_dict


def fasta_reader(fasta_file_path):
    protein_seq_dict = {}
    with open(fasta_file_path, 'r') as file_open:
        file_split = file_open.read().split('\n>')
        for each in file_split:
            split_line = each.split('\n')
            protein_seq_dict[split_line[0].split('|')[1]] = ''.join(split_line[1:])
    return protein_seq_dict

def read_description_into_dict(fasta_file):
    ID_description_dict = {}
    with open(fasta_file, 'r') as file_open:
        for line in file_open:
            if line.startswith('>'):
                ID_description_dict[line.split('|')[1]] = line.split('|')[2].rstrip('\n')
    return ID_description_dict

def read_description(fasta_file):
    # dictionary format: {'ID': ('sp', 'description')}
    ID_description_dict ={}
    with open(fasta_file, 'r') as file_open:
        file_split = file_open.read().split('\n>')
        for each in file_split:
            split_line = each.split('\n')
            ID_description_dict[split_line[0].split('|')[1]] = (split_line[0].split('|')[0],split_line[0].split('|')[2])
    return ID_description_dict

def contaminant_converter(contaminant_file, normalized_contaminant_file):
    contaminant_description_dict = {}
    with open(contaminant_file, 'r') as file_open:
        ID = ''
        value = ''
        description = ''
        for line in file_open:
            if line.startswith('>contaminant'):
                contaminant_description_dict[ID] = (value, description)
                value = ''
                ID = line.split('|')[0]
                description = line.split('|')[1]
            else:
                value += line.rstrip('\n')
    if ID not in contaminant_description_dict.keys():
        contaminant_description_dict[ID] = (value, description)
    with open(normalized_contaminant_file, 'w') as new_file:
        for ID in contaminant_description_dict.keys()[1:]:
            new_file.write('>CT|'+ID.lstrip('>')+'|'+contaminant_description_dict[ID][1])
            interval = np.arange(0, len(contaminant_description_dict[ID][0]), 60)
            interval = np.append(interval, len(contaminant_description_dict[ID][0]))
            for i in range(len(interval)-1):
                new_file.write(contaminant_description_dict[ID][0][interval[i]:interval[i+1]]+'\n')
    return  new_file




def aa_frequency(filename):
    fasta_name = 'C:/uic/lab/data/xinhao_data1/uniprot-proteome_UP000005640.fasta'
    aa_count_dic = {}
    identified_count_dict = {}
    ID = ''
    pep_dict = {}
    pep_seq_list = []
    protein_dic = read_fasta_into_dict(fasta_name)
    Total_protein_len = 0
    percentage_sum = 0
    Identified_total_protein_len = 0
    Total_identified_zero_count = 0
    with open(filename) as fileopen:
        for i in range(29):
            next(fileopen)
        Reverse_start = 0
        for line in fileopen:
            if line.startswith('Reverse_'):
                Reverse_start = 1
            elif line.startswith('sp'):
                Reverse_start = 0
                pep_dict[ID] = pep_seq_list
                pep_seq_list = ""
                ID = line.split('\t')[0].split('|')[1]


            elif len(line.split('\t')) == 15 and Reverse_start == 0:
                peptide = line.split('\t')[-1].split('.')[1]

                pep_seq_list += peptide



    if ID not in pep_dict:
        pep_dict[ID] = pep_seq_list

        for ID in protein_dic.keys():
            sequence = protein_dic[ID]
            aa_count_dic[ID] = [0]*len(sequence)
            Total_protein_len += len(sequence)
            if ID in pep_dict:
                identified_count_dict[ID] = [0]*len(sequence)
                Identified_total_protein_len += len(sequence)

                for pep in pep_dict[ID]:
                    start_position = sequence.find(pep)
                    end_position = start_position + len(pep)
                    aa_count_dic[ID][start_position:end_position] = [x+1 for x in aa_count_dic[ID][start_position:end_position]]
                    identified_count_dict[ID][start_position:end_position] = [x+1 for x in identified_count_dict[ID][start_position:end_position]]

        for ID in identified_count_dict.keys()[1:]:
            identified_zero_count = identified_count_dict[ID].count(0)
            Total_identified_zero_count += identified_zero_count
            total_len_for_one_ID = len(identified_count_dict[ID])
            percentage = float(total_len_for_one_ID-identified_zero_count)/total_len_for_one_ID
            percentage_sum += percentage
    average = float(percentage_sum/len(identified_count_dict.keys()[1:]))
    return average, Total_protein_len, Identified_total_protein_len, Total_identified_zero_count

if __name__ == "__main":
    Alldta_total_protein_len = 0
    Alldta_total_identified_protein_len = 0
    Alldta_total_identified_zero_count = 0
    Alldta_average_list = []
    path = 'C:/uic/lab/data/xinhao_data1/'
    datafiles = glob.glob(path+'*.dta')
    for datafile in datafiles:
        Alldta_total_protein_len += aa_frequency(datafile)[1]
        Alldta_total_identified_protein_len += aa_frequency(datafile)[2]
        Alldta_total_identified_zero_count += aa_frequency(datafile)[3]
        Alldta_average_list.append(aa_frequency(datafile)[0])
    overall_average = np.mean(Alldta_average_list)*100
    Overall_coverage = float(Alldta_total_identified_protein_len-Alldta_total_identified_zero_count)/Alldta_total_protein_len*100
    Overall_coverage_for_identified_protein = float(Alldta_total_identified_protein_len-Alldta_total_identified_zero_count)/Alldta_total_identified_protein_len*100

    print ('The average is %.2f%%\nThe overall coverage for all proteins is %.2f%%\nThe overall coverage for identified proteins is %.2f%%')\
        % (overall_average, Overall_coverage, Overall_coverage_for_identified_protein )

