from collections import defaultdict
import glob
import numpy as np


def read_fasta_info_dict2(filename):
    """
    -----
    read fasta file including reverse seq, exclude rev seq in returned dictionary
    -----
    :param filename:
    :return:
    """
    protein_seq_dict = {}
    with open(filename, 'r') as file_open:
        file_split = file_open.read().split('\n>')
        # read first entry
        protein_seq_dict[file_split[0].split('\n')[0].split('|')[1]] = ''.join(file_split[0].split('\n')[1:])
        # read the rest entry
        for each in file_split[1:]:
            if each.startswith('Rev'):
                continue
            else:
                split_line = each.split('\n')
                protein_seq_dict[split_line[0].split('|')[1]] = ''.join(split_line[1:])

    return protein_seq_dict


def fasta_reader(fasta_file_path):

    with open(fasta_file_path, 'r') as file_open:
        file_split = file_open.read().split('\n>')

    return {each.split('\n')[0].split('|')[1]: ''.join(each.split('\n')[1:]) for each in file_split}


def fasta_reader2(fasta_file_path):
    """
    read fasta file without reverse seq
    :param fasta_file_path:
    :return:
    """
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
    """
       read description and prefix in fasta file
       :param fasta_file:
       :return:
       """

    # dictionary format: {'ID': ('sp', 'description')}
    ID_description_dict = {}
    with open(fasta_file, 'r') as file_open:
        file_split = file_open.read().split('\n>')
        # first entry special case
        ID_description_dict[file_split[0].split('\n')[0].split('|')[1]] = (
        file_split[0].split('\n')[0].split('|')[0][1:], file_split[0].split('\n')[0].split('|')[2])
        for each in file_split[1:]:
            split_line = each.split('\n')
            ID_description_dict[split_line[0].split('|')[1]] = (
            split_line[0].split('|')[0], split_line[0].split('|')[2])
    return ID_description_dict


def fasta_reverse_generator(fasta_file_in, fasta_file_out):
    print ('reverse_algorithm = null, or other protease')
    # read protein sequence into dic
    protein_dict = fasta_reader(fasta_file_in)

    # read description into dic
    ID_descrip_dict = read_description(fasta_file_in)  # dictionary format: {'ID': ('sp'or'tr', description)}

    # write id and reverse sequence into fasta_file_out
    with open(fasta_file_out, 'w', newline='\n') as file_open:
        for id in protein_dict:
            forward_seq = protein_dict[id]

            rev_seq = protein_dict[id][::-1]

            block = range(0,len(forward_seq)+60,60)

            # write forward
            file_open.write('>'+ ID_descrip_dict[id][0]+'|'+id+'|'+ID_descrip_dict[id][1]+'\n')
            for i in range(len(block)-1):
                file_open.write(forward_seq[block[i]:block[i+1]]+'\n')

            # write reverse
            file_open.write('>Rev_'+ ID_descrip_dict[id][0]+'|'+id+'|'+ID_descrip_dict[id][1]+'\n')
            for i in range(len(block)-1):
                file_open.write(rev_seq[block[i]:block[i+1]]+'\n')
    return fasta_file_out


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


if __name__ == "__main__":
    fasta_file_input = 'D:/data/Naba_deep_matrisome/mouse_ecm_costom_proteome_db.fasta'
    fasta_file_out = 'D:/data/Naba_deep_matrisome/mouse_ecm_costom_proteome_db_rev.fasta'
    fasta_reverse_generator(fasta_file_input,fasta_file_out)

