# The goal of this script is to generate decoy protein fasta file
# The method here is segmented reverse, which reverse all the sequence between every K and R to ensure the size of digested proteome is the same.
import re
import random
from expasy_rules import expasy_rules
from collections import Counter
import os
from collections import Counter



default_params = {
    'wrap_length': 60,
    'cleavage_rule': expasy_rules['trypsin_simple'],
    'specificity': 2,
    'miss_cleavage': 2,
    'min_len': 7,
    'max_len': 50,
}


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def chunk_string(string, length):
    return (string[0 + i:length + i] for i in range(0, len(string), length))


def fasta_protein_info_reader(fasta_file_path):
    protein_info_dict = {}
    length = []
    with open(fasta_file_path, 'r') as file_open:
        file_split = file_open.read().split('\n>')
        # print file_split[0]
        if file_split[0].startswith('>'):
            file_split[0] = file_split[0].replace('>', '')
            # print file_split[0]
        else:
            print("This fasta file contains other info in first line")
        for each in file_split:
            if '\n' in each:
                split_line = each.split('\n')
                protein_info_dict[split_line[0]] = ''.join(split_line[1:])
                if split_line[1] != '':
                    length.append(len(split_line[1]))
            else:
                print(each + "has no sequence attached")
        # print protein_info_dict
        global wrap_length
        wrap_length = Counter(length).most_common(1)[0][0]
        return protein_info_dict


def fasta_reverse_generator(fasta_file_in, fasta_file_out, reverse_mode=1, cut_pattern=expasy_rules['trypsin_simple']):
    """
    Generate fasta file with reversed sequence at the end,
    reverse_mode=1 -> whole sequence reverse,
    reverse_mode=2 -> random shuffle,
    reverse_mode=3 -> cut_pattern based segment reverse
    :param fasta_file_in: filename input
    :param fasta_file_out: filename output
    :param reverse_mode:
    :param cut_pattern: only works when reverse_mode==3
    :return : None
    """
    try:
        protein_info_dict = fasta_protein_info_reader(fasta_file_in)
        reversed_dict = {}
        for key in protein_info_dict:
            if reverse_mode == 1:
                reversed_dict[key] = protein_info_dict[key][::-1]
            elif reverse_mode == 2:
                shuffle_seq = list(protein_info_dict[key])
                random.shuffle(shuffle_seq)
                reversed_dict[key] = ''.join(shuffle_seq)
            else:
                reversed_dict[key] = ''.join([i[::-1] for i in re.split(cut_pattern, protein_info_dict[key])])

        with open(fasta_file_out, 'w', newline='\n') as fw:
            for key in protein_info_dict:
                fw.write('>%s\n' % key)
                fw.write('\n'.join(chunk_string(protein_info_dict[key], wrap_length)))
                fw.write('\n')
            for key in reversed_dict:
                fw.write('>Rev_%s\n' % key)
                fw.write('\n'.join(chunk_string(reversed_dict[key], wrap_length)))
                fw.write('\n')
    except Exception as e:
        print("Error parsing protein fasta file, failed to generate reversed fasta, error %s" % e)
        return 1
    return 0


def find_reverse_tag(filename):
    temp_list = []
    protein_info_dict = fasta_protein_info_reader(filename)
    protein_name_list = list(protein_info_dict.keys())
    length = len(protein_name_list)
    for i in reversed(range(3, 10)):
        prefix_list = [each[:i] for each in protein_name_list]
        temp_list += [(i[0], i[1] / length) for i in Counter(prefix_list).most_common(5) if abs(i[1] / length - 0.5) < 0.1 and i[0].endswith('_')]
    prefix_stat = sorted(temp_list, key=lambda x: x[1], reverse=True)
    return prefix_stat, protein_info_dict


def correct_fasta_rev_tag(filename):
    prefix_list, protein_info_dict = find_reverse_tag(filename)
    if len(prefix_list) == 0:
        fasta_reverse_generator(filename, filename, reverse_mode=1)
        prefix_list, protein_info_dict = find_reverse_tag(filename)
    longest_prefix = prefix_list[0][0]
    print(filename, longest_prefix, prefix_list)
    if longest_prefix != 'Rev_':
        with open(filename, 'w', newline='\n') as fw:
            for key in protein_info_dict:
                if key.startswith(longest_prefix):
                    fw.write('>%s\n' % key.replace(longest_prefix, 'Rev_'))
                else:
                    fw.write('>%s\n' % key)
                fw.write('\n'.join(chunk_string(protein_info_dict[key], wrap_length)))
                fw.write('\n')
    return os.path.abspath(filename)


def cleave_protein_sequence(seq):
    cleavage_rule = default_params['cleavage_rule']
    specificity = default_params['specificity']
    miss_cleavage = default_params['miss_cleavage']
    min_len = default_params['min_len']
    max_len = default_params['max_len']

    peptide_list = [m.end() for m in re.finditer(cleavage_rule, seq)]  # Create a list of cuts in the protein - does not include the start and end
    peptide_list.insert(0, 0)
    peptide_list.append(len(seq))  # add start and end to the cuts list
    print (peptide_list)
    peptides= set()
    if specificity == 1:  # if semi is selected in parameters
        """this generates a list of start and end points for peptides in the sequence in the format (start,end) 
        GOING FORWARDS from tryptic ends (starts at the first cut NOT 0)"""
        peptide_coordinates_forwards = [(peptide_list[j], peptide_list[j] + k)
                                        # format (start,end) iterating through k up to the next cut point
                                        for i in range(1, miss_cleavage + 2)  # (or one after up to the number of missed)
                                        for j in range(1, len(peptide_list) - i)
                                        # j corresponds to the current cut location in the cuts list and iterates up to the end - the number of missed cuts)
                                        for k in range((peptide_list[j + i - 1]) - peptide_list[j] + 1, (peptide_list[j + i] + 1) - peptide_list[j])
                                        # iterates though k from 1 to the next cut point - the range is modified to not repeat the ones from lower values of j
                                        if min_len <= k <= max_len]  # will only calculate lengths between the min and max length
        """this generates a list of start and end points for peptides in the sequence in the format (start,end) 
               GOING BACKWARDS from tryptic ends (starts at the first cut removes iterated values of k and stops at the last cut)"""
        peptide_coordinates_backwards = [((peptide_list[j + i] - k), peptide_list[j + i]) for i in range(1, miss_cleavage + 2) for j in
                                         range(0, (len(peptide_list) - i - 1)) for k in
                                         range(peptide_list[j + i] - peptide_list[j + 1] + 1, peptide_list[j + i] - (peptide_list[j]) + 1)
                                         if min_len <= k <= max_len]
        peptide_coordinates = peptide_coordinates_backwards + peptide_coordinates_forwards  # combines the Forward and Back coordinates
        peptides = {seq[i[0]:i[1]] for i in peptide_coordinates}  # uses coordinates to generate a dictionary of peptides from the sequence
    elif specificity == 2:  # if full is selected in parameters
        peptides = {seq[peptide_list[i]:peptide_list[i + j]] for j in range(1, miss_cleavage + 2) for i in range(0, len(peptide_list) - j) if
                    min_len <= (peptide_list[i + j] - peptide_list[i]) <= max_len}  # Creates list of coordinates similar to above, but without the iteration with k

    # elif specificity == 2:
    #     for j in range(1,miss_cleavage+2):
    #         for i in range(0,len(peptide_list)-j):
    #             print (peptide_list[i], peptide_list[i+j])
    #
    #             if min_len <= (peptide_list[i + j] - peptide_list[i]) <= max_len:
    #                 peptides.add(seq[peptide_list[i]:peptide_list[i + j]])

    else:  # NON tryptic
        peptide_coordinates_forwards = {seq[0:i] for i in range(1, len(seq) + 1) if min_len <= i <= max_len}
        peptide_coordinates_backwards = {seq[i:len(seq)] for i in range(1, len(seq))
                                         if min_len <= (len(seq) - i) <= max_len}
        peptides = peptide_coordinates_backwards | peptide_coordinates_forwards
    return peptides


def cleave_chunk(shared_list:list, spec_chunk:list):
    peptide_set=set()
    for each in spec_chunk:
        protein_id, seq = each
        peptide_set|=cleave_protein_sequence(seq)
    shared_list.append(peptide_set)
    return shared_list


def theoretical_cleavage_mp(fasta_dict:dict) -> list :
    import multiprocessing
    manager = multiprocessing.Manager()
    shared_list = manager.list()
    processes = [multiprocessing.Process(target=cleave_chunk, args=(shared_list, spec_chunk)) for spec_chunk in chunks(list(fasta_dict.items()), 500)]

    for p in processes:
        p.start()

    # Exit the completed processes
    for p in processes:
        p.join()

    #peptide_list = list(set(shared_list))
    #print(shared_list)
    return shared_list


if __name__ == '__main__':
    import time
    start = time.time()
    fasta_file_in = 'D:/data/proteome_fasta/new 4.fasta.txt'
    # fasta_file_out = r'D:\E\C\crux-3.0.Windows.AMD64\test\uniprot-human.fasta'
    # fasta_reverse_generator(fasta_file_in,'test_reverse.fasta')
    # fasta_reverse_generator(fasta_file_in, fasta_file_out, reverse_mode=1)
    # print(time.time() - start)
    fasta_dict= fasta_protein_info_reader(fasta_file_in)
    total_peptide_set=set().union(*theoretical_cleavage_mp(fasta_dict))
    print(len(total_peptide_set))
    print (total_peptide_set)
    if 'KTLLSNLEEAKK' in total_peptide_set:
        print (True)
    else:
        print (False)
    # print(Counter([len(i) for i in total_peptide_set]))
    #print(correct_fasta_rev_tag(fasta_file_out))
    print(time.time() - start)
