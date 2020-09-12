import ahocorasick
#from test4 import read_peptide
from protein_coverage import read_fasta_info_dict2
import glob
import time
import multiprocessing_naive_algorithym
#import calculations_and_plot
import seaborn as sns
from matplotlib import pyplot as plt
import pickle as ppp

def automaton_trie(peptide_list):
    A = ahocorasick.Automaton()
    for idx, peptide in enumerate(peptide_list):
        A.add_word(peptide, (idx, peptide))
    A.make_automaton()
    return A

def automaton_matching(A, seq_line):
    result = []
    for end_idx, (insert_order, original_value) in A.iter(seq_line):
        start_idx = end_idx - len(original_value) + 1
        result.append((start_idx, end_idx, original_value))
        assert seq_line[start_idx:start_idx+len(original_value)] == original_value
    return result


if __name__=='__main__':

    path = 'C:/uic/lab/data/dta_result/'
    dtafiles = glob.glob(path + '*.dta')
    start = time.clock()
    peptide_list = read_peptide(dtafiles)
    peptide_list_set = set(peptide_list)
    peptide_list_unique = list(peptide_list_set)
    print (time.clock()-start)
    automaton = automaton_trie(peptide_list_unique)
    print (time.clock()-start)

    filename = 'C:/uic/lab/data/xinhao_data1/uniprot-proteome_UP000005640.fasta'
    protein_dict = read_fasta_into_dict(filename)[0]

    ID_list, seq_list = multiprocessing_naive_algorithym.extract_UNID_and_seq(protein_dict)

    seq_line = multiprocessing_naive_algorithym.creat_total_seq_line(seq_list)
    sep_pos_array = multiprocessing_naive_algorithym.separator_pos(seq_line)
    zero_line = multiprocessing_naive_algorithym.zero_line_for_seq(seq_line)
    print (time.clock()-start)
    result = automaton_matching(automaton, seq_line)
    print (time.clock()-start)
    for pos in result:
        zero_line[pos[0]:pos[1]+1] += 1 # pos[0] is start position, pos[1] is end position



    #read zero_line, seq_line, ID_list, seq_list, sep_pos_array into pickle file.
    ppp.dump(zero_line, open('number_line_5_5.p', 'wb'), protocol=-1)
    ppp.dump(seq_line, open('seq_line_5_5.p', 'wb'), protocol=-1)
    ppp.dump(ID_list, open('ID_list_5_5.p', 'wb'), protocol=-1)
    ppp.dump(seq_list, open('seq_list_5_5.p', 'wb'), protocol=-1)
    ppp.dump(sep_pos_array, open('sep_pos_5_5.p', 'wb'), protocol=-1)
    ppp.dump(result, open('aho_result_5_5.p', 'wb'), protocol=-1)


    # do calculations in the following two lines.
    overall_percentage, identified_coverage_array, coverage_list_ordered, protein_coverage_dict, average_coverage, identified_average_coverage, len_list = \
    calculations_and_plot.coverage_calculation(zero_line, sep_pos_array, ID_list)

    print ('the overall average is %s, the average coverage is %s, the identified average coverage is %s' \
          % (overall_percentage, average_coverage, identified_average_coverage))
    ID_ratio_dict = calculations_and_plot.ratio_of_coverage_to_sum(zero_line, sep_pos_array, ID_list)
    #print ID_ratio_dict['Q8WZ42-13'], ID_ratio_dict['Q8WXI7'], ID_ratio_dict['Q8NF91']
    ID_length_cov_dict = {}
    for ID, length, cov in zip(ID_list, len_list, coverage_list_ordered):
        ID_length_cov_dict[ID] = (length, cov)
    for ID in ID_length_cov_dict:
       if 80.0<ID_length_cov_dict[ID][1] and ID_length_cov_dict[ID][0]>5500.0:
            print (ID, ID_length_cov_dict[ID])
    df1 = calculations_and_plot.read_into_dataframe(protein_coverage_dict, identified_coverage_array)[2]

    g = sns.FacetGrid(df1, size=6)
    g = g.map(plt.plot, 'coverage', 'number_of_identified_proteins')
    plt.show()
    '''
    
    df = calculations_and_plot.protein_length_vs_coverage(len_list, coverage_list_ordered)[1]
  
    sns.set()
    sns.scatterplot(x='protein_length', y='coverage', data=df)
    plt.show()
    '''