# enzyme specificity
import numpy as np
enzyme_spec = {'trypsin':'KR','chymotrypsin':'FWY'}

aa_interger_dict = {
    'A':0,
    'R':1,
    'N':2,
    'D':3,
    'C':4,
    'Q':5,
    'E':6,
    'G':7,
    'H':8,
    'I':9,
    'L':10,
    'K':11,
    'M':12,
    'F':13,
    'P':14,
    'S':15,
    'T':16,
    'W':17,
    'Y':18,
    'V':19,
    'X':20,
    'Z':21
}
# https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php
hydrophobicity_score_dict = {
    'A':47,
    'R':-26,
    'N':-41,
    'D':-18,
    'C':52,
    'Q':-18,
    'E':8,
    'G':0,
    'H':-42,
    'I':100,
    'L':100,
    'K':-37,
    'M':74,
    'F':92,
    'P':-46,
    'S':-7,
    'T':13,
    'W':84,
    'Y':49,
    'V':79,
    'X':0,
    'Z':0
}

import pickle as ppp
from stat_models import matrix_target


def custom_ohe(matrix):
    """
    customized one hot encoder
    :param matrix: [[31mer1],[31mer2]...], 2d-array
    :return:
    """

    import numpy as np
    encoding_matrix = []
    for each_mer in matrix:
        array = np.zeros(31*22)
        for ind,val in enumerate(each_mer):
            array[ind*22+aa_interger_dict[val]]=1
        encoding_matrix.append(array)
    encoding_matrix = np.array(encoding_matrix)
    return encoding_matrix


def array_to_seq(array):
    """

    :param array: one-d numpy array
    :return:
    """
    aa_seq = ""
    int_aa_dict = {aa_interger_dict[each]:each for each in aa_interger_dict}
    array_bin = range(0,len(array)+22,22)
    for i in range(len(array_bin)-1):
        sub_array = array[array_bin[i]:array_bin[i+1]]
        aa = int_aa_dict[np.nonzero(sub_array)[0][0]]
        aa_seq+=aa
    return aa_seq


def hydrophobicity_cal(matrix):
    """
    calculate a hydrophobicity score for each 31mer, take the middle 10 aa as input
    :param matrix: [[31mer1],[31mer2]...], 2d-array
    :return:
    """
    import numpy as np
    hydro_matrix = []

    # max hydrophobicity score is 11*100=1100, min hydrophobicity score is 11*-46 = -506
    bins = np.linspace(-506,1100,10)
    #print (bins)
    for each_mer in matrix:
        hy_array = np.zeros(10)
        # 11-mer
        ele_mer = each_mer[10:21]
        hydro_score = sum([hydrophobicity_score_dict[each] for each in ele_mer])
        # one hot encoding of hydro score
        hy_array[np.digitize(hydro_score, bins)] = 1
        hydro_matrix.append(hy_array)
    hydro_matrix = np.array(hydro_matrix)
    return hydro_matrix


def matrix_addup(one_hot_matrix,hydropho_matrix):
    """
    concatenate two numpy 2d arrays
    :param one_hot_matrix: returned from custome_ohe
    :param hydropho_matrix: returned from hydrophobicity cal
    :return: a combined numpy 2d arrays, shape[0] would be same, shape[1] would be added up
    """
    return np.concatenate((one_hot_matrix,hydropho_matrix),axis=1)

if __name__ == '__main__':
    t_37C_240min_dict = ppp.load(open('tryp_24h_label_dict_11_8.p','rb'))
    matrix,target = matrix_target(t_37C_240min_dict)
    one_hot_matrix = custom_ohe(matrix)
    hydro_matrix = hydrophobicity_cal(matrix)
    new_matrix = matrix_addup(one_hot_matrix,hydro_matrix)

    # matrix = matrix[0]
    # print (matrix)
    #
    # print (array_to_seq(custom_ohe([matrix])[0]))