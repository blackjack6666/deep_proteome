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

if __name__ == '__main__':
    t_37C_240min_dict = ppp.load(open('tryp_24h_label_dict_11_8.p','rb'))
    matrix,target = matrix_target(t_37C_240min_dict)
    matrix = matrix[0]
    print (matrix)

    print (array_to_seq(custom_ohe([matrix])[0]))