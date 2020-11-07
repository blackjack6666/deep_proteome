# enzyme specificity

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
    :param matrix: [[31mer1],[31mer2]...]
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

t_37C_240min_dict = ppp.load(open('D:/data/deep_proteome/pickle_file/20200915_tryp_37C_1440min_new.p','rb'))
matrix,target = matrix_target(t_37C_240min_dict)
matrix = matrix[:5]

print (custom_ohe(matrix))