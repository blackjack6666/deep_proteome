import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers, models,backend
from tensorflow.keras.preprocessing.text import one_hot
from stat_models import df_dummy_getter, matrix_target_getter, train_test_data_split,dump_data, matrix_target,ohe
import pickle as ppp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report, precision_score, recall_score, f1_score, cohen_kappa_score, roc_curve,auc, confusion_matrix
import time
from scipy.sparse import csr_matrix
from parameters import custom_ohe
from sklearn.utils import class_weight
from tcn import TCN
from neural_net import predict
from tsv_reader import logomaker_from_stop_codon, heatmap_gen


# test the built model on different dataset

model = keras.models.load_model('D:/data/deep_proteome/deep_learning_models/tryp_30_thermo_30',custom_objects={'TCN': TCN})

test_data = ppp.load(open('D:/data/deep_proteome/non_specfic_search/tryp_30_thermo_30.p','rb'))

# load keras model and training dataset
matrix_seq,target = matrix_target(test_data)

matrix = custom_ohe(matrix_seq,polymer_len=15)


# results = model.evaluate(matrix, target, batch_size=64)

### class analysis
yhat_classes = np.argmax(predict(matrix,model=model), axis=-1)  # use when model produce 1D output (sigmoid, units=1, probability for class 1)
yhat_class_1_index =np.nonzero(yhat_classes)[0]

print (yhat_classes.shape)
print (yhat_classes)
print (yhat_class_1_index)
seq_class1_list = matrix_seq[yhat_class_1_index, :]
seq_class1_list = [''.join(each) for each in seq_class1_list]
filtered_seq_class1_list = [each for each in seq_class1_list if each[7] == 'E']
print (len(seq_class1_list))
prob_matrix = logomaker_from_stop_codon(seq_class1_list)
prob_matrix = prob_matrix.T*100
# heatmap_gen(prob_matrix)

# show logo for prediction of missed cleavaged polymer
"""
yhat_class_0_index = np.where(yhat_classes==0)[0]
seq_class0_list = matrix_seq[yhat_class_0_index, :]
seq_class0_list = [''.join(each) for each in seq_class0_list]
filtered_seq_class0_list = [each for each in seq_class0_list if each[8]=='F']
prob_matrix = logomaker_from_stop_codon(filtered_seq_class0_list)
prob_matrix = prob_matrix.T*100
ax = heatmap_gen(prob_matrix)
# print (model.predict(X_test))
# print (model.predict(test_maxtrix))
"""




### statistical plotting
# yhat_score = model.predict(matrix)[:,-1]  # probability of class 1 for each instance
# print (yhat_score)
# fpr_keras, tpr_keras, thresholds_keras = roc_curve(target, yhat_score)
# auc_keras = auc(fpr_keras, tpr_keras)
# print ('AUC: %f' % auc_keras)
# # precision tp / (tp + fp)
# precision = precision_score(target, yhat_classes)
# print('Precision: %f' % precision)
# # recall: tp / (tp + fn)
# recall = recall_score(target, yhat_classes)
# print('Recall: %f' % recall)
# # f1: 2 tp / (2 tp + fp + fn)
# f1 = f1_score(target, yhat_classes)
# print('F1 score: %f' % f1)
#
# classify_report = classification_report(target,yhat_classes)
# print (classify_report)
#
# plt.plot(fpr_keras, tpr_keras, color='darkorange',
#              lw=2, label='ROC curve (area = %0.3f)' % auc_keras)
# plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
# plt.xlim([0.0, 1.0])
# plt.ylim([0.0, 1.05])
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('ROC')
# plt.legend(loc="lower right")
# plt.show()