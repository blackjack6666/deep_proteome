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


# test the built model on different dataset

model = keras.models.load_model('D:/data/deep_proteome/deep_learning_models/tcn_pombe_gluc_tryps_15epoch',custom_objects={'TCN': TCN})

test_data = ppp.load(open('D:/data/deep_proteome/non_specfic_search/pombe_gluc_tryps_rep2_31mer.p','rb'))

# load keras model and training dataset
matrix,target = matrix_target(test_data)
matrix = custom_ohe(matrix,polymer_len=15)


results = model.evaluate(matrix, target, batch_size=64)

yhat_classes = np.argmax(predict(matrix,model=model), axis=-1)  # use when model produce 1D output (sigmoid, units=1, probability for class 1)

print (yhat_classes.shape)
# print (model.predict(X_test))
# print (model.predict(test_maxtrix))
yhat_score = model.predict(matrix)[:,-1]  # probability of class 1 for each instance


print (yhat_score)
fpr_keras, tpr_keras, thresholds_keras = roc_curve(target, yhat_score)
auc_keras = auc(fpr_keras, tpr_keras)
print ('AUC: %f' % auc_keras)
# precision tp / (tp + fp)
precision = precision_score(target, yhat_classes)
print('Precision: %f' % precision)
# recall: tp / (tp + fn)
recall = recall_score(target, yhat_classes)
print('Recall: %f' % recall)
# f1: 2 tp / (2 tp + fp + fn)
f1 = f1_score(target, yhat_classes)
print('F1 score: %f' % f1)

classify_report = classification_report(target,yhat_classes)
print (classify_report)

plt.plot(fpr_keras, tpr_keras, color='darkorange',
             lw=2, label='ROC curve (area = %0.3f)' % auc_keras)
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC')
plt.legend(loc="lower right")
plt.show()