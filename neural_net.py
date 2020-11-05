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

def sequential_nn():
    model = keras.Sequential()
    model.add(layers.Embedding(21, 651, input_length=31))  # 21 integers encode, output 64 dimension, input 31 dimension
    model.add(layers.Flatten())
    model.add(layers.Dense(64, activation="relu", input_shape=(645,)))
    model.add(layers.Dense(64, activation="relu"))
    # model.add(layers.Dense(64, activation="relu"))
    model.add(layers.Dense(2, activation="softmax"))
    return model


def cnn():
    model = models.Sequential()
    model.add(layers.Embedding(30, 636, input_length=31))  # 21 integers encode, output 64 dimension, input 31 dimension
    model.add(layers.Reshape((-1,31,636)))
    model.add(layers.Conv2D(filters=64,kernel_size=1, activation='relu', input_shape=(15,43,1)))

    # model.add(layers.MaxPooling2D((2, 2)))
    # model.add(layers.Conv2D(filters=64,kernel_size=1, activation='relu'))
    # model.add(layers.MaxPooling2D((2, 2)))
    # model.add(layers.Conv2D(filters=64,kernel_size=1, activation='relu'))
    model.add(layers.Flatten())
    model.add(layers.Dense(64, activation='relu'))
    model.add(layers.Dense(2,activation='softmax'))
    return model


def lstm():
    model = models.Sequential()
    # model.add(layers.Embedding(30, 64, input_length=31))  # 21 integers encode, output 64 dimension, input 31 dimension
    model.add(layers.LSTM(64, input_shape=(1,646)))
    model.add(layers.Dense(64, activation='relu'))
    model.add(layers.Dense(2, activation='softmax'))
    return model


def compile_model(un_compiled_model):
    un_compiled_model.compile(
        optimizer=keras.optimizers.Adam(lr=0.0001),  # Optimizer
        # Loss function to minimize
        loss=keras.losses.sparse_categorical_crossentropy,
        # List of metrics to monitor
        metrics=[keras.metrics.SparseCategoricalAccuracy()],
    )
    complied_model = un_compiled_model
    print (complied_model.summary())
    return complied_model


# def auc(y_true, y_pred):
#     auc = tf.metrics.auc(y_true, y_pred)[1]
#     backend.get_session().run(tf.local_variables_initializer())
#     return auc
# model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])

t_37C_240min_dict = ppp.load(open('D:/data/deep_proteome/pickle_file/20200915_tryp_37C_1440min_new.p','rb'))
test_dataset_dict = ppp.load(open('P62908_polymer_dict.p','rb'))
predict_matrix_array = ppp.load(open('P62908_matrix_2d_array.p', 'rb'))

# print (Counter([t_37C_240min_dict[each] for each in t_37C_240min_dict]))
# df_dummy = df_dummy_getter(t_37C_240min_dict)
# matrix, target = matrix_target_getter(df_dummy)
# matrix, target = dump_data(t_37C_240min_dict)
# print (matrix)

matrix,target = matrix_target(t_37C_240min_dict)
encoder,matrix = ohe(matrix)
# matrix = csr_matrix.toarray(matrix)
print (matrix.shape)
# test set from different dataset
test_maxtrix, test_target = matrix_target(test_dataset_dict)
test_maxtrix = encoder.transform(test_maxtrix)
predict_matrix = encoder.transform(predict_matrix_array)
# test_maxtrix = csr_matrix.toarray(test_maxtrix)

# one-hot encode, optional
# vocab_size = 21
# encoded_docs = [one_hot(seq, vocab_size) for seq in matrix]
# encoded_docs = np.array(encoded_docs)
# print (encoded_docs)

X_train, X_test, target_train, target_test = train_test_data_split(matrix,target) # matrix is 2d
# print (type(X_train))

# convert into 4d array for CNN
# X_train, X_test = X_train.values.reshape(X_train.shape[0],15,43,1), X_test.values.reshape(X_test.shape[0],15,43,1)

# convert in to 3d array for lstm
X_train, X_test, test_maxtrix,predict_matrix = X_train.reshape(X_train.shape[0],1,X_train.shape[1]), \
                                X_test.reshape(X_test.shape[0],1,X_test.shape[1]), \
                                test_maxtrix.reshape(test_maxtrix.shape[0],1,test_maxtrix.shape[1]), \
                                predict_matrix.reshape(predict_matrix.shape[0],1,predict_matrix.shape[1])


print (X_train.shape)
X_val = X_train[-300:]
y_val = target_train[-300:]

X_train = X_train[:-300]
y_train = target_train[:-300]

# class weights consideration
# from sklearn.utils import class_weight
# class_weights = class_weight.compute_class_weight('balanced',
#                                                  np.unique(y_train),
#                                                  y_train)
# class_weight_dict = dict(enumerate(class_weight))


model = compile_model(lstm())

print("Fit model on training data")
start = time.time()
history = model.fit(
    X_train,
    y_train,
    batch_size=64,
    epochs=20,
    # We pass some validation for
    # monitoring validation loss and metrics
    # at the end of each epoch
    # an error is calculated at the end of each batch, the error is used to improve the model, eg. move down the
    # error gradient
    validation_split=0.2,
    shuffle=True
)
print ("model trained time: ", time.time()-start)
print("Evaluate on test data")
results = model.evaluate(test_maxtrix, test_target, batch_size=64)
print("test loss, test acc:", results)
print (np.argmax(model.predict(predict_matrix), axis=-1))

# print (history.history)
# plt.plot(history.history['val_loss'], label='val_loss')
# plt.plot(history.history['val_sparse_categorical_accuracy'], label = 'val_accuracy')
# plt.xlabel('Epoch')
# plt.ylabel('Accuracy')
# # plt.ylim([0.5, 1])
# plt.legend(loc='lower right')
# plt.show()

# plot loss during training
plt.subplot(211)
plt.title('Loss')
plt.plot(history.history['loss'], label='train')
plt.plot(history.history['val_loss'], label='validation')
plt.legend(loc='upper right')
plt.xlim(0,20)
plt.ylim(0,1)
# plot accuracy during training
plt.subplot(212)
plt.title('Accuracy')
plt.plot(history.history['sparse_categorical_accuracy'], label='train')
plt.plot(history.history['val_sparse_categorical_accuracy'], label='validation')
plt.legend(loc='lower right')
plt.xlim(0,20)
plt.ylim(0,1)
plt.xlabel('Epoch')
plt.show()


# predict probabilities for test set
# yhat_probs = model.predict(X_test, verbose=0)
# print (yhat_probs)
# predict crisp classes for test set
yhat_classes = np.argmax(model.predict(test_maxtrix), axis=-1)
# print (yhat_classes.shape)
# print (model.predict(X_test))
yhat_score = model.predict(test_maxtrix)[:,-1]
# print (yhat_score)
fpr_keras, tpr_keras, thresholds_keras = roc_curve(test_target, yhat_score)
auc_keras = auc(fpr_keras, tpr_keras)
print ('AUC: %f' % auc_keras)
# precision tp / (tp + fp)
precision = precision_score(test_target, yhat_classes)
print('Precision: %f' % precision)
# recall: tp / (tp + fn)
recall = recall_score(test_target, yhat_classes)
print('Recall: %f' % recall)
# f1: 2 tp / (2 tp + fp + fn)
f1 = f1_score(test_target, yhat_classes)
print('F1 score: %f' % f1)

classify_report = classification_report(test_target,yhat_classes)
print (classify_report)

plt.plot(fpr_keras, tpr_keras, color='darkorange',
             lw=2, label='ROC curve (area = %0.2f)' % auc_keras)
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC')
plt.legend(loc="lower right")
plt.show()

