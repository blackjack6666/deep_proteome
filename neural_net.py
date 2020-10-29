import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers, models
from stat_models import df_dummy_getter, matrix_target_getter, train_test_data_split
import pickle as ppp
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, cohen_kappa_score, roc_auc_score, confusion_matrix


def sequential_nn():
    model = keras.Sequential()
    model.add(layers.Dense(128, activation="relu", input_shape=(645,)))
    model.add(layers.Dense(64, activation="relu"))
    model.add(layers.Dense(64, activation="relu"))
    model.add(layers.Dense(16, activation="softmax"))
    return model


def cnn():
    model = models.Sequential()
    model.add(layers.Conv2D(filters=64,kernel_size=1, activation='relu', input_shape=(15,43,1)))
    # model.add(layers.MaxPooling2D((2, 2)))
    # model.add(layers.Conv2D(filters=64,kernel_size=1, activation='relu'))
    # model.add(layers.MaxPooling2D((2, 2)))
    # model.add(layers.Conv2D(filters=64,kernel_size=1, activation='relu'))
    model.add(layers.Flatten())
    model.add(layers.Dense(64, activation='relu'))
    model.add(layers.Dense(32,activation='softmax'))
    return model


def lstm():
    model = models.Sequential()
    model.add(layers.LSTM(64, input_shape=(1,645)))
    model.add(layers.Dense(64, activation='relu'))
    model.add(layers.Dense(32, activation='softmax'))
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
    return complied_model

# model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])

t_37C_240min_dict = ppp.load(open('tryps_37C_240min_cleavage_label.p','rb'))
# print (Counter([t_37C_240min_dict[each] for each in t_37C_240min_dict]))
df_dummy = df_dummy_getter(t_37C_240min_dict)
matrix, target = matrix_target_getter(df_dummy)
X_train, X_test, target_train, target_test = train_test_data_split(matrix,target) # matrix is 2d

# convert into 4d array for CNN
# X_train, X_test = X_train.values.reshape(X_train.shape[0],15,43,1), X_test.values.reshape(X_test.shape[0],15,43,1)

# convert in to 3d array for lstm
X_train, X_test = X_train.values.reshape(X_train.shape[0],1,X_train.shape[1]), X_test.values.reshape(X_test.shape[0],1,X_test.shape[1])

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
history = model.fit(
    X_train,
    y_train,
    batch_size=32,
    epochs=20,
    # We pass some validation for
    # monitoring validation loss and metrics
    # at the end of each epoch
    # an error is calculated at the end of each batch, the error is used to improve the model, eg. move down the
    # error gradient
    validation_split=0.2,
    shuffle=True
)

print("Evaluate on test data")
results = model.evaluate(X_test, target_test, batch_size=32)
print("test loss, test acc:", results)

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
plt.ylim(0,4)
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
yhat_classes = np.argmax(model.predict(X_test), axis=-1)
# print (yhat_classes)
# reduce to 1d array
# yhat_probs = yhat_probs[:, 0]

# precision tp / (tp + fp)
precision = precision_score(target_test, yhat_classes)
print('Precision: %f' % precision)
# recall: tp / (tp + fn)
recall = recall_score(target_test, yhat_classes)
print('Recall: %f' % recall)
# f1: 2 tp / (2 tp + fp + fn)
f1 = f1_score(target_test, yhat_classes)
print('F1 score: %f' % f1)
