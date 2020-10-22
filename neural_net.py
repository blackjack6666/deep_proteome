import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers, models
from stat_models import df_dummy_getter, matrix_target_getter, train_test_data_split
import pickle as ppp
import matplotlib.pyplot as plt


def sequential_nn():
    model = keras.Sequential()
    model.add(layers.Dense(64, activation="relu", input_shape=(645,)))
    model.add(layers.Dense(128, activation="relu"))
    model.add(layers.Dense(64, activation="relu"))
    model.add(layers.Dense(32, activation="softmax"))
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
    model.add(layers.Dense(10,activation='softmax'))
    return model


def compile_model(un_compiled_model):
    un_compiled_model.compile(
        optimizer=keras.optimizers.Adam(),  # Optimizer
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
X_train, X_test = X_train.values.reshape(X_train.shape[0],15,43,1), X_test.values.reshape(X_test.shape[0],15,43,1)


X_val = X_train[-500:]
y_val = target_train[-500:]

X_train = X_train[:-500]
y_train = target_train[:-500]

model = compile_model(cnn())

print("Fit model on training data")
history = model.fit(
    X_train,
    y_train,
    batch_size=32,
    epochs=10,
    # We pass some validation for
    # monitoring validation loss and metrics
    # at the end of each epoch
    validation_data=(X_val, y_val),
)

print("Evaluate on test data")
results = model.evaluate(X_test, target_test, batch_size=32)
print("test loss, test acc:", results)

print (history.history)
plt.plot(history.history['accuracy'], label='accuracy')
plt.plot(history.history['val_accuracy'], label = 'val_accuracy')
plt.xlabel('Epoch')
plt.ylabel('Accuracy')
plt.ylim([0.5, 1])
plt.legend(loc='lower right')