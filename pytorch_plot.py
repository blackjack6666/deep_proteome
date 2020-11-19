import matplotlib.pyplot as plt
import pickle as ppp

history = ppp.load(open('pytorch_history.p','rb'))
train_loss, val_loss, acc = [], [], []
for each in history:
    train_loss.append(each[1])
    val_loss.append(each[2])
    acc.append(each[3])
plt.subplot(211)
plt.plot(train_loss, label='train_loss')
plt.plot(val_loss, label='validation_loss')
plt.legend(loc='upper right')
plt.subplot(212)
plt.plot(acc, label='accuracy')
plt.legend(loc='upper right')
plt.show()