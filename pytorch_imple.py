import pickle as ppp
from stat_models import matrix_target, train_test_data_split
from parameters import custom_ohe
from torch.utils.data import TensorDataset, DataLoader
import torch
from torch import nn
import numpy as np
from torch.optim.lr_scheduler import _LRScheduler
from torch.nn import functional as F


class LSTMClassifier(nn.Module):
    """Very simple implementation of LSTM-based time-series classifier."""

    def __init__(self, input_dim, hidden_dim, layer_dim, output_dim):
        super().__init__()
        self.hidden_dim = hidden_dim
        self.layer_dim = layer_dim
        self.rnn = nn.LSTM(input_dim, hidden_dim, layer_dim, batch_first=True)
        self.fc = nn.Linear(hidden_dim, output_dim)
        self.batch_size = None
        self.hidden = None

    def forward(self, x):
        h0, c0 = self.init_hidden(x)
        out, (hn, cn) = self.rnn(x, (h0, c0))
        out = self.fc(out[:, -1, :])
        return out

    def init_hidden(self, x):
        h0 = torch.zeros(self.layer_dim, x.size(0), self.hidden_dim)
        c0 = torch.zeros(self.layer_dim, x.size(0), self.hidden_dim)
        return [t.cuda() for t in (h0, c0)]


class CyclicLR(_LRScheduler):

    def __init__(self, optimizer, schedule, last_epoch=-1):
        assert callable(schedule)
        self.schedule = schedule
        super().__init__(optimizer, last_epoch)

    def get_lr(self):
        return [self.schedule(self.last_epoch, lr) for lr in self.base_lrs]


def cosine(t_max, eta_min=0):
    def scheduler(epoch, base_lr):
        t = epoch % t_max
        return eta_min + (base_lr - eta_min) * (1 + np.cos(np.pi * t / t_max)) / 2

    return scheduler

t_37C_240min_dict = ppp.load(open('D:/uic/lab/data/pickle_files/ecoli_non_specific_search_poly_dict.p','rb'))
matrix,target = matrix_target(t_37C_240min_dict)
matrix = custom_ohe(matrix)
X_train, X_test, target_train, target_test = train_test_data_split(matrix,target)
# X_train, X_test = X_train.reshape(X_train.shape[0],31,int(X_train.shape[1]/31)), \
#                         X_test.reshape(X_test.shape[0], 31, int(X_test.shape[1]/31))

# convert numpy array data into tensor
X_train_tensor,target_train_tensor = torch.tensor(X_train,dtype=torch.float32), torch.tensor(target_train,dtype=torch.long)
X_test_tensor, target_test_tensor = torch.tensor(X_test,dtype=torch.float32), torch.tensor(target_test,dtype=torch.long)
X_train_tensor = torch.reshape(X_train_tensor,(X_train_tensor.shape[0],31,22)) # samples, timesteps, number of features
X_test_tensor = torch.reshape(X_test_tensor,(X_test_tensor.shape[0],31,22))


train_tensor = TensorDataset(X_train_tensor,target_train_tensor)
train_loader = DataLoader(dataset=train_tensor, batch_size=64, shuffle=True)
test_tensor = TensorDataset(X_test_tensor,target_test_tensor)
test_loader = DataLoader(dataset=test_tensor, batch_size=64, shuffle=True)


input_dim = 22  # length of one-hot encoding 31mer
hidden_dim = 16
layer_dim = 2
output_dim = 2
seq_dim = 128

lr = 0.0005
n_epochs = 100
iterations_per_epoch = len(train_loader)
best_acc = 0
patience, trials = 100, 0
history = []

model = LSTMClassifier(input_dim, hidden_dim, layer_dim, output_dim)
model = model.cuda()
criterion = nn.CrossEntropyLoss()
opt = torch.optim.RMSprop(model.parameters(), lr=lr)
sched = CyclicLR(opt, cosine(t_max=iterations_per_epoch * 2, eta_min=lr/100))

print('Start model training')

for epoch in range(1, n_epochs + 1):

    for i, (x_batch, y_batch) in enumerate(train_loader):
        model.train()
        x_batch = x_batch.cuda()
        y_batch = y_batch.cuda()
        sched.step()
        opt.zero_grad()
        out = model(x_batch)
        loss = criterion(out, y_batch)
        loss.backward()
        opt.step()

    model.eval()
    correct, total = 0, 0
    for x_val, y_val in test_loader :
        x_val, y_val = [t.cuda() for t in (x_val, y_val)]
        out = model(x_val)
        val_loss = criterion(out,y_val)
        preds = F.log_softmax(out, dim=1).argmax(dim=1)
        total += y_val.size(0)
        correct += (preds == y_val).sum().item()

    acc = correct / total
    history.append((epoch,loss.item(),val_loss.item(), acc))
    if epoch % 5 == 0:
        print(f'Epoch: {epoch:3d}. Loss: {loss.item():.4f}. Val_loss: {val_loss.item():.4f}. Acc.: {acc:2.2%}')

    if acc > best_acc:
        trials = 0
        best_acc = acc
        torch.save(model.state_dict(), 'best.pth')
        print(f'Epoch {epoch} best model saved with accuracy: {best_acc:2.2%}')
    else:
        trials += 1
        # if try 100 times and still no better accuracy, terminate the training
        if trials >= patience:
            print(f'Early stopping on epoch {epoch}')
            break

ppp.dump(history, open('pytorch_history.p', 'wb'))
# model = LSTMClassifier(input_dim, hidden_dim, layer_dim, output_dim)
# model.load_state_dict(torch.load('best.pth'))
# model.eval()
