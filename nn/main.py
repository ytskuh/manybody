#%%
import torch
from torch import nn
import matplotlib.pyplot as plt

device = "cuda"

NEURON_NUM = 64
TRAIN_DATA = 16

class NeuralNetwork(nn.Module):
    def __init__(self, N):
        super().__init__()
        self.width = N
        self.hidden_layer = nn.Linear(1, N)
        self.output_layer = nn.Linear(N, 1, bias=False)

    def forward(self, x):
        hidden = torch.sigmoid(self.hidden_layer(x))
        output = self.output_layer(hidden)
        return output/self.width
    
net = NeuralNetwork(NEURON_NUM)
net.to(device)

x_train = torch.rand(TRAIN_DATA, 1, device=device)
y_train = torch.sin(3*x_train) + torch.randn(TRAIN_DATA, 1, device=device)*0.1

criterion = nn.MSELoss()
optimizer = torch.optim.SGD(net.parameters(), lr=10)


BURN_IN = 100000
TOTAL_EPOCH = 200000

#%%
a=0.01
hidden_layer_weight = []
hidden_layer_bias = []
output_layer_weight = []
for epoch in range(TOTAL_EPOCH):
    y_pred = net(x_train)
    loss = criterion(y_pred, y_train)

    optimizer.zero_grad()
    loss.backward()
    optimizer.step()

    net.hidden_layer.weight.data += torch.randn(NEURON_NUM, 1, device=device)*a
    net.hidden_layer.bias.data += torch.randn(NEURON_NUM, device=device)*a
    net.output_layer.weight.data += torch.randn(1, NEURON_NUM, device=device)*a

    if epoch > BURN_IN and epoch%100 ==0:
        hidden_layer_weight.append(net.hidden_layer.weight.data.clone())
        hidden_layer_bias.append(net.hidden_layer.bias.data.clone())
        output_layer_weight.append(net.output_layer.weight.data.clone())
    
FINAL_WIDTH = NEURON_NUM * len(hidden_layer_bias)
net_mean = NeuralNetwork(FINAL_WIDTH)
net_mean.to(device)
net_mean.hidden_layer.weight.data = torch.cat(hidden_layer_weight, dim=0)
net_mean.hidden_layer.bias.data = torch.cat(hidden_layer_bias, dim=0)
net_mean.output_layer.weight.data = torch.cat(output_layer_weight, dim=1)


#%%
x_test = torch.linspace(0, 1, 100, device=device).unsqueeze(1)
y_pred = net_mean(x_test)
plt.scatter(x_train.cpu(), y_train.cpu())
plt.plot(x_test.cpu(), y_pred.detach().cpu())
plt.show()

#%%
def plot(i):
    net = NeuralNetwork(NEURON_NUM)
    net.hidden_layer.weight.data = hidden_layer_weight[i]
    net.hidden_layer.bias.data = hidden_layer_bias[i]
    net.output_layer.weight.data = output_layer_weight[i]
    x_test = torch.linspace(0, 1, 100, device=device).unsqueeze(1)
    y_pred = net(x_test)
    plt.scatter(x_train.cpu(), y_train.cpu())
    plt.plot(x_test.cpu(), y_pred.detach().cpu())
    plt.show()
    return

# %%
