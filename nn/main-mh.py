#%%
import torch
from torch import nn
import matplotlib.pyplot as plt
import random

dev = "cuda"

DIM = 1

class nnm(nn.Module):
    def __init__(self, N):
        super().__init__()
        self.width = N
        self.hidden_layer = nn.Linear(DIM, N)
        self.output_layer = nn.Linear(N, 1, bias=False)

    def forward(self, x):
        hidden = torch.sigmoid(self.hidden_layer(x))
        output = self.output_layer(hidden)
        return output/self.width
    
    def activated(self, x):
        return self.output_layer.weight.data.reshape(-1)*torch.sigmoid(self.hidden_layer(x))
    
WIDTH = 64

net = nnm(WIDTH)
net.to(dev)

#%%
NUM_DATA = 100000
DELTA = 0.5
y_train = (torch.rand(NUM_DATA, device=dev).round()*2-1).unsqueeze(1)
x_train = torch.randn((NUM_DATA, DIM), device=dev)*(1+DELTA)*(y_train==1) + torch.randn((NUM_DATA, DIM), device=dev)*(1-DELTA)*(y_train==-1)

# %%
a = 1
beta = 1
hidden_layer_weight = []
hidden_layer_bias = []
output_layer_weight = []
BURN_IN = 50000

for i in range(NUM_DATA-20):
    x = x_train[i]
    y = y_train[i:i+20]
    ax = net.activated(x).data.clone()

    rand1 = torch.randn(WIDTH, DIM, device=dev)*a
    rand2 = torch.randn(WIDTH, device=dev)*a
    rand3 = torch.randn(1, WIDTH, device=dev)*a

    net.hidden_layer.weight.data += rand1
    net.hidden_layer.bias.data += rand2
    net.output_layer.weight.data += rand3
    ax_new = net.activated(x).data

    alpha=torch.exp(-beta * (-2/WIDTH*(y.reshape(-1)*(ax_new.sum(axis=-1)-ax.sum(axis=-1))).mean() + 1/(2*WIDTH**2)*(ax_new.pow(2).sum(axis=-1)-ax.pow(2).sum(axis=-1)).mean() + 1/(2*WIDTH**2)*((ax_new.reshape(20, 1, -1)*ax_new.reshape(20, -1, 1)).sum(axis=(-1, -2))-(ax.reshape(20, 1, -1)*ax.reshape(20, -1, 1)).sum(axis=(-1, -2))).mean()))

    t = random.random()
    if t > alpha:
        net.hidden_layer.weight.data -= rand1
        net.hidden_layer.bias.data -= rand2
        net.output_layer.weight.data -= rand3

    if i > BURN_IN and i%100 ==0:
        hidden_layer_weight.append(net.hidden_layer.weight.data.clone())
        hidden_layer_bias.append(net.hidden_layer.bias.data.clone())
        output_layer_weight.append(net.output_layer.weight.data.clone())

FINAL_WIDTH = WIDTH * len(hidden_layer_bias)
net_mean = nnm(FINAL_WIDTH)
net_mean.to(dev)
net_mean.hidden_layer.weight.data = torch.cat(hidden_layer_weight, dim=0)
net_mean.hidden_layer.bias.data = torch.cat(hidden_layer_bias, dim=0)
net_mean.output_layer.weight.data = torch.cat(output_layer_weight, dim=1)
# %%
x_test = torch.linspace(-6, 6, 100, device=dev).unsqueeze(1)
y_pred = net_mean(x_test)
plt.scatter(x_train.cpu(), y_train.cpu())
plt.plot(x_test.cpu(), y_pred.detach().cpu())
plt.show()
# %%
