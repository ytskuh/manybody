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
        return self.output_layer.weight.data.reshape(-1)*torch.sigmoid(self.hidden_layer(x)).data
    
WIDTH = 64

net = nnm(WIDTH)
net.to(dev)

#%%
NUM_DATA = 1000000
DELTA = 0.5
x_train = torch.rand(NUM_DATA, device=dev).unsqueeze(1)
y_train = torch.sin(3*x_train)+torch.randn(NUM_DATA, device=dev).unsqueeze(1)*0.2

# %%
a = 0.1
beta = 1
hidden_layer_weight = []
hidden_layer_bias = []
output_layer_weight = []
BURN_IN = 990000

for i in range(NUM_DATA):
    x = x_train[i]
    y = y_train[i]
    ax = net.activated(x).clone()

    rand1 = torch.randn(WIDTH, DIM, device=dev)*a
    rand2 = torch.randn(WIDTH, device=dev)*a
    rand3 = torch.randn(1, WIDTH, device=dev)*a

    net.hidden_layer.weight.data += rand1
    net.hidden_layer.bias.data += rand2
    net.output_layer.weight.data += rand3
    ax_new = net.activated(x)

    alpha=torch.exp(-beta * (-2/WIDTH*y*(ax_new.sum()-ax.sum())+1/(2*WIDTH**2)*(ax_new.pow(2).sum()-ax.pow(2).sum())+1/(2*WIDTH**2)*((ax_new*ax_new.reshape(-1, 1)).sum()-(ax*ax.reshape(-1, 1)).sum())))

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
x_test = torch.linspace(-0, 1, 100, device=dev).unsqueeze(1)
y_pred = net_mean(x_test)
plt.scatter(x_train.cpu()[0:300], y_train.cpu()[0:300], s=1)
plt.plot(x_test.cpu(), y_pred.detach().cpu())
plt.show()
# %%
