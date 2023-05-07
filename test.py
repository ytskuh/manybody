import torch
import time

a=time.time()
data=torch.randn((100, 100), device='cuda')
b=time.time()

print(b-a)

a=time.time()
data=torch.randn((100, 100), device='cpu')
b=time.time()
print(b-a)

a=time.time()
data=torch.randn((100, 100), device='cuda')
b=time.time()

print(b-a)