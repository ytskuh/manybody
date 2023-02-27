import numpy as np
import pandas as pd
import scipy.integrate as integrate
from math import sqrt, pi
import matplotlib.pyplot as plt

import os

def execute(num):
    command = f'target/release/dyson-browian --output lab-dbn/dyson-browian-{str(num)}.csv --particle-num {str(num)} --iterations 10000000 --step-time 0.0001'
#    os.system(command)
    data = pd.read_csv(f'lab-dbn/dyson-browian-{str(num)}.csv')
    return data.to_numpy()

num_list = [10, 20, 40, 80, 160, 320, 640]
error = []
for num in num_list:
    bins_num = 20
    bin_range = (-sqrt(2), sqrt(2))
    h = (bin_range[1]-bin_range[0])/bins_num
    hist, bin_edges = np.histogram(execute(num)[3000000:], bins = bins_num, range = bin_range, density=True)
    expected_hist = []
    for i in range(bins_num):
        expected_hist.append(integrate.quad(lambda x: 1/pi*sqrt(1-x**2/4), bin_edges[i], bin_edges[i+1])[0]/h)
    sum = 0
    for i in range(bins_num):
#        print(abs(hist[i]-expected_hist[i])*h)
        sum += abs(hist[i]-expected_hist[i])*h
    error.append(sum)

for i in range(len(error)-1):
    print(error[i+1]-error[i])

plt.plot(num_list, error)
plt.scatter(num_list, error)
plt.xscale('log')
plt.savefig('lab-dbn/result.png')
plt.plot()
