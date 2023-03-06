"""Analyze the relationship between $L_1$ error and particle numbers of RMBC
method on the Dyson Browian Motion problem
"""


import numpy as np
import pandas as pd
import scipy.integrate as integrate
from math import sqrt, pi
import matplotlib.pyplot as plt

import os

def command_composer(program, args):
    command = program
    for key in args:
        if args[key]:
            command += f" --{key}={args[key]}"
        else:
            command +=f" --{key}"
    return command

def execute(num):
    command = command_composer(
        'target/release/dyson-browian',
        {
            "output": f"lab-dbn/dyson-browian-{str(num)}.csv",
            "particle-num": str(num),
            "iterations": "10000000",
            "step-time": "0.001",
            "distribution": None,
            "low": "-2.0",
            "high": "2.0",
            "interval-num": "100",
            "burnin-step": "3000000"
        }
    )
    os.system(command)
    data = pd.read_csv(f'lab-dbn/dyson-browian-{str(num)}.csv')
    return data.to_numpy()

def f(x):
    if x<-sqrt(2) or x>sqrt(2): return 0
    return 1/pi*sqrt(2-x**2)

num_list = [10, 20, 40, 80, 160, 320, 640]
error = []
for num in num_list:
    bins_num = 100
    bin_range = (-2, 2)
    h = (bin_range[1]-bin_range[0])/bins_num
    hist, bin_edges = np.histogram(np.zeros(1), bins = bins_num, range = bin_range, density=True)
    hist = np.array(execute(num), dtype="float64")
    hist/=hist.sum()*h
    expected_hist = []
    for i in range(bins_num):
        expected_hist.append(integrate.quad(f, bin_edges[i], bin_edges[i+1])[0]/h)
    sum = 0
    for i in range(bins_num):
#        print(abs(hist[i]-expected_hist[i])*h)
        sum += abs(hist[i]-expected_hist[i])*h
    error.append(sum)

plt.plot(num_list, error)
plt.scatter(num_list, error)
plt.xscale('log')
plt.yscale('log')
plt.savefig('lab-dbn/result6.png')
plt.plot()
