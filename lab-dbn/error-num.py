"""Analyze the relationship between $L_1$ error and particle numbers of RMBC
method on the Dyson Browian Motion problem
"""

import numpy as np
import pandas as pd
import scipy.integrate as integrate
from math import sqrt, pi
import matplotlib.pyplot as plt
import os
import sys

sys.path.append("/home/fred/projects/manybody")
import wrapper

bin_range = (-2, 2)
bin_num = 100

def execute(num):
    filename = f"lab-dbn/dyson-browian-{str(num)}.txt"
    command = wrapper.command_composer(
        'target/release/dyson-browian',
        {
            "output": filename,
            "particle-num": str(num),
            "iteration": str(80000*num),
            "step-time": "0.00001",
            "p": "4",
            "low": str(bin_range[0]),
            "high": str(bin_range[1]),
            "interval-num": str(bin_num),
            "burn-in": str(40000*num)
        }
    )
#    os.system(command)
    with open(filename, 'r') as f:
        data_str = f.read()
    data_list = data_str.strip().rstrip(']').lstrip('[').split(',')
    data = np.array(data_list, dtype=float)
    return data

def f(x):
    if x<-sqrt(2) or x>sqrt(2): return 0
    return 1/pi*sqrt(2-x**2)

num_list = [10, 20, 40, 80, 160, 320, 640]
error = []
for num in num_list:
    hist = execute(num)
    error.append(wrapper.l1divergence(hist, f, bin_range[0], bin_range[1], bin_num))

plt.plot(num_list, error)
plt.scatter(num_list, error)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number')
plt.ylabel(r'$L_1$-divergence')
plt.savefig('lab-dbn/result5.png')
plt.plot()

