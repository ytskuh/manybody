#!/bin/python

from math import sqrt, pi
import numpy as np
import matplotlib.pyplot as plt
import polars as pl

import os
import sys

top = os.path.dirname(
    os.path.dirname(
        os.path.abspath(sys.argv[0])
    )
)
sys.path.append(top)

import wrapper

bin_range = (-2, 2)
bin_num = 100

def execute1(num):
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
    os.system(command)
    with open(filename, 'r') as f:
        data_str = f.read()
    data_list = data_str.strip().rstrip(']').lstrip('[').split(',')
    data = np.array(data_list, dtype=float)
    return data

def execute2(num):
    filename = f"lab-dbn/dyson-browian-mh-{str(num)}.txt"
    command = wrapper.command_composer(
        'target/release/dyson-browian-mh',
        {
            "output": filename,
            "particle-num": str(num),
            "iteration": str(80000*num),
        }
    )
    os.system(command)
    data = pl.read_csv(filename)
    return data.to_numpy()

def f(x):
    if x<-sqrt(2) or x>sqrt(2): return 0
    return 1/pi*sqrt(2-x**2)

_, bins = np.histogram(np.empty(0), density=True, bins=bin_num, range = bin_range)

counts = execute1(100)

x = execute2(100)
counts2, bins = np.histogram(x, density=True, bins=bin_num, range = bin_range)

plt.stairs(counts, bins)
plt.stairs(counts2, bins)
X = np.linspace(-np.sqrt(2) , np.sqrt(2), 100)
Y = np.sqrt(2-X**2)/np.pi
plt.plot(X, Y)
plt.xlim((-2, 2))
plt.show()
