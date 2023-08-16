#!/bin/python

from math import sqrt, pi
import numpy as np
import matplotlib.pyplot as plt
import polars as pl

import os
import sys
import time

top = os.path.dirname(
    os.path.dirname(
        os.path.abspath(sys.argv[0])
    )
)
sys.path.append(top)

import wrapper

bin_range = (-2, 2)
bin_num = 40

def execute1(num2):
    num=500
    filename = f"lab-dbn/dyson-browian-{str(num)}.txt"
    command = wrapper.command_composer(
        'target/release/dyson-browian',
        {
            "output": filename,
            "particle-num": str(num),
            "iteration": str(num2),
            "step-time": "0.0001",
            "p": "2",
            "raw": "",
            "burn-in": str(num2//2)
        }
    )
    a=time.time()
    os.system(command)
    b=time.time()
    data = pl.read_csv(filename)
    return (data.to_numpy(), b-a)

def execute2(num2):
    num=500
    filename = f"lab-dbn/dyson-browian-mh-{str(num)}.txt"
    command = wrapper.command_composer(
        'target/release/dyson-browian-mh',
        {
            "output": filename,
            "particle-num": str(num),
            "iteration": str(num2),
        }
    )
    a=time.time()
    os.system(command)
    b=time.time()
    data = pl.read_csv(filename)
    return (data.to_numpy(), b-a)

def f(x):
    if x<-sqrt(2) or x>sqrt(2): return 0
    return 1/pi*sqrt(2-x**2)

e1=[]
e2=[]
t1l=[]
t2l=[]
for ns in [10000, 20000, 40000, 80000, 160000, 320000, 640000, 1280000]:
    x, t1 = execute1(ns)
    counts, bins = np.histogram(x, density=True, bins=bin_num, range = bin_range)
    x, t2 = execute2(ns)
    counts2, bins = np.histogram(x, density=True, bins=bin_num, range = bin_range)

    e1.append(wrapper.l1divergence(counts, f, bin_range[0], bin_range[1], bin_num))
    e2.append(wrapper.l1divergence(counts2, f, bin_range[0], bin_range[1], bin_num))
    t1l.append(t1)
    t2l.append(t2)

print(e1)
print(e2)
print(t1l)
print(t2l)

plt.scatter(t1l, e1, marker='o', label='RBMC')
plt.plot(t1l, e1)
plt.scatter(t2l, e2, marker='^', label='MH')
plt.plot(t2l, e2)
plt.xlabel('time consumed')
plt.ylabel(r'$L_1$ divergence')
plt.yscale('log')
plt.show()

# X = np.linspace(-np.sqrt(2) , np.sqrt(2), 100)
# Y = np.sqrt(2-X**2)/np.pi
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
# ax1.stairs(counts, bins, label='RBMC')
# ax1.plot(X, Y, label='limit')
# ax1.set_title('RBMC vs limit')
# ax2.stairs(counts2, bins, label='MH')
# ax2.plot(X, Y, label='limit')
# ax2.set_title('MH vs limit')
# ax1.set_xlim((-2, 2))
# ax2.set_xlim((-2, 2))
# fig.savefig('testresult.png')
# fig.show()
