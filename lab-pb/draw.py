#!/bin/python

from math import sqrt, pi
import numpy as np
import matplotlib.pyplot as plt

import subprocess
import os
import sys

top = os.path.dirname(
    os.path.dirname(
        os.path.abspath(sys.argv[0])
    )
)
sys.path.append(top)
bin_range=(0.0, 20.0)
bin_num = 40

result = subprocess.run([
    "target/release/possion-boltzmann",
    "--low", str(bin_range[0]),
    "--high", str(bin_range[1]),
    "--interval-num", str(bin_num)
    ], 
    stdout=subprocess.PIPE)

output = result.stdout.decode('utf-8')
print(output)

str1, str2, _, _=output.split('\n')
c=np.array(eval(str1))
d=np.array(eval(str2))

bins = np.linspace(bin_range[0], bin_range[1], bin_num+1)

c /= 4*pi/3*(bins[1:]**3 - bins[:-1]**3)
d /= 4*pi/3*(bins[1:]**3 - bins[:-1]**3)
plt.stairs(c/2, bins)
plt.stairs(d, bins)
plt.show()
