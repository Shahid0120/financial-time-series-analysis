"""

"""
import numpy as np
from numpy import array, linspace, ones, cos, sin, pi
import matplotlib.pyplot as plt

# a & b  
a0 = 1
a = np.array([1, -1, 0, 1/2])
b = np.array([0, 2, -1])


T = 2
n = 3

f = 1/T

# Fine grid of plot points
tlo = -2
thi = 6
tplt = linspace(tlo, thi, 801)

# Trignometric polynomial. Remeber a_k, b_k stores as a[k-1], b[k-1]

yplt = a0 * np.ones(tplt.size)

for k in range(1, n+1):
    theta = 2 * pi * k * tplt / T
    yplt = yplt + a[k-1] * cos(theta) + b[k-1] * sin(theta)


plt.figure(1)
plt.clf()
plt.plot(tplt, yplt)
plt.grid(True)
plt.xlabel('t')
plt.title(f'Trigonometric polynomial of degree {n}, period {T}')
