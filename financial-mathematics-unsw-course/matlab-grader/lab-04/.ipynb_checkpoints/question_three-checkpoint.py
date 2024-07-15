"""
Computing DFT using FFT
Shahid Hussain
"""
m = 6
from numpy.random import randn, seed
seed(1)
ym = randn(2**m)
print(ym)