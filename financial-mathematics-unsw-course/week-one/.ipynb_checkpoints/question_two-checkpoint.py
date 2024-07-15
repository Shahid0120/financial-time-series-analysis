import numpy as np 
import matplotlib.pyplot as plt 
from numpy import log,pi
from numpy import exp
from numpy import finfo, float32
from numpy import float64
from numpy import arange, sin
from numpy import finfo, log, sqrt
from numpy import empty, exp


eps = finfo(float).eps
a = 7 * log(pi)

# Question 2. 
print('solution ',  a)
print(pi ** 7 == exp(a))
print('difference: ', abs(pi ** 7 - exp(a)))
print('expected error : ' ,pi ** 7 * eps)

## What is epsilon ?
## Comapre epsilon_1 to epsilon_2 andtwhy epislon_1 is chosen 
## why is is epsilon_0 not used? the different between 0 to the next machien number? why is it not used?

## Question 3.
print('Q3.a') 
print(log(0))
print(log(-1))

print('Q3.b')
z = float64(0)
print(1/z)
print(-1/z)

print('Q3.c')
# print(0/0)
print(z/z)
print(np.inf - np.inf)
print(0 * np.inf)
print( 0 * np.nan) 
print(exp(-np.inf))

## Question 4. 
## overflow and underflow 

def f(x):
        return sin(pi * x)
x = arange(101)
y = f(x) 
plt.plot(x,y)


# Create an array of k values from 0 to 100
k = np.arange(0, 101)

# Calculate f(k) using sin(pi * k)
f_k = np.sin(pi * k)

# Print the results
print(f_k)

# Question 5 
# is real min the smallest number ? no underflow, but the smallest floapoint normalised number

## Question 4
# (a)
#epsilon_single = np.finfo(float).eps
#print(epsilon_single)

## Question 6
#J[empty] = 1 - exp(-1) 
#for i in range(1 ,11):
#   J[i] = 1 - i*J[i - 1] 
#print(J)
# print(exp(-1))

## Question 7 
# Speed (flops/sec) = Clock (cycles/sec) × Cores × flops/cycle

# a 2.5 cycles/nanseconds 

# b 

## Question 9 
# why is 2^3 - 1? where is the minus 1 come from? since we  0 0 0 doesnt have a FPN
