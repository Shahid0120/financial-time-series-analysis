import numpy as np 
alpha = 0 
beta = 1/2
N = 10 ** 6 
T = 10
num =  50000
U0 = 0 
U1 = np.divide(T, np.multiply(2, N)) + 0
U_n_minus_one = U0
U_n = U1
U_n_coefficient = -2 * (num ** 2)/ (T / N)**2 - num / (T/N) - 1
U_n_minus_one_coefficient = num**2 / (T/N)**2
U_n_one_coefficient = num**2/(T/N)**2 + num/(T/N)


iter = 
while iter < num + 1:
    U_n_plus_one = (- np.multiply(U_n_minus_one, U_n_minus_one_coefficient) - np.multiply(U_n,U_n_coefficient)) / U_n_one_coefficient
    U_n_minus_one = U_n
    U_n =  U_n_plus_one
    iter += 1 

print(iter)
print(U_n_plus_one)