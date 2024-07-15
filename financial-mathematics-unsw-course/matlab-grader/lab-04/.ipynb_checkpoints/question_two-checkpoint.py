"""
Lab04 Question 2 
Shahid Hussain
"""
import numpy as np 
import scipy 

def main():
    N = 10 ** 6 
    T = 10
    n = [1000, 30000, 60000, 90000]
    alpha = 0 
    beta = 1/2
    # sx calculation 
    sx = sXCalculation(N, T, n)

    # U1 calculation 
    U1 = np.divide(T, np.multiply(2, N)) + 0

    
    # vi
    num = 1000
    U1000 = functionU(N,T, num, 0, 1/2, U1, 0)
    #print(U1000)


    num = 10000
    U10000 = functionU(N,T, num, 0, 1/2, U1, 0)
    #print(U10000)

    num =  50000
    U50000 = functionU(N,T, num, 0, 1/2, U1, 0)
    print(U50000)

    num = 100000
    U100000 = functionU(N,T, num, 0, 1/2, U1, 0)
    #print(U100000)

    return 

def sXCalculation(N,T, n:list):
    sx = 0
    for val in n:
        sx += np.divide(np.multiply(T,val), N)

    return sx

def functionU(N,T, num, alpha, beta, U1, U0):

    U_n_minus_one = U0
    U_n = U1
    U_n_coefficient = -2 * (num ** 2)/ (T/N)**2 - num / (T/N) - 1
    U_n_minus_one_coefficient = num**2 / (T/N)**2
    U_n_one_coefficient = num**2/(T/N)**2 + num/(T/N)


    iter = 1
    while iter < num + 1:
        U_n_plus_one = (- np.multiply(U_n_minus_one, U_n_minus_one_coefficient) - np.multiply(U_n,U_n_coefficient)) / U_n_one_coefficient
        U_n_minus_one = U_n
        U_n =  U_n_plus_one
        iter += 1 

    return U_n_plus_one


if __name__ == "__main__":
    main() 
