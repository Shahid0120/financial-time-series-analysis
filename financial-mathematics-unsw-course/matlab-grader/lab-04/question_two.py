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
    #num = 1000
    #U1000 = functionU(N,T, num, 0, 1/2, U1, 0)
    #print(U1000)

    #U10000 = functionU(N,T,10000, 0, 1/2, 5e-06, 0)
    #print(U10000)


    #U100000 = functionU(N,T, 100000, 0, 1/2, U1, 0)
    #print(U100000)

    #U50000 = functionU(N,T, 50000, 0, 1/2, U1, 0)
    #print(U50000)

    #num = 100000
    #U100000 = functionU(N,T, num, 0, 1/2, U1, 0)
    #print(U100000)

    x = []
    for val in range(0, N + 1):
        x.append(np.divide(np.multiply(T,val), N))

    U = compute_U(N, T, 0, U1)
    print(U[100000])
    bessel_values = scipy.special.jv(1, x)  
    normbessel1 = np.linalg.norm(U - bessel_values, np.inf) 
    R5normbessel1 = round(normbessel1, 5)
    
    print(f"normbessel1: {normbessel1}")
    print(f"R5normbessel1: {R5normbessel1}")
    return 

def sXCalculation(N,T, n:list):
    sx = 0
    for val in n:
        sx += np.divide(np.multiply(T,val), N)

    return sx

def compute_U(N, T, alpha, beta):
    Usol = np.zeros(N + 1, dtype=np.float64)
    Usol[0] = alpha
    Usol[1] = beta

    delta = T / N

    for n in range(1, N):
        U_n_minus_one = Usol[n - 1]
        U_n = Usol[n]
        U_n_coefficient = -2 * (n ** 2) / (delta ** 2) - n / delta - 1
        U_n_minus_one_coefficient = (n ** 2) / (delta ** 2)
        U_n_one_coefficient = (n ** 2) / (delta ** 2) + n / delta

        if U_n_one_coefficient == 0:
            # Handle division by zero or very small division
            U_n_plus_one = np.inf  # or any suitable value
        else:
            U_n_plus_one = (- U_n_minus_one * U_n_minus_one_coefficient - U_n * U_n_coefficient) / U_n_one_coefficient

        Usol[n + 1] = U_n_plus_one

    return Usol


def functionURetry(Nsol,T, num, alpha, beta, U1, U0):
    Usol = np.zeros(Nsol)
    Usol[0] = U0
    Usol[1] = U1
    x = num 
    H = T / Nsol
    p = x**2/ H**2
    q = (-2* x**2 - x*H + H**2 * x**2 - H**2) / H**2
    r = (x**2 + H*x)/ H**2
    for nsol in range(2, Nsol+1):
        Usol[nsol] = ((-Usol[nsol - 2] * p) - Usol[nsol - 1] * q) / r

    return Usol

if __name__ == "__main__":
    main() 

"""
% (a)
Nsol = 10^2;
% (b)
Usol = zeros(Nsol+1,1);
Usol(2) = 2/Nsol;
for nsol = 2:Nsol
    Usol(nsol+1) = -Usol(nsol-1)+Usol(nsol)*(2-Nsol^(-2)) + (2 + nsol/Nsol + nsol^2/Nsol^2)/Nsol^2;
end

% Solution
xsol = linspace(0,1,Nsol+1);
usol = sin(xsol) + xsol + xsol.^2;

% Error
Uerrsol = norm(usol-Usol');

%
Nlargesol = 10^6;
Ulargesol = zeros(Nlargesol+1,1);
Ulargesol(2) = 2/Nlargesol;
for nsol = 2:Nlargesol
    Ulargesol(nsol+1) = -Ulargesol(nsol-1)+Ulargesol(nsol)*(2-Nlargesol^(-2)) + (2 + nsol/Nlargesol + nsol^2/Nlargesol^2)/Nlargesol^2;
end
Ulargesol=Ulargesol(:);

"""
