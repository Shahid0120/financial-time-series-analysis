# Question One 

import numpy as np 
import scipy
import matplotlib.pyplot as plt
import datetime
from datetime import datetime, date
from scipy.optimize import root_scalar
from numpy import finfo, float64, sqrt, log, exp
from scipy.stats import norm


def blackScholesModel(S, K, r, sigma, Tmt):

    # Inputs for Black-Scholes Equation
    d1 = ((np.log(S/K) + r*np.sqrt(Tmt))/ sigma*np.sqrt(Tmt)) + sigma*np.sqrt(Tmt)
    d2 = sigma*np.sqrt(Tmt)

    # Black-Scholes Equation
    c = S*norm.cdf(d1) + np.exp(-r*Tmt)*scipy.stats.norm.cdf(d2)

    # Derivative of Black-Scholes w.r.t volatility
    deriv = S* scipy.stats.norm.pdf(d1) * np.sqrt(Tmt)


    return c, deriv



# Inputs for Black-Scholes Model
S = 34.810
r = 0.025 
t0 = date.fromisoformat('2014-08-27')
tf = date.fromisoformat('2014-09-25')
Tmt = (tf-t0).days / 365

# Vector of strikes
K = [33.5, 34, 34.5, 35, 35.5, 36]
# Corresponding vector of call prices
cmkt = [1.460, 1.040, 0.640, 0.370, 0.155, 0.070]

siglo = finfo(float64).eps
sighi = 1.0
sigint = [siglo, sighi]

def blackScholesValue(sigma, K):
    # Inputs for Black-Scholes Model
    S = 34.810
    r = 0.025 
    t0 = date.fromisoformat('2014-08-27')
    tf = date.fromisoformat('2014-09-25')
    Tmt = (tf-t0).days / 365

    f, deriv = blackScholesModel(S, K, r, sigma, Tmt)
    return f 

# Iterate over all data stirck and calls
for j, Kj in enumerate(K):
    # function for where sigma is variable unknown
    # making function = 0, f = Call_Value(BS model) - Market_Call_price
    fj = lambda sigma: blackScholesValue(sigma, Kj) - cmkt[j]
    
    # Finding Implied Volatility 
    soln = root_scalar(fj, bracket=sigint, method='brentq')
    implied_volatility = soln.root
    
    print(f"strike:{Kj}    Market Call: {cmkt[j]}   NumberOfIterations: {soln.iterations}    Implied Volatility:{implied_volatility}")
