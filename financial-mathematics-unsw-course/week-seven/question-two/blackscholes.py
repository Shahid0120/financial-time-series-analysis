# MATH3311/MATH5335: File = blackscholes.py

from numpy import sqrt, log, exp
from scipy.stats import norm

def blackscholes(S, K, r, sigma, Tmt):
    """
    Black-Scholes formula for the value of a call option.

    Also returns the rate of change of the option value with respect 
    to the volatility sigma.  
    Input arguments:

    S   = underlying asset price
    K   = strike price
    r   = risk-free interest rate
    Tmt = time to maturity = T - t where T = expiry

    If sigma is a vector of volatilities, then both the
    call value and its derivatives are vectors of the same size."""

    s = sigma * sqrt(Tmt)
    d1 = ( log(S/K) + ( r + sigma**2/2) * Tmt ) / s
    d2 = d1 -  s

    c = S * norm.cdf(d1) - K * exp(-r*Tmt) * norm.cdf(d2)
    # Derivative of call vlaue w.r.t. volatility sigma
    dcds = S * norm.pdf(d1) * sqrt(Tmt)

    return c, dcds

def bsput(S, X, r, sigma, Tmt):
    """
    Black and Scholes formula for the value of a European put option
    and its derivative with respect to volatility sigma
    If sigma is a vector of volatilities, then both the
    call vlaue in its derivatives are vectors of the same size."""

    s = sigma * sqrt(Tmt)

    d1 = ( log(S/X) + ( r + sigma**2 / 2 )*Tmt ) / s
    d2 = d1 -  s

    c = S * norm.cdf(d1) - X * exp(-r*Tmt) * norm.cdf(d2)
    p = X * exp(-r*Tmt) * norm.cdf(-d2) - S * norm.cdf(-d1)

    # From put-call parity c + X*exp(-r*Tmt) = p + S
    # we have p'(sigma) = c'(sigma)
    dpds = S * norm.pdf(d1) * sqrt(Tmt)

    return p, dpds