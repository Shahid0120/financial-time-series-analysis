function [c, dcds] = blackscholes(S, K, r, sigma, Tmt)
% [c, dcds] = blackscholes(S, K, r, sigma, Tmt)
% MATH3311/MATH5335: File = blackscholes.m
% Black and Scholes formula for the value of a call option
% and its derivative with respect to volatility sigma
% S   = underlying asset price
% K   = strike price
% r   = risk-free interest rate
% Tmt = time to maturity = T - t where T = expiry
% If sigma is a vector of volatilities, then both the
% call value and its derivatives are vectors of the same size.
%
% Uses normpdf and normcdf from Statistics toolbox.

s = sigma * sqrt(Tmt);

d1 = ( log(S/K) + ( r + sigma.^2/2)*(Tmt) ) ./ s;
d2 = d1 -  s;

% Use normpdf and normcdf from Statistics toolbox
c = S .* normcdf(d1) - K * exp(-r*Tmt) * normcdf(d2);

% Derivative of call vlaue w.r.t. volatility sigma
dcds = S .* normpdf(d1) * sqrt(Tmt);
