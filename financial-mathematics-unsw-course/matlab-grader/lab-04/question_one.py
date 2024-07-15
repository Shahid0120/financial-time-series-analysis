"""
Lab 04 Question Newtons Method
Shahid Hussain
"""
import numpy as np 
import scipy
import matplotlib.pyplot as plt 
from numpy import exp, log
from scipy.optimize import fsolve

def main():
    dfx = []
    fx = []

    for index in range(0, 21):
        index = index / 10
        value_fz = function(index)
        value_dfz = derivativeFunction(index)
        fx.append(value_fz)
        dfx.append(value_dfz)

    z0 = 19 / 10
    zn, function_vals = newtonMethod(fx, z0, 20)
    zm, func_values = fsolveMethod(fx, z0, 20)
    diffnm = 1.487962065498177 - 1.4879620654981773
    print(f"the difference {diffnm} for diffnm for both methods")
    return 

def function(z):
    return exp(z) - 2 * (z ** 2)

def derivativeFunction(z):
    return exp(z) - 4 * z



def fsolveMethod(fx, z0, num):
    zm = []
    func_values = []
    while num > 0:
        root = fsolve(function, z0)
        print(root)
        zm.append(root[0])
        func_values.append(function(z0))
        z0 = root
        iter -= 1 

    return zm, func_values 

def newtonMethod(fx, z0, iter):
    zn = [] 
    function_vals = []
    while iter > 0:
        func_value = function(z0)
        deriv_value = derivativeFunction(z0)
        z_next = z0 - (func_value / deriv_value)
        zn.append(z_next)
        function_val_z_next =  function(z_next)
        function_vals.append(function_val_z_next)

        z0 = z_next
        iter -= 1 

    return zn, function_vals

if __name__ == "__main__":
    main()
