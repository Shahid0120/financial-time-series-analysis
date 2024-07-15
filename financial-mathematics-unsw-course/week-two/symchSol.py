#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Checking if Matrix is Symetric

@author: shahid
"""
import numpy as np 
from numpy import all,any, maximum, finfo, ndarray, float32, float64, transpose, array, abs

def symchk(A, tol=None):
    # check if tol is reference
    if tol is None:
        tol = 10 * len(A)**2 * finfo(np.float64).eps
    elif tol < 0:
        raise ValueError("tolerence is a negative number")
    
    # Checks Square
    if A.shape[0] != A.shape[1]:
        return "This is not a square matrix"
    
    # Checks if it is a Square matrix 
    n = A.shape[0]
    for i in range(n):
        for j in range(n):
            if abs(A[i, j] - A[j, i]) > tol * np.maximum(A[i, j], A[j, i]):
                return "This is not a symmetric but square matrix"
            
    return "This is a symmetric and square matrix"

