#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Testing Symetric Matrices

@author: shahid
"""
import numpy as np 
from symchSol import symchk
from numpy import array


   
A = np.array([[1,2], [2,1]])
sent = symchk(A, -3)
print(A)
print(sent)
