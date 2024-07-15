import numpy as np 
import scipy 
from numpy import loadtxt,  ones, empty, linspace, exp 
from numpy import zeros_like, pi, cos, sin, round, floor
from scipy.linalg import lstsq
from scipy.interpolate import InterpolatedUnivariateSpline, CubicSpline
from scipy.fft import fft
from matplotlib.pyplot import plot, figure, legend, grid, xlabel

data = loadtxt('gauss.dat')