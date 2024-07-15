from scipy.optimize import linprog
from numpy import zeros, empty, ones

def lilp(A, b):
    """MATH3311/MATH5335: File = lilp.py
    Find the best l_\infty solution to A*x = b by

    Minimize || A*x - b ||_\infty

    Converts l_\infty problem into a linear programming problem
    and uses Python's linprog routine to solve the LP

    Input required
    A  : m by n coefficient matrix
    b  : m vector of data values

    Output
    x  : n vector of best l_\infty estimator
    r  : m vector of residuals r = A*x - b
    f  : || r ||\infty_"""

    m, n = A.shape
    c = zeros(n+1)
    c[n] = 1.0
    A_ub = -ones((2*m,n+1))
    A_ub[:m,:n] = -A.copy()
    A_ub[m:,:n] =  A.copy()
    b_ub = empty(2*m)
    b_ub[:m] = -b.copy()
    b_ub[m:] =  b.copy()
    A_eq = None
    b_eq = None

    res = linprog(c, A_ub, b_ub, A_eq, b_eq)
    x = res.x[:n]
    r = A @ x - b
    f = abs(r).max()
    return x, r, f, res

#if __name__ == '__main__':
#    from numpy import linspace
#    from numpy.random import randn
#    from matplotlib.pyplot import figure, plot, grid, title, legend, show

#    m = 30
#    t = linspace(0, 3, m)
#    y = 1 + 0.5 * t + 0.2 * randn(m)
#    A = ones((m, 2))
#    A[:,1] = t

#    x, r, f, res = lilpy(A, y)
#    print(res.message)

#    # Plot data and l1 fit
#    figure(1)
#    plot(t, y, 'o', t[[0,-1]], x[0] + x[1] * t[[0,-1]])
#    grid(True)
#    title('L-infinity fit')
#    legend(('Data', f'{x[0]:0.3f} + {x[1]:0.3f} t'), loc = 'upper left')
#    show()