from scipy.optimize import linprog
from numpy import zeros, ones, eye, inf, concatenate

def l1lp(A, b):
    """
    l1lp(A, b)

    MATH3311/MATH5335: file = l1lp.py
    Find the best l_1 solution to A*x = b by

    Minimize || A*x - b ||_1

    Converts l_1 problem into a linear programming problem
    and uses Python's linprog routine to solve the LP

    Input required
    A  : m by n coefficient matrix
    b  : m vector of data values

    Output
    x  : n vector of best l_1 estimator
    r  : m vector of residuals r = A*x - b
    f  : || r ||_1
    """
    m, n = A.shape
    c = zeros(n+m)
    c[n:] = ones(m)
    A_ub = zeros((2*m, n+m))
    A_ub[0:m,:] = concatenate((A.copy(), -eye(m)),axis=1)
    A_ub[m:,:] = concatenate((-A.copy(),-eye(m)),axis=1)
    b_ub =concatenate((b.copy(),-b.copy()),axis=None)
    bounds = []
    for j in range(n):
        bounds.append((-inf,inf))
    for j in range(m):
        bounds.append((0.0,inf))

    A_eq = None
    b_eq = None
    res = linprog(c, A_ub, b_ub, A_eq, b_eq, bounds)

    x = res.x[:n]
    r = A @ x - b
    f = abs(r).sum()
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

#    x, r, f, res = l1lpy(A, y)
#    print(res.message)

#   # Plot data and l1 fit
#    figure(1)
#    plot(t, y, 'o', t[[0,-1]], x[0] + x[1] * t[[0,-1]])
#    grid(True)
#    title('L1 fit')
#    legend(('Data', f'{x[0]:0.3f} + {x[1]:0.3f} t'), loc = 'upper left')
#    show()