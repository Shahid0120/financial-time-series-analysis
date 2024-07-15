import numpy as np
import scipy


def main():
    N = 10**6  # Number of grid points
    T = 10      # Time interval
    n_values = [1000, 30000, 60000, 90000]  # List of steps for calculations
    alpha = 0   # Boundary parameter (assumed constant here)
    beta = 1/2   # Boundary parameter (assumed constant here)

    # Initial value calculation (replace with your actual calculation)
    U1 = np.divide(T, np.multiply(2, N))  # Example initial value

    # Use vectorized solution for efficiency
    U = functionURetry(N, T, n_values, alpha, beta, U1, 0)

    # Calculate and print the infinity norm of the error
    q = []
    for val in range(0, N + 1):
        q.append(np.divide(np.multiply(T,val), N)) # Generate x values
    bessel_values = scipy.special.jv(1, q)  # Compute Bessel function values
    error_norm = np.linalg.norm(U - bessel_values, np.inf)  # Infinity norm of the difference

    print(f"Error norm (Linf): {error_norm}")

    return


def functionURetry(Nsol, T, num, alpha, beta, U1, U0):
    """
    Solves the differential equation using a finite difference method (vectorized).
    """
    Usol = np.zeros(Nsol)
    Usol[0] = U0
    Usol[1] = U1
    x = num
    H = T / Nsol
    p = x**2 / H**2
    q = (-2* x**2 - x*H + H**2 * x**2 - H**2) / H**2
    r = (x**2 + H*x) / H**2
    Usol[2:] = ((-Usol[:-2] * p) - Usol[1:-1] * q) / r

    return Usol


if __name__ == "__main__":
    main()
