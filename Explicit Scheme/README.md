# Explicit Scheme for Solving the Heat Equation

## Overview
This project implements an explicit finite difference scheme to numerically solve the one-dimensional heat equation:

\[ \frac{\partial u}{\partial t} = \frac{\partial^2 u}{\partial x^2} \]

The numerical solution is compared with an analytical solution obtained using the Fourier series expansion.

## Methodology

### Analytical Solution
The analytical solution is constructed using the Fourier series expansion:
\[
    u(x,t) = \sum_{i=0}^{k} \frac{16(-1)^i}{(2i+1)\pi} \cos\left(\frac{(2i+1)x}{4}\right) e^{-t(2i+1)^2}
\]
where \( k \) controls the number of terms in the expansion.

### Numerical Solution
The explicit finite difference scheme is applied to discretize the equation:
\[
    u_i^{n+1} = \left(1 - 32 \frac{\Delta t}{\Delta x^2}\right) u_i^n + 16 \frac{\Delta t}{\Delta x^2} (u_{i+1}^n + u_{i-1}^n)
\]
where:
- \( \Delta x = L/N \) is the spatial step size,
- \( \Delta t = \Delta x^2 / 32 \) is the time step size (chosen for numerical stability),
- \( u_i^n \) represents the function value at position \( i \) and time step \( n \).

### Boundary and Initial Conditions
- **Boundary conditions:**
  - \( u(2\pi, t) = 0 \)
  - Neumann boundary condition \( u'(0,t) = 0 \), leading to \( u[1][t] = u[0][t] \)
- **Initial condition:** \( u(x, 0) = 4 \) for all \( 0 \leq x \leq 2\pi \).

## Implementation
- The analytical solution is computed using the Fourier series expansion.
- The numerical scheme iterates over discrete time steps to approximate the solution.
- The results are visualized by plotting the analytical and numerical solutions for a selected time step.

## Results
- The numerical solution is compared with the analytical solution at a specific time step.
- The plots illustrate the accuracy of the explicit scheme in solving the heat equation.
- The explicit method is efficient but requires a small time step for numerical stability.

## Future Improvements
- Implement an implicit method (e.g., Crank-Nicholson) for improved stability.
- Compute and visualize the error between numerical and analytical solutions.
- Extend to higher dimensions or different boundary conditions.

## Dependencies
- Python
- NumPy
- Matplotlib


