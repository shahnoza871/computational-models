# Poisson Equation Solver Using the Implicit Thomas Algorithm

## Overview
This project implements a numerical solution to the Poisson equation using the implicit Thomas algorithm. The solver achieves an error order of \(O(dx^2)\), ensuring a second-order accurate approximation. The computed numerical solution is compared against the analytical solution, and the error is evaluated.

## Problem Statement
The Poisson equation being solved is:
\[ \frac{d^2u}{dx^2} = f(x) \]
where:
\[ f(x) = -(x^2 + 2x + 3) \]
and the analytical solution is:
\[ u(x) = -\left( \frac{x^4}{12} + \frac{x^3}{3} + \frac{3x^2}{2} - \frac{35}{12} \right) \]

## Numerical Method
The solution is obtained using the **Thomas algorithm**, a specialized form of Gaussian elimination for tridiagonal systems. The discretization is performed using the **finite difference method**, where the second derivative is approximated as:
\[ \frac{d^2u}{dx^2} \approx \frac{u_{i-1} - 2u_i + u_{i+1}}{dx^2} \]
The resulting system of linear equations is then solved using forward and backward sweeps.

### Steps Involved:
1. Discretize the domain with \( n \) points and step size \( dx = \frac{1}{n} \).
2. Construct the coefficient matrix with values:
   - \( A = 1.0 \) (subdiagonal)
   - \( B = -2.0 \) (diagonal)
   - \( C = 1.0 \) (superdiagonal)
3. Use the Thomas algorithm to solve for \( u(x) \).
4. Compare the numerical and analytical solutions.
5. Compute the absolute error.
6. Visualize results using **matplotlib**.

## Results
- The numerical solution converges to the analytical solution with an error \( O(dx^2) \).
- The **maximum error** is displayed on the plot.
- The visualization compares both solutions for validation.

