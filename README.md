# Numerical Solution of 2D Laplace Equation

This repository contains three Python implementations of numerical solvers for the 2D Laplace equation using the finite difference method (FDM). The three iterative methods implemented are:

- **Jacobi Method**
- **Gauss-Seidel Method**
- **Successive Over-Relaxation (SOR) Method**

## Problem Description
The 2D Laplace equation is given by:

\[ \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} = 0 \]

This equation is solved over a square plate (\(n \times n\)) with Dirichlet boundary conditions. A heated boundary is defined where:
- \(u = 1\) along a section of the left boundary and the bottom boundary.
- \(u = 0\) elsewhere initially.

The numerical solution iterates until the maximum error between consecutive iterations is below a specified tolerance (\(10^{-4}\)).

## Code Descriptions
### 1. Jacobi Method (`jacobi.py`)
The Jacobi method updates each grid point independently using the values from the previous iteration:
\[ u_{i,j}^{(k+1)} = \frac{1}{4} ( u_{i+1,j}^{(k)} + u_{i-1,j}^{(k)} + u_{i,j+1}^{(k)} + u_{i,j-1}^{(k)} ) \]

✅ **Advantages**:
- Simple and easy to parallelize.
- All updates are computed independently before applying them to the grid.

❌ **Disadvantages**:
- Convergence is slow compared to Gauss-Seidel and SOR.
- Requires additional memory for storing the updated values.

### 2. Gauss-Seidel Method (`gauss_seidel.py`)
The Gauss-Seidel method updates each grid point sequentially, using already updated values in the current iteration:
\[ u_{i,j}^{(k+1)} = \frac{1}{4} ( u_{i+1,j}^{(k)} + u_{i-1,j}^{(k+1)} + u_{i,j+1}^{(k)} + u_{i,j-1}^{(k+1)} ) \]

✅ **Advantages**:
- Faster convergence than Jacobi since updates propagate more quickly.
- Requires less memory as it updates in place.

❌ **Disadvantages**:
- Harder to parallelize due to dependencies between updates.
- Still relatively slow for large grids.

### 3. Successive Over-Relaxation (SOR) Method (`relaxation.py`)
The SOR method enhances Gauss-Seidel by introducing a relaxation factor \( \omega \) to accelerate convergence:
\[ u_{i,j}^{(k+1)} = \omega \left( \frac{1}{4} ( u_{i+1,j}^{(k)} + u_{i-1,j}^{(k+1)} + u_{i,j+1}^{(k)} + u_{i,j-1}^{(k+1)} ) \right) + (1 - \omega) u_{i,j}^{(k)} \]

✅ **Advantages**:
- Significantly faster than both Jacobi and Gauss-Seidel when \(\omega\) is optimally chosen.

❌ **Disadvantages**:
- Choosing an optimal relaxation factor \( \omega \) requires trial and error.
- Not easily parallelizable.

## Dependencies
- Python 3.x
- `matplotlib` (for visualization)

To install the required package:
```sh
pip install matplotlib
```

## Running the Codes
Each method is implemented in a separate script. To run a specific solver, use:
```sh
python jacobi.py
python gauss_seidel.py
python relaxation.py
```
During execution, the program will print "still running" every 300 iterations to indicate progress. The final solution is displayed using contour plots.

## Results
Each method generates a contour plot of the steady-state temperature distribution across the square plate. The number of iterations required varies:
- Jacobi (slowest convergence)
- Gauss-Seidel (faster than Jacobi)
- SOR (fastest when \( \omega \) is optimally chosen, default \( \omega = 1.5 \))

## Author
This project was created to demonstrate different iterative numerical methods for solving PDEs. For any questions, feel free to reach out!

---
**Note:** The relaxation factor \( \omega \) in SOR can be tuned for better performance. The default value of 1.5 is a commonly used heuristic but may need adjustment for specific problems.

