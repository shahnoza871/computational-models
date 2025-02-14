import numpy as np
import matplotlib.pyplot as plt


def analytical(x):
    return -((x**4)/12+(x**3)/3+(3*x**2)/2-35/12)

def f(x):
    return (-1)*(x**2+2*x+3)


def D(i, dx):
    return dx**2*f(i*dx)


n=100
dx=1/n
xlist = [i*dx for i in range(n)]

A = 1.0
B = -2.0
C = A
alphas = [1]
# alphas = [0]
betas = [0]
for i in range(n):
    alphas.append(-A/(B+C*alphas[i]))
    betas.append((D(i, dx) - C * betas[i]) / (B + C * alphas[i]))

P = np.zeros(n)
# P[0] = (1+(np.pi**3-np.pi*(3+2*np.e))/(np.pi**2+1)**2)+(-1)*(2*np.pi/(1+np.pi**2)**2)
P[n-1]=1
for i in range(n-2, -1, -1):
    P[i]=alphas[i+1]*P[i+1]+betas[i+1]

analytical_solution = []
for i in range(len(xlist)):
    analytical_solution.append(analytical(i*dx))

errors = []
for i in range(len(xlist)):
    errors.append(abs(analytical_solution[i]-P[i]))
# print(errors)

plt.plot(xlist, analytical_solution, label="Analytical Solution")
plt.plot(xlist, P, label="Numerical Solution", color='purple')
plt.text(0.1, 0.01, f"Max Error={max(errors)*100}%" )
plt.grid(True)
plt.legend()
plt.show()