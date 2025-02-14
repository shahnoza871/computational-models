import numpy as np
import matplotlib.pyplot as plt

# fourier series expansion 
def analytical(k, x, t):
    res = []
    for i in range(0, k):
        res.append(((16*(-1)**i)/((2*i+1)*np.pi))*(np.cos((2*i+1)*x/4))*np.exp(-t*(2*i+1)**2))
    return sum(res)

# numerical explicit scheme
def numerical(u, dx, dt):
    N  = len(u)
    for n in range(N-1):
        for i in range(1,N-1):
            u[n+1][i]=(1-32*dt/dx**2)*u[n][i]+(16*dt/dx**2)*(u[n][i+1]+u[n][i-1])
    return u


L = 2*np.pi
N = 100
dx = L/N
dt = dx**2/32

u = np.zeros((N, N))
# bc u(2pi, t)=0 
# bc u'(0, t)=0 => u[1][t]=u[0][t] hence i make assumption that u[0][t]=4
for i in range(N):
    u[i][0]=4
# ic u(0<=x<=2*pi, 0)=4
for i in range(int(N)):
    u[0][i]=4
u = numerical(u, dx, dt)

# plotting the graph
x=[]
t_ind = 99
analyt = []
for i in range(N):
    x.append(i*dx)
    analyt.append(analytical(k=50, x=i*dx, t=dt*t_ind))

# errors = []
# for i in range(N):
#     errors.append(abs(analyt[i]-u[t_ind+1][i]))

plt.plot(x, analyt, color = 'blue')
plt.plot(x, u[t_ind], color='green')
plt.xlabel("Space Coordinates")
plt.ylabel(f"Values of Heat Function at time {dt*t_ind}")
# plt.text(0.1, 0.01, f"Max Error={np.max(errors)*100}%" )
plt.grid()
plt.show()
