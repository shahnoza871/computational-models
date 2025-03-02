import matplotlib.pyplot as plt
n=100
dx=dy=1/n
dt=0.001
itt=0

# introducing relaxation coefficient
w=1.5

u=[[]]

# initial state of u
for i in range(n):
    u_row = []
    for j in range(n):
        if j==0 and (n/2<=i<=n-1) or (i==n-1 and 0<=j<=n/3):
            u_row.append(1)
        else:
            u_row.append(0)
    u[0].append(u_row)

xlist=[i*dx for i in range(n)]
ylist=[j*dy for j in range(n)]

while True: 
    error = 0
    u.append([])
    itt+=1
    for i in range(n):
        u_row = []
        for j in range(n):
            if j==0 and (n/2<=i<=n-1) or (i==n-1 and 0<=j<=n/3):
                u_row.append(1)
            else:
                u_row.append(0)
        u[itt].append(u_row)
    for i in range(1, n-1):
        for j in range(1, n-1):
            u[itt][i][j]=(w/4)*(u[itt-1][i+1][j]+u[itt][i-1][j]+u[itt][i][j-1]+u[itt-1][i][j+1]-4*(1-1/w)*u[itt-1][i][j])
            error = max(error, abs(u[itt][i][j]-u[itt-1][i][j]))

    if itt%300==0:
        print('still running')
        plt.contourf(xlist, ylist, u[-1])
        plt.show()
    
    if error<0.0001:
        break

print(itt)
plt.contourf(xlist, ylist, u[-1])
plt.show()