import matplotlib.pyplot as plt
n=100
dx=dy=1/n # nxn square plate that is heated
dt=0.001
itt=0

u=[]

# initial state of u
for i in range(n):
    u_row = []
    for j in range(n):
        if j==0 and (n/2<=i<=n-1) or (i==n-1 and 0<=j<=n/3):
            u_row.append(1)
        else:
            u_row.append(0)
    u.append(u_row)

xlist=[i*dx for i in range(n)]
ylist=[j*dy for j in range(n)]

while True: 
    # iterations 
    un = []
    error = 0
    
    # initial state of u for each iteration itt
    for i in range(n):
        u_row = []
        for j in range(n):
            if j==0 and (n/2<=i<=n-1) or (i==n-1 and 0<=j<=n/3):
                u_row.append(1)
            else:
                u_row.append(0)
        un.append(u_row)
    for i in range(1, n-1):
        for j in range(1, n-1):
            un[i][j]=(1/4)*(u[i+1][j]+u[i-1][j]+u[i][j-1]+u[i][j+1])
            error = max(error, abs(un[i][j]-u[i][j]))
    itt+=1

    # to iterate repeatedly using prev data in order to find more precise u
    u=un

    if itt%300==0:
        print('still running')
        plt.contourf(xlist, ylist, u)
        plt.show()
    
    if error<0.0001:
        break

print(itt)
plt.contourf(xlist, ylist, u)
plt.show()