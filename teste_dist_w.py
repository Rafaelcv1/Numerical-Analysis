import numpy as np
import math
import matplotlib.pyplot as plt

####################
#Gauss-seidel com SOR
####################

def genMatrix(n):
    output = np.zeros((n,n))
    output[0][0] = 2
    output[0][1]= -1
    for i in range(1, n-1):
        output[i][i-1]= -1
        output[i][i]=2
        output[i][i+1]= -1
    output[n-1][n-1] = 2
    output[n-1][n-2]= -1
    return output

def w_opt(n):
    A = genMatrix(n)
    D = np.diag(np.diag(A))
    LpU = A-D
    iter_jac = -np.linalg.inv(D)@LpU
    raio_spec = np.max(np.abs(np.linalg.eigvals(iter_jac)))
    out = 2/(1+math.sqrt(1-raio_spec**2))
    return out

def gauss_seidel_SOR(n, w, calcOpt):

    if calcOpt:
        A = genMatrix(n)
        D = np.diag(np.diag(A))
        LpU = A-D
        iter_jac = -np.linalg.inv(D)@LpU
        raio_spec = np.max(np.abs(np.linalg.eigvals(iter_jac)))
        w = 2/(1+math.sqrt(1-raio_spec**2))

    x0 = np.zeros(n)
    x = np.zeros(n)
    b = np.zeros(n)
    b[0] = 1
    b[n-1] = 1

    it = 0
    itmax = 10000
    tol = 1.0e-8
    dif = 1.0

    while (dif > tol and it < itmax):
        x[0] = (1-w)*x0[0] + w*(b[0]+x0[1])/2
        for i in range(1,n-1):
            x[i] = (1-w)*x0[i] + w*(b[i]+x[i-1]+x0[i+1])/2
        x[n-1] = (1-w)*x0[n-1] + w*(b[n-1]+x[n-2])/2
        dif = np.max(np.abs(x-x0))
        x0 = np.copy(x)
        it+=1

    return it

arr = np.zeros(1000)
x = np.zeros(1000)
for k in range(1,1001):
    arr[k-1] = gauss_seidel_SOR(5, 0.002*k, False)
    x[k-1] = 0.002*k


# testar w local

plt.plot(x[600:801], arr[600:801], marker='o')
plt.title(f"w otimo: {w_opt(5)}")
plt.xlabel('w')
plt.ylabel('numero de iteracoes')
plt.show()
