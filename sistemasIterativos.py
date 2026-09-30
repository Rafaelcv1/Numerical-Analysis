import numpy as np

n = 4
x0 = np.zeros(n)
x = np.zeros(n)
b = np.array([6,25,-11,15])

it = 1
itmax = 100
tol = 1.0e-6
dif = 1

print()
print("Jacobi")
while (dif > tol and it < itmax):
    x[0] = (b[0]+x0[1]-2*x0[2])/10
    x[1] = (b[1]+x0[0]+x0[2]-3*x0[3])/11
    x[2] = (b[2]-2*x0[0]+x0[1]+x0[3])/10
    x[3] = (b[3]-3*x0[1]+x0[2])/8
    dif = max(abs(x-x0))
    print(it, dif, x)
    x0 = np.copy(x)
    it+=1

print()
print("Gauss-seidell")

x0 = np.zeros(n)
x = np.zeros(n)

it = 1
dif = 1

while (dif > tol and it < itmax):
  x[0] = (b[0]+x0[1]-2*x0[2])/10
  x[1] = (b[1]+x[0]+x0[2]-3*x0[3])/11
  x[2] = (b[2]-2*x[0]+x[1]+x0[3])/10
  x[3] = (b[3]-3*x[1]+x[2])/8
  dif = max(abs(x-x0))
  print(it, dif, x)
  x0 = np.copy(x)
  it+=1

print()
print("sistema iterativo")
n = 5
x0 = np.zeros(n)
x = np.zeros(n)
b = np.zeros(n)
b[0] = 0
b[n-1] = 1

it = 1
itmax = 100
tol = 1.0e-6
dif = 1

print()
print("Jacobi")
while (dif > tol and it < itmax):
    x[0] = (b[0]+x0[1])/2
    for i in range(1,n-1):
        x[i] = (b[i]+x0[i-1]+x0[i+1])/2
    x[n-1] = (b[n-1]+x0[n-2])/2

    dif = max(abs(x-x0))
    print(it, dif, x)
    x0 = np.copy(x)
    it+=1

x0 = np.zeros(n)
x = np.zeros(n)

it = 1
dif = 1
w = 1

print()
print("Gauss w = 1")
while (dif > tol and it < itmax):
  x[0] = (1-w)*x0[0] + w*(b[0]+x0[1])/2
  for i in range(1,n-1):
    x[i] = (1-w)*x0[i] + w*(b[i]+x0[i-1]+x[i+1])/2
  x[n-1] = (1-w)*x0[0] + w*(b[n-1]+x[n-2])/2
  dif = max(abs(x-x0))
  print(it, dif, x)
  x0 = np.copy(x)
  it+=1
