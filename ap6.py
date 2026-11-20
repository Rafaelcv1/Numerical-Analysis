import numpy as np

r = np.array([2.7, 2, 1.61, 1.2, 1.02])
theta = np.array([48, 67, 83, 108, 126])

n = len(r)
D = np.zeros((n, n))

for i in range(n):
    D[i, 0] = theta[i]

for j in range(1, n):
    for i in range(n - j):
        D[i, j] = (D[i + 1, j - 1] - D[i, j - 1]) / (r[i + j] - r[i])

coeficientes = np.zeros(n)
for k in range(n):
    coeficientes[k] = D[0, k]

termos = []
for k in range(n):
    c = coeficientes[k]
    if k == 0:
        termos.append(f"{c:.6f}")
    else:
        fatores = "".join([f"(r - {r[i]:g})" for i in range(k)])
        sinal = " + " if c >= 0 else " - "
        termos.append(f"{sinal}{abs(c):.6f}{fatores}")

polinomio_str = "P(r) = " + "".join(termos)

print("Coeficientes")
for k in range(n):
    print(f"a_{k} = {coeficientes[k]:.6f}")

print("\nPolinômio Interpolador de Newton")
print(polinomio_str)
