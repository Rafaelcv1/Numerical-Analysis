import matplotlib.pyplot as plt
import numpy as np

# a)

x = np.array([0.5, 1, 2, 3, 4, 5, 6])
y = np.array([32, 47, 65, 92, 132, 190, 275])

plt.scatter(x, y)
plt.plot(x, y, label="dados discretos")

plt.title("Número de Bactérias por Tempo")
plt.xlabel("Tempo (horas)")
plt.ylabel("Bactérias")

# O gráfico de dispersão apresenta um crescimento crescente e não linear do número de bactérias ao longo do tempo
# A inclinação aumenta à medida que o tempo passa, indicando que o crescimento é cada vez mais rápido
# Comportamento compatível com o de uma função exponencial, embora o gráfico não comprove que esse seja o modelo exato

# b)

# linearizar y = ab^x
# Y = a0+a1X
# ln(y) = ln(ab^x); ln(y) = ln(a)+ln(b^x); ln(y) = ln(a)+ln(b)*x
# Y = ln(y); a0 = ln(a); a1 = ln(b); X = x
# a = e^a0; b = e^a1
# g(x) = 29.6*1,4)^x

# linearizar y = ax^b
# Y = a0+a1X
# ln(y) = ln(ax^b); ln(y) = ln(a)+ln(x^b); ln(y) = ln(a)+b*ln(x)
# Y = ln(y); a0 = ln(a); a1 = b; X = ln(x)
# a = e^a0; b = a1; x = e^X
# g(x) = 47*x^0.8

# c)

# métodos dos mínimos quadrados:
N = len(x)

# y = ab^x
X_exp = x
Y_exp = np.log(y)

M_exp = np.array([
    [N, np.sum(X_exp)],
    [np.sum(X_exp), np.sum(X_exp**2)]
])

v_exp = np.array([
    np.sum(Y_exp),
    np.sum(X_exp*Y_exp)
])

a0_exp, a1_exp = np.linalg.solve(M_exp, v_exp)

a_exp = np.exp(a0_exp)
b_exp = np.exp(a1_exp)

y_exp = a_exp * b_exp**x

# residuos e soma do quadrados dos resíduos
res_exp = y - y_exp
r_exp = np.sum(res_exp**2)

# y = ax^b
X_pot = np.log(x)
Y_pot = np.log(y)

M_pot = np.array([
    [N, np.sum(X_pot)],
    [np.sum(X_pot), np.sum(X_pot**2)]
])

v_pot = np.array([
    np.sum(Y_pot),
    np.sum(X_pot*Y_pot)
])

a0_pot, a1_pot = np.linalg.solve(M_pot, v_pot)

a_pot = np.exp(a0_pot)
b_pot = a1_pot

y_pot = a_pot * x**b_pot

# residuos e soma do quadrados dos resíduos
res_pot = y - y_pot
r_pot= np.sum(res_pot**2)

# plotar os graficos:
x_curva = np.linspace(min(x), max(x), 300)

# y=ab^x
y_exp_curva = a_exp * b_exp**x_curva
plt.plot(x_curva, y_exp_curva, label="Exponencial: ab^x")

# plotar y=ax^b
y_pot_curva = a_pot * x_curva**b_pot
plt.plot(x_curva, y_pot_curva, label="Potência: ax^b")

# plotar os erros calculados
plt.scatter(x, res_exp, marker="o", label="Resíduos ab^x")
plt.scatter(x, res_pot, marker="o", label="Resíduos ax^b")
plt.axhline(0, linestyle="--", color="black")

plt.grid(True)
plt.legend()
plt.show()

print("Soma dos quadrados dos resíduos:")
print(f"Exponencial: {r_exp:.4f}")
print(f"Potência: {r_pot:.4f}")

# O modelo exponencial apresenta resíduos individuais relativamente pequenos,
# enquanto o modelo de potência apresenta diferenças maiores, principalmente nos últimos tempos da tabela.
# Isso reforça a escolha do modelo exponencial para representar o crescimento bacteriano.
