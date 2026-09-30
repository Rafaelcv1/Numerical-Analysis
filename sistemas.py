import numpy as np

# substituição inversa (matriz diagonal superior)
def inverseSubstitution(matrix_a: np.ndarray, matrix_b: np.ndarray) -> np.ndarray:
    system_lenght : int = len(matrix_a)
    answer = np.zeros(system_lenght)

    answer[-1] = matrix_b[-1]/matrix_a[-1][-1]
    s = 0

    for k in range(system_lenght-1, -1, -1):
        s = matrix_b[k]
        for j in range(k+1, system_lenght, 1):
            s = s - matrix_a[k][j]*answer[j]
        answer[k] = s/matrix_a[k][k]
    return answer

# substituição direta (matriz diagonal inferior)
def directsubstitution(matrix_a: list[list[float]], matrix_b: list[float]) -> np.ndarray:
    system_lenght : int = len(matrix_a)
    answer = np.zeros(system_lenght)
    answer[0] = matrix_b[0]/matrix_a[0][0]
    s = 0

    for k in range(1, system_lenght, 1):
        s = matrix_b[k]
        for j in range(0, k, 1):
            s = s - matrix_a[k][j]*answer[j]
        answer[k] = s/matrix_a[k][k]
    return answer

# eliminação de Gauss
class GaussAnswer:
    def __init__(self, matrix_a: np.ndarray, matrix_b: np.ndarray, system_length: int):
        self.matrix_a = matrix_a
        self.matrix_b = matrix_b
        self.system_lenght = system_length

    system_lenght: int = 0
    matrix_a: np.ndarray= np.zeros((system_lenght, system_lenght))
    matrix_b: np.ndarray= np.zeros((system_lenght, system_lenght))

def gaussElimination(matrix_a: np.ndarray, matrix_b: np.ndarray) -> GaussAnswer:
    system_lenght : int = len(matrix_a)
    matrix_m = np.zeros((system_lenght, system_lenght))

    for k in range(0, system_lenght-1, 1):
        for i in range(k+1, system_lenght, 1):
            matrix_m[i][k] = matrix_a[i][k]/matrix_a[k][k]
            for j in range(k+1, system_lenght, 1):
                matrix_a[i][j] = matrix_a[i][j] - matrix_m[i][k]*matrix_a[k][j]
            matrix_b[i] = matrix_b[i] - matrix_m[i][k]*matrix_b[k]

    return GaussAnswer(matrix_a, matrix_b, system_lenght)

# eliminação tridiagonal de Gauss
def triadiagonalGaussElimination(d: list[float], a: list[float], c: list[float], b: list[float]) -> np.ndarray:
    system_lenght: int = len(d)
    answer: np.ndarray = np.zeros(system_lenght)

    for k in range(0, system_lenght-1, 1):
        d[k+1]=d[k+1]-a[k+1]/d[k]*c[k]
        b[k+1]=b[k+1]-a[k+1]/d[k]*b[k]

    answer[system_lenght-1] = b[system_lenght-1]/d[system_lenght-1]

    for k in range(system_lenght-2, -1, -1):
        answer[k] = (b[k]-c[k]*answer[k+1])/d[k]
    return answer

# geração da matriz tridiagonal
def genTridiagonalBarMatrix(t0: float, tn: float, n: int) -> np.ndarray:
    answer = np.zeros((4, n-1))
    answer[0] = np.full(n-1, 2)
    answer[1] = np.full(n-1, -1)
    answer[2] = np.full(n-1, -1)
    answer[3][0]=t0
    answer[3][n-2]=tn

    return answer

### TESTING AREA

# tesntando substituição inversa
# triang_sup_a = [
#     [1.0, 3.0, -3.0],
#     [0.0, 3.0, 4.0],
#     [0.0, 0.0, 9.0]
# ]
# triang_sup_b = [3.0, 5.0, 7.0]
# ans = inverseSubstitution(np.array(triang_sup_a), np.array(triang_sup_b)).tolist()
# for i in range(3):
#     print(ans[i])

# testando substituição direta
# triang_inf_a = [
#     [1.0, 0.0, 0.0],
#     [2.0, 3.0, 0.0],
#     [3.0, 9.0, 9.0]
# ]
# triang_inf_b = [3.0, 5.0, 7.0]
# ans = directsubstitution(triang_inf_a, triang_inf_b).tolist()
# for i in range(3):
#     print(ans[i])

# testando eliminação gaussiana
# gaussianElim_a = [
#     [1.0, 2.0, 11.0],
#     [2.0, 3.0, -2.0],
#     [3.0, 9.0, 9.0]
# ]
# gaussianElim_b = [3.0, 5.0, 7.0]
# scaled = gaussElimination(np.array(gaussianElim_a), np.array(gaussianElim_b))
# sol = inverseSubstitution(scaled.matrix_a, scaled.matrix_b).tolist()
#
# for i in range(3):
#     print(sol[i])

# testando eliminação gaussiana em esparsas

t0 = 0.0
tn = 50.0
n = 5

problem = genTridiagonalBarMatrix(t0, tn, n).tolist()
answer = triadiagonalGaussElimination(problem[0], problem[1], problem[2], problem[3]).tolist()
print(f"T0: {t0}")
for i in range(len(answer)):
    print(f"T{i+1}: {answer[i]}")
print(f"TN: {tn}")
