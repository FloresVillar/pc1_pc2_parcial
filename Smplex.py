import numpy as np
from scipy.optimize import linprog
class Simplex:
    def __init__(self, A, b, c):
        self.m = len(b)  # Número de restricciones
        self.n = len(c)  # Número de variables
        self.tableau = np.zeros((self.m + 1, self.n + self.m + 1))

        # Llenar la tabla con las restricciones y la función objetivo
        self.tableau[:-1, :self.n] = A  # Coeficientes de restricciones
        self.tableau[:-1, self.n:self.n + self.m] = np.eye(self.m)  # Variables de holgura
        self.tableau[:-1, -1] = b  # Términos independientes
        self.tableau[-1, :self.n] = -np.array(c)  # Función objetivo (c)

    def pivot(self, row, col):
        # Hacer que el elemento pivote sea 1
        self.tableau[row] /= self.tableau[row, col]
        
        # Hacer ceros en la columna pivote
        for r in range(self.m + 1):
            if r != row:
                self.tableau[r] -= self.tableau[r, col] * self.tableau[row]

    def solve(self):
        while True:
            # Buscar la columna pivote
            col = np.argmin(self.tableau[-1, :-1])
            if self.tableau[-1, col] >= 0:
                break  # Óptimo alcanzado
            
            # Buscar la fila pivote
            ratios = self.tableau[:-1, -1] / self.tableau[:-1, col]
            ratios[ratios <= 0] = np.inf  # Ignorar valores negativos o cero
            row = np.argmin(ratios)

            self.pivot(row, col)

        # Obtener la solución óptima
        solution = np.zeros(self.n)
        for i in range(self.m):
            basic_var = np.where(self.tableau[i, :self.n] == 1)[0]
            if len(basic_var) > 0 and np.sum(self.tableau[i, :self.n]) == 1:
                solution[basic_var[0]] = self.tableau[i, -1]

        return solution, -self.tableau[-1, -1]  # Retornar solución y valor óptimo

if __name__ == "__main__":
    """# Definir las restricciones (A), términos independientes (b) y función objetivo (c)
    A = np.array([[2, 1, 3, 4],
                  [3, 2, 4, 2],
                  [1, 1, 2, 3]])
    b = np.array([200, 300, 150])
    c = np.array([50, 40, 70, 60])

    simplex = Simplex(A, b, c)
    solution, optimal_value = simplex.solve()

    print("Solución óptima:")
    for i in range(len(solution)):
        print(f"x{i + 1} = {solution[i]}")
    print("Valor óptimo de la función objetivo:", optimal_value)"""
    A = np.array([[2, 1, 3, 4],
                  [3, 2, 4, 2],
                  [1, 1, 2, 3]])
    b = np.array([200, 300, 150])
    c = np.array([50, 40, 70, 60])
    res = linprog(-c, A_ub=A, b_ub=b, method='highs')   
    if res.success:
        print("Resultados:")
        print(f"x1 (Mesas): {res.x[0]:.2f}")
        print(f"x2 (Sillas): {res.x[1]:.2f}")
        print(f"x3 (Armarios): {res.x[2]:.2f}")
        print(f"x4 (Estanterías): {res.x[3]:.2f}")
        print(f"Beneficio máximo Z: {-res.fun:.2f}")
    else:
        print("No se pudo encontrar una solución.")  
