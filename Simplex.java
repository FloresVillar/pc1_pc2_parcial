public class Simplex {

    private double[][] tabla; // Tabla simplex
    private int numConstraints; // Número de restricciones
    private int numVariables;   // Número de variables

    // Constructor: inicializa la tabla simplex con restricciones y variables
    public Simplex(double[][] A, double[] b, double[] c) {
        numConstraints = b.length;
        numVariables = c.length;
        tabla = new double[numConstraints + 1][numVariables + numConstraints + 1];

        // Inicializar tabla: llenar con coeficientes de restricciones (A), término independiente (b), y función objetivo (c)
        for (int i = 0; i < numConstraints; i++) {
            for (int j = 0; j < numVariables; j++) {
                tabla[i][j] = A[i][j];
            }
            tabla[i][numVariables + i] = 1.0; // Variable de holgura
            tabla[i][tabla[0].length - 1] = b[i]; // Término independiente
        }
        for (int j = 0; j < numVariables; j++) {
            tabla[numConstraints][j] = -c[j]; // Función objetivo en la última fila
        }
    }

    // Método para encontrar la columna pivote (variable entrante)
    private int findPivotColumn() {
        int pivotColumn = 0;
        for (int j = 1; j < tabla[0].length - 1; j++) {
            if (tabla[numConstraints][j] < tabla[numConstraints][pivotColumn]) {
                pivotColumn = j;
            }
        }
        return tabla[numConstraints][pivotColumn] < 0 ? pivotColumn : -1;
    }

    // Método para encontrar la fila pivote (variable saliente)
    private int findPivotRow(int pivotColumn) {
        int pivotRow = -1;
        double minRatio = Double.POSITIVE_INFINITY;
        for (int i = 0; i < numConstraints; i++) {
            double ratio = tabla[i][tabla[0].length - 1] / tabla[i][pivotColumn];
            if (tabla[i][pivotColumn] > 0 && ratio < minRatio) {
                minRatio = ratio;
                pivotRow = i;
            }
        }
        return pivotRow;
    }

    // Realiza la operación de pivote en la tabla
    private void performPivot(int pivotRow, int pivotColumn) {
        double pivotValue = tabla[pivotRow][pivotColumn];
        
        // Divide la fila pivote por el valor del pivote para hacer que el pivote sea 1
        for (int j = 0; j < tabla[0].length; j++) {
            tabla[pivotRow][j] /= pivotValue;
        }

        // Hacer ceros en las otras filas en la columna pivote
        for (int i = 0; i <= numConstraints; i++) {
            if (i != pivotRow) {
                double factor = tabla[i][pivotColumn];
                for (int j = 0; j < tabla[0].length; j++) {
                    tabla[i][j] -= factor * tabla[pivotRow][j];
                }
            }
        }

        // Imprimir el estado de la tabla después de cada pivote
        printTable();
    }

    // Método para imprimir el estado de la tabla
    private void printTable() {
        System.out.println("Estado de la tabla:");
        for (double[] fila : tabla) {
            for (double valor : fila) {
                System.out.printf("%10.4f ", valor);
            }
            System.out.println();
        }
        System.out.println();
    }

    // Método para ejecutar el algoritmo Simplex
    public double[] solve() {
        while (true) {
            int pivotColumn = findPivotColumn();
            if (pivotColumn == -1) break; // Óptimo alcanzado

            int pivotRow = findPivotRow(pivotColumn);
            if (pivotRow == -1) throw new ArithmeticException("Solución ilimitada");

            performPivot(pivotRow, pivotColumn);
        }

        // Obtener solución óptima
        double[] solution = new double[numVariables];
        for (int j = 0; j < numVariables; j++) {
            boolean isBasic = false;
            double value = 0;
            for (int i = 0; i < numConstraints; i++) {
                if (tabla[i][j] == 1) {
                    if (isBasic) {
                        isBasic = false;
                        break;
                    }
                    isBasic = true;
                    value = tabla[i][tabla[0].length - 1];
                } else if (tabla[i][j] != 0) {
                    isBasic = false;
                    break;
                }
            }
            if (isBasic) {
                solution[j] = value;
            }
        }
        return solution;
    }

    public static void main(String[] args) {
        double[][] A = {
            {2, 1, 3, 4},
            {3, 2, 4, 2},
            {1, 1, 2, 3}
        };
        double[] b = {200, 300, 150};
        double[] c = {50, 40, 70, 60};

        Simplex simplex = new Simplex(A, b, c);
        double[] solution = simplex.solve();

        System.out.println("Solución óptima:");
        for (int i = 0; i < solution.length; i++) {
            System.out.println("x" + (i + 1) + " = " + solution[i]);
        }
        System.out.println("Valor óptimo de la función objetivo: " + (-simplex.tabla[simplex.numConstraints][simplex.tabla[0].length - 1])); // Ajuste aquí
    }
}
