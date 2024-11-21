import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.LinkedList;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

public class PC3Preg1InterpolacionDeLagrange {
    private static String FILENAME = "DATAPC3InterpolacionDeLagrange.TXT"; // Archivo de datos
    private static int N = 6;               // Número de filas
    private static int H = 4;               // Número de hilos
    private static double CADENA;
    private static int BLOCK = 5;           // Tamaño de cada dato BLOCK -1
    private static byte[] RECORD = new byte[BLOCK]; // Para la lectura de cada dato
    private static double[][] AA = new double[N][2]; // Matriz de N filas y 2 columnas
    private static LinkedList<Thread> hilos = new LinkedList<>();
    private static double resultadoFinal = 0.0; // Resultado final para la suma de los términos

    private static BlockingQueue<Integer> colaTareas = new ArrayBlockingQueue<>(N);

    public static void Generar() {
        WriteData();
        AsignarDatosMatriz();
    }

    private static double ObtenerElemento() {
        String CAD = "";
        for (int i = 0; i < BLOCK - 1; i++) {
            CAD = CAD + (char) (RECORD[i]);
        }
        return Double.parseDouble(CAD);
    }

    // Crear la data aleatoria y escribirla en el archivo
    public static void WriteData() {
        double X;
        long num;
        try {
            FileWriter FW = new FileWriter(FILENAME);
            for (int i = 1; i <= N * 2; i++) { // Generar datos para N filas y 2 columnas
                X = Math.random() * (double) Math.pow(10, BLOCK - 2);
                num = (long) Math.pow(10, BLOCK - 2) + (long) X;
                FW.write(num + " ");
            }
            FW.close();
        } catch (IOException E) {
            System.out.print(E.getMessage());
        }
    }

    // Asignar datos leídos del archivo a la matriz AA
    private static void AsignarDatosMatriz() {
        long n, P, T;
        int k, i, j;
        try {
            RandomAccessFile RAF = new RandomAccessFile(FILENAME, "r");
            T = RAF.length();
            n = T / BLOCK;
            P = -1;
            k = 0;
            i = 0;
            j = 0;
            while ((k <= n - 1) && (P == -1)) {
                RAF.seek(k * BLOCK); // Posicionar puntero en el archivo
                RAF.read(RECORD);    // Leer el bloque de datos
                CADENA = ObtenerElemento();
                AA[i][j] = CADENA;   // Asignar a la matriz de N filas x 2 columnas
                if (j == 1) {        // Cambiar de fila cada vez que se llenan las 2 columnas
                    j = 0;
                    i++;
                } else {
                    j++;
                }
                k++;
            }
            RAF.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    // Imprimir la matriz AA
    public static void Imprimir(double[][] m) {
        int filas = m.length;
        System.out.println();
        for (int i = 0; i < filas; i++) {
            for (int j = 0; j < m[i].length; j++) {
                System.out.printf("%12.2f", m[i][j]);
            }
            System.out.println();
        }
    }

    //-----------------------------
    // Polinomio de Lagrange Serial
    public static double polinomioLagrangeSerial(double[][] puntos, double x) {
        double resultado = 0.0;
        int n = puntos.length;
        for (int i = 0; i < n; i++) {
            double xi = puntos[i][0];
            double yi = puntos[i][1];
            double terminoLagrange = yi;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    double xj = puntos[j][0];
                    terminoLagrange *= (x - xj) / (xi - xj);
                }
            }
            resultado += terminoLagrange;
        }
        return resultado;
    }

    // Polinomio de Lagrange Paralelo con Threads básicos
    public static double polinomioLagrangeParalelo(double[][] puntos, double x) throws InterruptedException {
        // Se usará la variable resultadoFinal que será modificada por cada hilo
        resultadoFinal = 0.0;

        // Crear los hilos y asignarles tareas
        for (int i = 0; i < N; i++) {
            final int index = i;
            Thread hilo = new Thread(() -> {
                double termino = calcularTerminoLagrange(index, puntos, x);
                // Sincronización para evitar acceso concurrente a la variable resultadoFinal
                synchronized (PC3Preg1InterpolacionDeLagrange.class) {
                    resultadoFinal += termino;
                }
            });
            hilos.add(hilo);
            hilo.start(); // Iniciar el hilo
        }

        // Esperar a que todos los hilos terminen
        for (Thread hilo : hilos) {
            hilo.join(); // Esperar que el hilo termine antes de continuar
        }

        return resultadoFinal;
    }

    // Calcular el término de Lagrange para un índice i
    private static double calcularTerminoLagrange(int i, double[][] puntos, double x) {
        double xi = puntos[i][0];
        double yi = puntos[i][1];
        double terminoLagrange = yi;

        // Calcular el término L_i(x)
        for (int j = 0; j < N; j++) {
            if (j != i) {
                double xj = puntos[j][0];
                terminoLagrange *= (x - xj) / (xi - xj);
            }
        }

        return terminoLagrange;
    }

    // Polinomio de Lagrange con BlockingQueue
    public static double polinomioLagrangeParaleloCola(double[][] puntos, double x) throws InterruptedException {
        resultadoFinal = 0.0;
        
        // Iniciar la cola con los índices de los puntos
        for (int i = 0; i < N; i++) {
            colaTareas.put(i); // El productor coloca los índices en la cola
        }

        // Crear hilos consumidores
        for (int i = 0; i < H; i++) {
            Thread consumidor = new Thread(() -> {
                try {
                    while (!colaTareas.isEmpty()) {
                        int index = colaTareas.take(); // Tomar un índice de la cola
                        double termino = calcularTerminoLagrange(index, puntos, x);
                        // Sincronización para evitar acceso concurrente a la variable resultadoFinal
                        synchronized (PC3Preg1InterpolacionDeLagrange.class) {
                            resultadoFinal += termino;
                        }
                    }
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt(); // Manejar la interrupción
                }
            });
            hilos.add(consumidor);
            consumidor.start(); // Iniciar el hilo consumidor
        }

        // Esperar a que todos los hilos terminen
        for (Thread hilo : hilos) {
            hilo.join(); // Esperar que el hilo termine antes de continuar
        }

        return resultadoFinal;
    }

    // Método principal
    public static void main(String[] args) {
        Generar();
        Imprimir(AA); // Imprimir los puntos

        double xEvaluar = 1234.0; // Valor en el que deseas evaluar el polinomio

        // Evaluación del polinomio de Lagrange de forma serial
        double resultadoSerial = polinomioLagrangeSerial(AA, xEvaluar);
        System.out.println("Polinomio de Lagrange (Serial) evaluado en x = " + xEvaluar + ": " + resultadoSerial);

        // Evaluación del polinomio de Lagrange de forma paralela
        try {
            double resultadoParalelo = polinomioLagrangeParalelo(AA, xEvaluar);
            System.out.println("Polinomio de Lagrange (Paralelo) evaluado en x = " + xEvaluar + ": " + resultadoParalelo);

            // Evaluación del polinomio usando BlockingQueue
            double resultadoCola = polinomioLagrangeParaleloCola(AA, xEvaluar);
            System.out.println("Polinomio de Lagrange (BlockingQueue) evaluado en x = " + xEvaluar + ": " + resultadoCola);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }
}
