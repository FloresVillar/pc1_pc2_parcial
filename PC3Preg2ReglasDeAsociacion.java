import java.io.*;
import java.util.*;
import java.util.concurrent.*;

public class PC3Preg2ReglasDeAsociacion {
    private static String FILENAME = "DATAPC3Preg2ReglasDeAsociacion.TXT";    
    private static int M = 20; // filas
    private static int N = 10; // columnas
    private static int H = 4;  // hilos
    private static double CADENA; // variable auxiliar
    private static int B = 5; // cantidad de dígitos =B-1 de cada dato 
    private static byte[] RECORD = new byte[B]; // para la lectura de cada dato 
    private static double[][] A = new double[M][N]; // matriz de datos 
    private static BlockingQueue<double[]> COLA;

    // Variables globales para las frecuencias acumuladas
    private static int frecuenciaTotalRegla1 = 0;
    private static int frecuenciaTotalRegla2 = 0;
    private static final Object lock = new Object(); // Para sincronizar acceso a las frecuencias

    // Obtener elemento de la lectura del archivo
    private static double ObtenerElemento() {
        String CAD = "";
        for (int i = 0; i < B - 1; i++) {
            CAD = CAD + (char) (RECORD[i]);
        }
        return Double.parseDouble(CAD);
    }

    // Asignar datos a la matriz A
    private static void AsignarDatosMatriz() {
        long n, P, T;
        int k = 0, i = 0, j = 0;
        try {
            RandomAccessFile RAF = new RandomAccessFile(FILENAME, "r");
            T = RAF.length();
            n = T / B;
            P = -1;
            while ((k <= n - 1) && (P == -1)) {
                RAF.seek(k * B); // Pone el puntero en la posición k
                RAF.read(RECORD);
                CADENA = ObtenerElemento();
                A[i][j] = CADENA;
                if (k == (N * (i + 1) - 1)) {
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

    // Crear la data de manera aleatoria
    public static void WriteData() {
        double X;
        long num;
        try {
            FileWriter FW = new FileWriter(FILENAME);
            for (int i = 1; i <= M * N; i++) {
                X = Math.random() * (double) Math.pow(10, B - 2);
                num = (long) Math.pow(10, B - 2) + (long) X;
                FW.write(num + " ");
            }
            FW.close();
        } catch (IOException E) {
            System.out.print(E.getMessage());
        }
    }

    // Imprimir la matriz
    public static void ImprimirMatriz(double[][] M) {
        int filas = M.length;
        int columnas = M[0].length;
        System.out.println();
        for (int i = 0; i < filas; i++) {
            for (int j = 0; j < columnas; j++) {
                System.out.printf("%12.2f", M[i][j]);
            }
            System.out.println();
        }
    }

    // Almacenar las filas en la cola
    public static void AlmacenarEnCola() {
        COLA = new LinkedBlockingQueue<>();
        for (int i = 0; i < M; i++) {
            try {
                COLA.put(A[i]);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
            }
        }
    }

    // Imprimir la cola
    public static void ImprimirDeCola() {
        COLA.forEach(row -> {
            for (double value : row) {
                System.out.print(value + " ");
            }
            System.out.println();
        });
    }

    // Función para encontrar reglas de asociación en un conjunto de transacciones (filas)
    public static void encontrarReglasDeAsociacion(List<double[]> transacciones) {
        int frecuenciaRegla1 = 0;
        int frecuenciaRegla2 = 0;

        for (double[] fila : transacciones) {
            int valorColumna1 = (int) fila[0];
            int valorColumna2 = (int) fila[1];

            // Regla 1: Ambos elementos son múltiplos de 5
            if (valorColumna1 % 5 == 0 && valorColumna2 % 5 == 0) {
                frecuenciaRegla1++;
            }

            // Regla 2: Columna 1 termina en 1 y Columna 2 termina en 7
            if (valorColumna1 % 10 == 1 && valorColumna2 % 10 == 7) {
                frecuenciaRegla2++;
            }
        }

        // Sincronizar y acumular frecuencias globales
        synchronized (lock) {
            frecuenciaTotalRegla1 += frecuenciaRegla1;
            frecuenciaTotalRegla2 += frecuenciaRegla2;
        }
    }

    // Crear hilos para procesar la cola
    public static void procesarEnHilos() throws InterruptedException {
        List<Thread> threads = new ArrayList<>();

        // Dividir las transacciones entre los hilos
        for (int i = 0; i < H; i++) {
            final int hiloIndex = i;
            Thread t = new Thread(() -> {
                List<double[]> transacciones = new ArrayList<>();
                int numTransaccionesPorHilo = M / H;

                for (int i1 = 0; i1 < numTransaccionesPorHilo; i1++) {
                    try {
                        transacciones.add(COLA.take());
                    } catch (InterruptedException e) {
                        Thread.currentThread().interrupt();
                    }
                }

                // Encontrar reglas de asociación para las transacciones de este hilo
                encontrarReglasDeAsociacion(transacciones);
            });
            threads.add(t);
            t.start();
        }

        // Esperar que todos los hilos terminen
        for (Thread t : threads) {
            t.join();
        }

        // Imprimir las frecuencias acumuladas después de que todos los hilos hayan terminado
        System.out.println("Regla 1: {Columna 1 y Columna 2 son múltiplos de 5} con frecuencia total: " + frecuenciaTotalRegla1);
        System.out.println("Regla 2: {Columna 1 termina en 1 y Columna 2 termina en 7} con frecuencia total: " + frecuenciaTotalRegla2);
    }

    // Método principal
    public static void main(String[] args) throws InterruptedException {
        WriteData();
        AsignarDatosMatriz();
        ImprimirMatriz(A);
        AlmacenarEnCola();
        ImprimirDeCola();
        // Procesar la cola en 4 hilos
        procesarEnHilos();
    }
}
