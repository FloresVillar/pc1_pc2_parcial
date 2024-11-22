import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.LinkedList;
import java.util.Random;

public class PC4Pregunta1DataSetObservaciones {
    private static String FILENAME1 = "DATAPC4Preg1Columna1.TXT";
    private static String FILENAME2 = "DATAPC4Preg1Columna2.TXT";
    private static String FILENAME3 = "DATAPC4Preg1Columna3.TXT"; //se crea la DATA y se guarda en FILENAME
    private static int N = 6; //se recomienda N>2000 pero N<<10000     
    private static int H = 4; //para N=10000 la laptop colapsa    
    private static double CADENA;                       
    private static int BLOCK1 = 3;    
    private static int BLOCK2 = 4;
    private static int BLOCK3 = 5;  
    private static byte[][] records = new byte[3][]; //tamaño de cada dato BLOCK-1     
    private static byte[] RECORD1 = new byte[BLOCK1];     
    private static byte[] RECORD2 = new byte[BLOCK2];
    private static byte[] RECORD3 = new byte[BLOCK3];
    private static double[][] AA = new double[N][N];
    private static LinkedList<Thread> hilos = new LinkedList<Thread>();

    //--------------------------------------------------------------------------------------
    public static void Generar() {
        AsignarDatosMatriz();
        Imprimir(AA);
    }
    //----------------------------------------------------------------------------------
    public static void Inserciones(){
        double[][]COPIASerial = Copiar(AA);
        double[] fila =generarArrayAleatorio(N);
        double[][]AAInsercionSerial= insertarFilaSecuencial(COPIASerial,fila,4);
        double [][]COPIAParalelo = Copiar(AA);
        double[][]AAInsercionParalela = insertarFilaParalela(COPIAParalelo,fila,2);
        Imprimir(AAInsercionSerial);
        Imprimir(AAInsercionParalela); 
    }
    //--------------------------------------------------------------------------------
    public static void Eliminaciones(){
        double[][]COPIASerial = Copiar(AA);
        double[][]AAEliminacionSerial= eliminarFilaSerial(COPIASerial,4);
        double [][]COPIAParalelo = Copiar(AA);
        double[][]AAEliminacionParalela = eliminarFilaParalela(COPIAParalelo,2);
        Imprimir(AAEliminacionSerial);
        Imprimir(AAEliminacionParalela);   
    }
    //-------------------------------------------------------------------------
    //----------------------------------------------------------------------
    public static double[][] eliminarFilaSerial(double[][] matriz, int index) {
        int n = matriz.length;
        int m = matriz[0].length;
        
        // Crear una nueva matriz con una fila menos
        double[][] nuevaMatriz = new double[n - 1][m];
        
        // Copiar las filas antes de la fila a eliminar
        for (int i = 0; i < index; i++) {
            nuevaMatriz[i] = matriz[i];
        }
        
        // Copiar las filas después de la fila a eliminar
        for (int i = index + 1; i < n; i++) {
            nuevaMatriz[i - 1] = matriz[i];
        }
        
        return nuevaMatriz;
    }
    //---------------------------------------------------------------------
    public static double[][] eliminarFilaParalela(double[][] matriz, int index) {
        int n = matriz.length;
        int m = matriz[0].length;
        
        // Crear una nueva matriz con una fila menos
        double[][] nuevaMatriz = new double[n - 1][m];
        
        // Hilo para copiar las filas antes de la fila a eliminar
        Thread hilo1 = new Thread(() -> {
            for (int i = 0; i < index; i++) {
                nuevaMatriz[i] = matriz[i];
            }
        });
        
        // Hilo para copiar las filas después de la fila a eliminar
        Thread hilo2 = new Thread(() -> {
            for (int i = index + 1; i < n; i++) {
                nuevaMatriz[i - 1] = matriz[i];
            }
        });
    
        // Iniciar los hilos
        hilo1.start();
        hilo2.start();
    
        try {
            // Esperar a que los hilos terminen
            hilo1.join();
            hilo2.join();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        
        return nuevaMatriz;
    }
    //---------------------------------------------------------------------------
    public static double[] generarArrayAleatorio(int N) {
        Random rand = new Random();
        double[] array = new double[N];

        // Llenar el array respetando el patrón
        for (int i = 0; i < N; i += 3) {
            // Generar valores aleatorios para los 3 bloques
            double valor1 = Math.round((rand.nextDouble() * 1000) * 100.0) / 100.0;  // Primer valor
            double valor2 = Math.round((rand.nextDouble() * 500) * 100.0) / 100.0;   // Segundo valor
            double valor3 = Math.round((rand.nextDouble() * 2000) * 100.0) / 100.0;  // Tercer valor

            // Asignar el mismo valor a las tres posiciones
            array[i] = valor1;
            array[i + 1] = valor2;
            array[i + 2] = valor3;
        }

        return array; // Devuelve el array generado
    }
    //----------------------------------------------------------------------
    //-----------------------------------------------------------------------
    public static void Ordenaciones(){
        double[][]COPIASerial = Copiar(AA);
        double[][]AAOrdenadoSerial= ordenarPorPrimeraColumna(COPIASerial);
        double [][]COPIAParalelo = Copiar(AA);
        double[][]AAOrdenadoParalela = ordenarPorPrimeraColumnaParalelo(COPIAParalelo);
        Imprimir(AAOrdenadoSerial);
        Imprimir(AAOrdenadoParalela);
    }
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    // Método para insertar filas secuencialmente
    public static double[][] insertarFilaSecuencial(double[][] matriz, double[] fila, int index) {
        int n = matriz.length;
        int m = matriz[0].length;
        double[][] nuevaMatriz = new double[n + 1][m];

        // Copiar filas antes de la fila a insertar
        for (int i = 0; i < index; i++) {
            nuevaMatriz[i] = matriz[i];
        }

        // Insertar la nueva fila en la posición indicada
        nuevaMatriz[index] = fila;

        // Copiar las filas después de la fila insertada
        for (int i = index; i < n; i++) {
            nuevaMatriz[i + 1] = matriz[i];
        }

        return nuevaMatriz;
    }

    // Método para insertar filas en paralelo
    public static double[][] insertarFilaParalela(double[][] matriz, double[] fila, int index) {
        int n = matriz.length;
        int m = matriz[0].length;
        double[][] nuevaMatriz = new double[n + 1][m];

        Thread hilo1 = new Thread(() -> {
            // Copiar filas antes de la fila a insertar
            for (int i = 0; i < index; i++) {
                nuevaMatriz[i] = matriz[i];
            }
        });

        Thread hilo2 = new Thread(() -> {
            // Insertar la nueva fila en la posición indicada
            nuevaMatriz[index] = fila;
        });

        Thread hilo3 = new Thread(() -> {
            // Copiar las filas después de la fila insertada
            for (int i = index; i < n; i++) {
                nuevaMatriz[i + 1] = matriz[i];
            }
        });

        // Iniciar los hilos
        hilo1.start();
        hilo2.start();
        hilo3.start();

        try {
            // Esperar a que todos los hilos terminen
            hilo1.join();
            hilo2.join();
            hilo3.join();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        return nuevaMatriz;
    }
    //-------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------
    public static double[][] Copiar(double[][]A){
		int n=A.length;
		double [][]COPIA=new double[n][n];
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				COPIA[i][j]=A[i][j];
			}
		}
		return COPIA;
	}
    //-----------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------
    // Método secuencial para ordenar por la primera columna
    public static double[][] ordenarPorPrimeraColumna(double[][] matriz) {
        int n = matriz.length;

        // Algoritmo de burbuja secuencial
        for (int i = 0; i < n - 1; i++) {
            for (int j = 0; j < n - 1 - i; j++) {
                // Comparar los elementos de la primera columna
                if (matriz[j][0] > matriz[j + 1][0]) {
                    // Intercambiar las filas
                    double[] temp = matriz[j];
                    matriz[j] = matriz[j + 1];
                    matriz[j + 1] = temp;
                }
            }
        }
        return matriz; // Devuelve la matriz ordenada
    }

    // Método paralelo para ordenar por la primera columna
    public static double[][] ordenarPorPrimeraColumnaParalelo(double[][] matriz) {
        int n = matriz.length;
        Thread[] hilos = new Thread[n - 1]; // Número de hilos (n-1 pasadas del algoritmo de burbuja)

        // Crear hilos
        for (int i = 0; i < n - 1; i++) {
            final int idx = i;  // Índice del hilo
            hilos[i] = new Thread(() -> {
                for (int j = 0; j < n - 1 - idx; j++) {
                    // Comparar los elementos de la primera columna
                    if (matriz[j][0] > matriz[j + 1][0]) {
                        // Intercambiar las filas
                        double[] temp = matriz[j];
                        matriz[j] = matriz[j + 1];
                        matriz[j + 1] = temp;
                    }
                }
            });
        }

        // Iniciar todos los hilos
        for (Thread hilo : hilos) {
            hilo.start();
        }

        // Esperar a que todos los hilos terminen
        for (Thread hilo : hilos) {
            try {
                hilo.join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
        return matriz; // Devuelve la matriz ordenada
    }
    //--------------------------------------------------------------------------------------
    @FunctionalInterface
    interface DoubleSupplier {
        double get();
    }

    //--------------------------------------------------------------------------------------
    private static double ObtenerElemento1() {
        String CAD = "";
        for (int i = 0; i < BLOCK1; i++) {  //modificado: leyendo BLOCK1 bytes
            CAD = CAD + (char)(RECORD1[i]);
        }
        return Double.parseDouble(CAD);
    }

    private static double ObtenerElemento2() {
        String CAD = "";
        for (int i = 0; i < BLOCK2; i++) {  //modificado: leyendo BLOCK2 bytes
            CAD = CAD + (char)(RECORD2[i]);
        }
        return Double.parseDouble(CAD);
    }

    private static double ObtenerElemento3() {
        String CAD = "";
        for (int i = 0; i < BLOCK3; i++) {  //modificado: leyendo BLOCK3 bytes
            CAD = CAD + (char)(RECORD3[i]);
        }
        return Double.parseDouble(CAD);
    }

    //--------------------------------------------------------------------------------------
    private static DoubleSupplier[] obtenerElementos = new DoubleSupplier[] {
        () -> ObtenerElemento1(),
        () -> ObtenerElemento2(),
        () -> ObtenerElemento3()
    };

    //--------------------------------------------------------------------------------------
    public static void WriteData1() { //crea la data 
        double X;
        long num;
        try {
            FileWriter FW = new FileWriter(FILENAME1);
            for (int i = 1; i <= N; i++) {
                X = Math.random() * (double) Math.pow(10, BLOCK1 - 2);
                num = (long) Math.pow(10, BLOCK1 - 2) + (long) X; //creando el dato cantidad de digitos=BLOCK-1
                FW.write(num + " ");
            }
            FW.close();
        } catch (IOException E) {
            System.out.print(E.getMessage());
        }
    }

    //--------------------------------------------------------------------------------------
    public static void WriteData2() { //crea la data 
        double X;
        long num;
        try {
            FileWriter FW = new FileWriter(FILENAME2);
            for (int i = 1; i <= N; i++) {
                X = Math.random() * (double) Math.pow(10, BLOCK2 - 2);
                num = (long) Math.pow(10, BLOCK2 - 2) + (long) X; //creando el dato cantidad de digitos=BLOCK-1
                FW.write(num + " ");
            }
            FW.close();
        } catch (IOException E) {
            System.out.print(E.getMessage());
        }
    }

    //--------------------------------------------------------------------------------------
    public static void WriteData3() { //crea la data 
        double X;
        long num;
        try {
            FileWriter FW = new FileWriter(FILENAME3);
            for (int i = 1; i <= N; i++) {
                X = Math.random() * (double) Math.pow(10, BLOCK3 - 2);
                num = (long) Math.pow(10, BLOCK3 - 2) + (long) X; //creando el dato cantidad de digitos=BLOCK-1
                FW.write(num + " ");
            }
            FW.close();
        } catch (IOException E) {
            System.out.print(E.getMessage());
        }
    }

    //--------------------------------------------------------------------------------------
    private static void AsignarDatosMatriz() {
        WriteData1();
        WriteData2();
        WriteData3();
        long n, P, T;
        int k, i;
        String[] files = {FILENAME1, FILENAME2, FILENAME3};
        int[] block = {BLOCK1, BLOCK2, BLOCK3};
        records[0] = RECORD1;
        records[1] = RECORD2;
        records[2] = RECORD3;
        for (int global = 0; global < N;) {
            for (int cont = 0; cont < 3; cont++) {
                try {
                    RandomAccessFile RAF = new RandomAccessFile(files[cont], "r");
                    T = RAF.length();
                    n = T / block[cont];
                    P = -1;
                    k = 0;
                    i = 0;
                    while ((k <= n - 1) && (P == -1)) {
                        RAF.seek(k * block[cont]);
                        RAF.read(records[cont]);
                        CADENA = obtenerElementos[cont].get();
                        AA[i][cont + global] = CADENA;
                        i++;
                        if (i == N) break;
                        k++;
                    }
                    RAF.close();
                } catch (IOException ex) {
                    ex.printStackTrace();
                }
            }
            global += 3;
        }
    }

    //--------------------------------------------------------------------------------------
    public static synchronized void Imprimir(double[][] m) {
        System.out.println();
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                System.out.printf("%12.2f", m[i][j]);
            }
            System.out.println();
        }
        System.out.println();
    }
    //-------------------------------------------------------------------
    // Método para modificar una fila de forma serial
public static double[][] modificarFilaSerial(double[][] matriz, double[] filaModificada, int index) {
    int n = matriz.length;
    int m = matriz[0].length;
    
    // Crear una nueva matriz
    double[][] nuevaMatriz = new double[n][m];
    
    // Copiar las filas anteriores a la fila modificada
    for (int i = 0; i < index; i++) {
        nuevaMatriz[i] = matriz[i];
    }
    
    // Reemplazar la fila en la posición index
    nuevaMatriz[index] = filaModificada;
    
    // Copiar las filas posteriores a la fila modificada
    for (int i = index + 1; i < n; i++) {
        nuevaMatriz[i] = matriz[i];
    }
    
    return nuevaMatriz;
}
//-----------------------------------------------------------------------
public static double[][] modificarFilaParalela(double[][] matriz, double[] filaModificada, int index) {
    int n = matriz.length;
    int m = matriz[0].length;
    
    // Crear una nueva matriz
    double[][] nuevaMatriz = new double[n][m];
    
    // Hilo para copiar las filas anteriores a la fila modificada
    Thread hilo1 = new Thread(() -> {
        for (int i = 0; i < index; i++) {
            nuevaMatriz[i] = matriz[i];
        }
    });
    
    // Hilo para reemplazar la fila en la posición indicada
    Thread hilo2 = new Thread(() -> {
        nuevaMatriz[index] = filaModificada;
    });
    
    // Hilo para copiar las filas posteriores a la fila modificada
    Thread hilo3 = new Thread(() -> {
        for (int i = index + 1; i < n; i++) {
            nuevaMatriz[i] = matriz[i];
        }
    });
    
    // Iniciar los hilos
    hilo1.start();
    hilo2.start();
    hilo3.start();

    try {
        // Esperar a que todos los hilos terminen
        hilo1.join();
        hilo2.join();
        hilo3.join();
    } catch (InterruptedException e) {
        e.printStackTrace();
    }
    
    return nuevaMatriz;
}
//--------------------------------------------------------------
public static void Modificaciones(){
    double[][]COPIASerial = Copiar(AA);
    double[] fila =generarArrayAleatorio(N);
    double[][]AAModificacionSerial= modificarFilaSerial(COPIASerial,fila,4);
    double [][]COPIAParalelo = Copiar(AA);
    double[][]AAModificacionParalela = modificarFilaParalela(COPIAParalelo,fila,2);
    Imprimir(AAModificacionSerial);
    Imprimir(AAModificacionParalela); 
}
    //--------------------------------------------------------------------------------------
    public static void main(String[] args) {
        Generar();
        Ordenaciones();
        Inserciones();
        Eliminaciones();
        Modificaciones();
    }
}
