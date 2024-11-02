import java.io.*;
//------este programa calcula los determinantes de una matriz rectangular MxN----------------------------------------------------------------------
//------los tamaños de las submatrices son n:2,3,4,5  estas submatrices son todas las matrices
//------cuadradas posibles de los tamaños mencionados---------------------------------------------
//------el ProcesoParalelo() no imprime las submatrices pues resulta un embrollo, unicamente calcula las submatrices y sus determinantes
//-----------------------------------------------------------------------------------------------------------
//------el ProcesoSerial() imprime las submatrices del siguiente modo:
//------submatriz nxn en (ubicacion del elemento[0][0]) desde donde se "expande" una matriz cuadrada 
//------si  se quiere imprimir dichas submatrices se podria descomentar imprimirMatriz(subM); en linea 87  --------------------------------------
//----------------------------------------------------------------dentro de DeterminantesSubMatricesSerial()---------------------------
//--------------------------------------------------------------------------------------------------------------
public class PC2Preg2DeterminantesSubMatrices{
    private static String FILENAME = "DATAPC2Preg2DeterminantesSubMatrices.TXT";    
    private static int    M = 200;           //filas   M>=N  y  N>=5 pues las submatrices dimension 2,3,4,5   
    private static int    N = 100;           //IMPORTANTE ↑      
    private static int H =4 ;                         
    private static double CADENA;                        
    private static int BLOCK = 5;                       //cantidad de digitos de cada dato BLOCK-1 (3465""3654""9685), pues el byte* separador es "" ,contando de en direccion →
    private static byte[] RECORD = new byte[BLOCK];     //para la lectura de cada dato 
    private static double [][] A=new double[M][N]; 
    private static Thread [] hilos = new Thread[H];
    //----------------------------------------------------------------------------------------------------------
    private static double ObtenerElemento() {
        String CAD;
            CAD = "";
            for(int i=0;i<BLOCK-1;i++) {
                CAD = CAD + (char)(RECORD[i]);
            }
            return Double.parseDouble(CAD);
        }
    //---------------------------------------------------------------------------------------------------------    
    public static double Determinante(double[][] a) {
        int n = a.length;
        double[][] m = new double[n][n];                            //la copia por si acaso
        for (int i = 0; i < n; i++) {
            System.arraycopy(a[i], 0, m[i], 0, n);  //java tiene este metodo eficiente y confiable
        }
        double det =1;
        for(int i=0;i<n-1;i++){
            if(m[i][i]==0){                     //el pivote es 0 throw new ArithmeticException("pivote 0");
                boolean flag = false;           //se cambia la fila para evitar lo anterior
                for(int j=i+1;j<n;j++){
                    if(m[j][i]!=0){
                        double []temp = m[i];
                        m[i] = m[j];
                        m[j] = temp;
                        det *=-1;
                        flag = true;
                    }
                }
                if(!flag){
                    return 0;
                }
            }
			for(int j=i+1;j<n;j++){         //mediante gauss
				double fij=m[j][i]/m[i][i];
				for(int k=0;k<n;k++){
					m[j][k] -=(fij*m[i][k]);
				}
			}
		}                     
        for(int k=0;k<m.length;k++){        //producto de diagonales=determinante
            det*=m[k][k];                                                
        } 
        return det;
    }
    //-----------------------------------------------------------------------------------------------
    // Metodo proorcionado de los codigos hechos en clase
    public static void imprimirMatriz(double[][]M){
        int filas=M.length;
        int columnas=M[0].length;
        System.out.println();
        for(int i=0;i<filas;i++){
                for(int j=0;j<columnas;j++){
                        System.out.printf("%12.2f",M[i][j]);
                }
                System.out.println();
        }
    }
    //----------------------------------------------------------------------------------------------------
    //--este metodo es el que emplea el ProcesoSerial() pues si imprime las submatrices
    public static void DeterminantesSubMatricesSerial(int n) {
        for (int filaInicio = 0; filaInicio <= M - n; filaInicio++) {
            for (int columnaInicio = 0; columnaInicio <= N - n; columnaInicio++) {
                double[][] subM = new double[n][n];
                for (int i = 0; i < n; i++) {
                    System.arraycopy(A[filaInicio + i], columnaInicio, subM[i], 0, n);
                }       //incluso con esta linea comentado ↓ , el tiempo serial es mayor,descomentar para ver detalles(matrices pequeñas)
                //imprimirMatriz(subM);       
                double determinante = Determinante(subM);
                System.out.println("(serial)Submatriz " + n + "x" + n + " en (" + filaInicio + "," + columnaInicio + ") - Determinante: " + determinante);
            }
        }
    }
    //--------------------------------------------------------------------------------------------------------
    //--como se menciona, no se imprimen las submatrices pues resulta confuso
    public static void DeterminantesSubMatricesParalelo(int n) {
        for (int filaInicio = 0; filaInicio <= M - n; filaInicio++) {
            for (int columnaInicio = 0; columnaInicio <= N - n; columnaInicio++) {
                double[][] subM = new double[n][n];
                for (int i = 0; i < n; i++) {
                    System.arraycopy(A[filaInicio + i], columnaInicio, subM[i], 0, n);
                }
                double determinante = Determinante(subM);
                synchronized(System.out){
                    System.out.println("(paralelo)Submatriz " + n + "x" + n + " en (" + filaInicio + "," + columnaInicio + ") - Determinante: " + determinante);}
            }
        }
    }
    //----------------------------------------------------------------------------------------------------
    public static void ProcesoSerial() {
        //hasta submatrices de tamaño 2,3,4,5 pues se tienen H=4 hilos 
        System.out.println("proceso Serial");
        for(int n=2;n<=5;n++){ 
            DeterminantesSubMatricesSerial(n);
        }
        System.out.println("fin de proceso Serial\n");
    }
    //-----------------------un hilo para cada tamaño n:2,3,4,5 -------------------------------------------------------------------------------
    public static void ProcesoParalelo() { 
        hilos[0] = new Thread(()->DeterminantesSubMatricesParalelo(2));
        hilos[1] = new Thread(()->DeterminantesSubMatricesParalelo(3));
        hilos[2] = new Thread(()->DeterminantesSubMatricesParalelo(4));
        hilos[3] = new Thread(()->DeterminantesSubMatricesParalelo(5));
        for(Thread hil:hilos){
            hil.start();
        }
        for (Thread hil:hilos){
            try{
                hil.join();
            }
            catch(InterruptedException e){
                e.printStackTrace();
            }
        }
    }
    //---------------------------------------------------------------------------------------------------------
    private static void AsignarDatosMatriz() {
    long n,P,T;
    int k,i,j;
        try {
             RandomAccessFile RAF = new RandomAccessFile(FILENAME,"r");
             T = RAF.length();
             n = T/BLOCK;
             P = -1;
             k = 0;
             i = 0;
             j = 0;
             while((k<=n-1)&&(P==-1)) {             //El método java.io.RandomAccessFile.seek(long pos) establece el desplazamiento del puntero de archivo
                 RAF.seek(k*BLOCK);                 //, coloca el puntero en esta posicion(posicion=es un entero)                          medido desde el comienzo de este archivo, en el que se produce la siguiente lectura o escritura.
                 RAF.read(RECORD);                  // lee la info de RAF → RECORD , si hay mas parametros llenara la leng de RECORD      Lee hasta b.lengthbytes de datos de este archivo en una matriz de bytes.
                 CADENA = ObtenerElemento();        //RECORD = bytes[BLOCK] osea va a asignar 5 posiciones a CADENA
                 A[i][j] =CADENA;                   //asignar elemento a matriz , un double de 4 cifras
                 if(k==(N*(i+1)-1)){                //            N=4      
                    j = 0;                          // 0 1 2 3 =N*(0+1) -1=3  hace el salto a la siguiente fila
                    i ++;                           //  mientras j 0→N-1    cuando hace el salto j=0 e i ++
                 }
                 else{
                    j++;                            //caso contrario solo j ++
                 }
                 k++;                               // k ++ para posicionar el puntero y acceder al inicio del siguiente
             }                                      //bloque(dato ) en el archivo
             RAF.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }    
    //--------------------------------------------------------------------------------------------------------------
    public static void WriteData() {       //crea la data 
    double X;
    long num;
        try {
            FileWriter FW = new FileWriter(FILENAME);
            for(int i=1;i<=M*N;i++) {
                X = Math.random()*(double)Math.pow(10,BLOCK-2);
                num = (long)Math.pow(10,BLOCK-2) + (long)X; //creando el dato cantidad de digitos=BLOCK-1
                FW.write(num + " ");
            } 
            FW.close();
        }
        catch (IOException E) {
            System.out.print(E.getMessage());
        }
    }
    //------------------------------------------------------------------------------------------------------------------
    public static void main(String[] args) {        //descomentar "//\\" para visualizar matrices pequeñas
        WriteData();
        AsignarDatosMatriz();
        //\\imprimirMatriz(A);
        long start = System.nanoTime();
        ProcesoSerial(); 
        long endSerial = (System.nanoTime()-start)/1000000;
        start = System.nanoTime();
        ProcesoParalelo(); 
        long endParalelo = (System.nanoTime()-start)/1000000;
        System.out.println("\nel tiempo serial es "+endSerial+" milisegundos");
        System.out.println("\nel tiempo paralelo es "+endParalelo+ " milisegundos");
    }
}
