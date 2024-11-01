import java.io.*;
import java.util.LinkedList;
//import java.util.ArrayList;
//import java.util.List;
//import java.util.Scanner;
//-----------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------
public class PC2Preg2SubDeterminantesSerialParalelo{
    private static String FILENAME = "DATADETERMINANTES.TXT";
    private static int M = 20;                                          //filas
    private static int N = 5;                                           //columnas necesariamente debe ser mucho menor que Nmax=5
    private static double CADENA;
    private static int BLOCK = 5;                                       //cantidad de digitos de cada dato BLOCK-1 (3465""3654""9685), pues el byte* separador es "" ,contando de en direccion →
    private static byte[] RECORD = new byte[BLOCK];                     //para la lectura de cada dato 
    private static LinkedList<Thread> hilos = new LinkedList<Thread>();
    private static double [][] A=new double[M][N];
    //--------------------------------------------------------------------------------------------------------------------------------
    private static double ObtenerElemento() {
        String CAD;
            CAD = "";
            for(int i=0;i<BLOCK-1;i++) {
                CAD = CAD + (char)(RECORD[i]);
            }
            return Double.parseDouble(CAD);
        }
    //---------------------------------------------------------------------------------------------------------------------------------
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
                 RAF.seek(k*BLOCK);                 //coloca el puntero en esta posicion(posicion=es un entero)                          medido desde el comienzo de este archivo, en el que se produce la siguiente lectura o escritura.
                 RAF.read(RECORD);                  // lee la info de RAF → RECORD , si hay mas parametros llenara laa leng de RECORD      Lee hasta b.lengthbytes de datos de este archivo en una matriz de bytes.
                 CADENA = ObtenerElemento();        //RECORD = bytes[BLOCK] osea va a asignar 5 posiciones a CADENA
                 A[i][j] =CADENA;                   //asignar elemento a matriz , un double de 4 cifras
                 if(k==(N*(i+1)-1)){                //          ejemplo si  N=4      
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
    //-----------------------------------------------------------------------------------------------------------------------------
    public static void Gauss(double [][]a){          //triangulariza una matriz , de modo que para calcular su determinante =
		int n =a.length;                            // producto de los elementos de la diagonal
        for(int i=0;i<n-1;i++){
			for(int j=i+1;j<n;j++){
				double fij=a[j][i]/a[i][i];
				for(int k=0;k<n;k++){
					a[j][k] -=(fij*a[i][k]);
				}
			}
		}
	}
    //-------------------------------------------------------------------------------------------------------------------------------
    public static void ImprimirMatriz(double[][]M){
        int filas=M.length;
        int columnas=M[0].length;
        for(int i=0;i<filas;i++){
                for(int j=0;j<columnas;j++){
                        System.out.printf("%12.2f",M[i][j]);
                }
                System.out.println();
        }
        System.out.println();
    }
    //-----------------------------------------------------------------------------------------------------------------------------------
    public static void WriteData(int N) {                              //crea la data 
    double X;
    long num;
        try {
            FileWriter FW = new FileWriter(FILENAME);
            for(int i=1;i<=M*N;i++) {
                X = Math.random()*(double)Math.pow(10,BLOCK-2);
                num = (long)Math.pow(10,BLOCK-2) + (long)X;         //creando el dato cantidad de digitos=BLOCK-1, 
                FW.write(num + " ");                                  //por ello la potencia de 10 es BLOCK-2
            } 
            FW.close();
        }
        catch (IOException E) {
            System.out.print(E.getMessage());
        }
    }
    //---------------------------------------------------------------------------------------------------------------------------------
    public static void ProcesoSerial(){                 //determinante de las submatrices, condicion necesaria: columnas <= filas
        int filas = A.length;
        int columnas = A[0].length;
        for(int n=2;n<columnas;n++){             //la logica es nxn=dimension de la submatriz, n_maximo=columnas  
            if(n==2){
                for(int f=0;f<=filas-n;f++){     //avance ↓ si M=5 hay 5 ,filas N=4 columnas n=2 las filas 
                    //asignar elementos  Ma←A    //avanzaran hasta 5-2=3 hasta la fila 3 desdde0 asi el 
                    int g=(int)columnas/n+1;    //ultimo bloque sera A  [ M-2=3][]
                    for(int c=0;c<g;c++){        //                     [ M-2+1][]  
                        double [][]Ma = new double[n][n];    // g es la cantidad de submatrices de tamaño n=2 que se tendra
                        for(int i=0;i<n;i++){               //N=4 4/2=2    g= 2+1=3    [] [] [] [] 
                            for(int j=0;j<n;j++){           //                           1  2  3 
                                Ma[i][j] = A[f+i][c*(g-2)+j];// es obvio que i-j asigna a la submtriz que se esta creando*
                            }                               //para las filas si f=0 i:0→1 pues n=2   el manejo de las columnas es mas especial
                        }                                   // g=3 total de grupos entonces g-1=2     si c=0 j=0  c*(2)+j=0 primera columna de A 
                        ImprimirMatriz(Ma);                  //                                                j=1 c*(2)+j=1  segunda columna de A
                        Gauss(Ma);                           // i=0    fila=f=0                         asignando  Ma[0,0] ←A[0,0]y Ma[0,1]←A[0,1]
                        ImprimirMatriz(Ma);                  // i ++ entonces i=1                          la misma logica                                  
                        //determinante                      //                                           Ma[1,0 ]←A[1,0]    Ma[1, 1]←A[1,1]                                       
                        double det = 1;                     // ahora asignamos el segundo bloque c=1   c*(2) +j = 2+j
                        for(int k=0;k<Ma.length;k++){        // la fila(f) no varia todavia  fila=f=0        j=0  A[0,2]   j=1 A[0,3]
                            det*=Ma[k][k];                   //                                          Ma[0,0] ←A[0,2]    Ma[0,1]←A[0,3]                            
                        }                                   //   i++                                      Ma[1,0]←A[1,2]    Ma[1,1]←A[1,3]   y asi para el ultimo bloque g=3
                        System.out.println("determinante de matriz");       //luego el bucle de 'f' avanza f=1
                        System.out.println(det);                                //asigmamos de la fila 1 y fila 2 segun i++
                        System.out.println();                                   // la logica de las columnas se mantiene
                    }                                                           //las filas avanzan hasta la penultima fila M-n=3
                
                }
            }else{
                if(n==3){
                    for(int f=0;f<=filas-n;f++){        //asignar elementos  Ma←A            nxn = 3x3
                        int g=(int)columnas/n+1;        //en este caso los grupos seran 2 N=4/3+1 =2
                        for(int c=0;c<g;c++){
                            double [][]Ma = new double[n][n];
                            for(int i=0;i<n;i++){
                                for(int j=0;j<n;j++){
                                    Ma[i][j] = A[f+i][c*(g-1)+j]; // la logica es la misma g-1=1                         j:0,1,2 < 3
                                }                                  // para  c=0 i=0  c*(1)+j   j=0     Ma[0,0]←A[0,0]  j=1  Ma[0,1]←A[0,1] j=2 Ma[0,2]←A[0,2]   
                            }                                      //             i=1          j=0     Ma[1,0]←A[1,0]  j=1  Ma[1,1]←A[1,1] j=2 Ma[1,2]←A[1,2]                     
                            ImprimirMatriz(Ma);                     //             i=2          j=0     Ma[2,0]←A[2,0]  j=1  Ma[2,1]←A[2,1] j=2 Ma[2,2]←A[2,2]
                            Gauss(Ma);                              // para c=1 i=0 c*(1)+j     j=0     Ma[0,0]←A[0,1]  j=1  Ma[0,1]←A[0,2]  j=2 Ma[0,2]←A[0,3]
                            ImprimirMatriz(Ma);                     //         i=1                            ←A[1,1]             ←A[1,2]            ←A[1,3]  
                            //determinante                         //         i=2                            ←A[2,1]             ←A[2,2]             ←A[2,3]
                            double det = 1;
                            for(int k=0;k<Ma.length;k++){        //avanza fila f++ y asi hasta M=5, n=3   M-n=2(fila) para que puede ir hasta A[f=2 +i=2,]
                                det*=Ma[k][k];                    //                                                                           A[4,]              
                            }
                            System.out.println("determinante de matriz");//para hallar el determinante se triangulariza via gauss
                            System.out.println(det);                //luego el producto de los elementos de la diagonal , 
                            System.out.println();
                        }
                    
                    }
                }
            }
        }
    }

    private static void ProcesoParalelo(){          //paralelizando para n=2 (hilo1)   para n=3 (hilo2) 
        int filas = A.length;
        int columnas = A[0].length;
        for(int n=2;n<columnas;n++){
            if(n==2){
                Thread hil1 = new Thread(new Runnable (){
                    public void run(){
                        final int n=2;
                        for(int f=0;f<=filas-n;f++){
                            //asignar elementos  M←A
                            int g=(int)columnas/n +1;
                            for(int c=0;c<g;c++){
                                double [][]M = new double[n][n];
                                for(int i=0;i<n;i++){
                                    for(int j=0;j<n;j++){
                                        M[i][j] = A[f+i][c*(g-2)+j];
                                    }
                                }
                                ImprimirMatriz(M);
                                Gauss(M);
                                ImprimirMatriz(M);
                                //determinante
                                double det = 1;
                                for(int k=0;k<M.length;k++){
                                    det*=M[k][k];
                                }
                                System.out.println("determinante de matriz");
                                System.out.println(det);
                                System.out.println();
                            }
                        
                        }
                    }
                });
                hilos.add(hil1);
                hil1.start();
            }else{
                if(n==3){
                    Thread hil2 = new Thread(new Runnable(){
                        public void run(){
                        final int n = 3;        
                    for(int f=0;f<=filas-n;f++){
                        //asignar elementos  M←A
                        int g=(int)columnas/n+1;
                        for(int c=0;c<g;c++){
                            double [][]M = new double[n][n];
                            for(int i=0;i<n;i++){
                                for(int j=0;j<n;j++){
                                    M[i][j] = A[f+i][c*(g-1)+j];
                                }
                            }
                            ImprimirMatriz(M);
                            Gauss(M);
                            ImprimirMatriz(M);
                            //determinante
                            double det = 1;
                            for(int k=0;k<M.length;k++){
                                det*=M[k][k];
                            }
                            System.out.println("determinante de matriz");
                            System.out.println(det);
                            System.out.println();
                        }
                    
                    }
                        }
                    });
                    hilos.add(hil2);
                    hil2.start();
                }
            }
            for (Thread hil:hilos) {
				try{
					hil.join();
				}catch(InterruptedException e){
					e.printStackTrace();
				}						
			}
        }
    }
    public static void main(String[] args) {
        WriteData(N);
        AsignarDatosMatriz();
        ImprimirMatriz(A);
        long start= System.nanoTime();
        ProcesoSerial();
        long tiempoSerial =(System.nanoTime()-start)/1000000;
        start=System.nanoTime();
        ProcesoParalelo();
        long tiempoParalelo=(System.nanoTime()-start)/1000000;
        System.out.println("\nlos tiempos serial: "+tiempoSerial+" milisegundos"+"\tparalelo: "+tiempoParalelo+" milisegundos");   
    }
}
