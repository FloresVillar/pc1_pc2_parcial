import java.io.*;
//import java.util.ArrayList;
import java.util.LinkedList;
//import java.util.List;
//import java.util.Scanner;
public class SubDeterminantesParalalelo{
    private static String FILENAME = "PARALELO_DETERMINANTES.TXT";
    private static int M = 10;  //filas
    private static int N = 4;//columnas
    private static double CADENA;
    private static int BLOCK = 5; //cantidad de digitos de cada dato BLOCK-1 (3465""3654""9685), pues el byte* separador es "" ,contando de en direccion →
    private static byte[] RECORD = new byte[BLOCK]; //para la lectura de cada dato 
    private static LinkedList<Thread> hilos = new LinkedList<Thread>();
    private static double [][] A=new double[M][N];
    
    private static double ObtenerElemento() {
        String CAD;
            CAD = "";
            for(int i=0;i<BLOCK-1;i++) {
                CAD = CAD + (char)(RECORD[i]);
            }
            return Double.parseDouble(CAD);
        }

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
             //System.out.println(n);
             while((k<=n-1)&&(P==-1)) {         //El método java.io.RandomAccessFile.seek(long pos) establece el desplazamiento del puntero de archivo
                 RAF.seek(k*BLOCK);             //, coloca el puntero en esta posicion(posicion=es un entero)                          medido desde el comienzo de este archivo, en el que se produce la siguiente lectura o escritura.
                 RAF.read(RECORD);              // lee la info de RAF → RECORD , si hay mas parametros llenara laa leng de RECORD      Lee hasta b.lengthbytes de datos de este archivo en una matriz de bytes.
                 CADENA = ObtenerElemento();    //RECORD = bytes[BLOCK] osea va a asignar 5 posiciones a CADENA
                 A[i][j] =CADENA;               //asignar elemento a matriz , un double de 4 cifras
                 //System.out.println(CADENA+"k="+k +" i="+i+" j="+j);
                 if(k==(N*(i+1)-1)){
                    j = 0;
                    i ++;
                 }
                 else{
                    j++;
                 }
                 k++; 
             }
             RAF.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    public static void Gauss(double [][]a){
		int n =a.length;
        for(int i=0;i<n-1;i++){
			for(int j=i+1;j<n;j++){
				double fij=a[j][i]/a[i][i];
				for(int k=0;k<n;k++){
					a[j][k] -=(fij*a[i][k]);
				}
			}
		}
	}
    

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

    public static void WriteData(int N) {       //crea la data 
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

    private static void Pregunta2(){//determinante de las submatrices
        WriteData(N);
        AsignarDatosMatriz();
        ImprimirMatriz(A);
        //int menor = M >= N ? N : M;
        if(M < N){
            //transponer matriz
        }
        //de modo que columnas siempre menor o igual a filas
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
      Pregunta2();   
    }
}
