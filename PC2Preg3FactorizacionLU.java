//import java.io.*;
//import javax.swing.*;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.*;
//------																				  A11 A12								
//-----calculo de la fatorizacion LU serial en doolitle()------LU paralelo en LuBloques() A21 A22 
//-----N TIENE QUE PAR pues la matriz de va a dividir en 4 partes y reemzamblar luego de trabajar los bloques--------
//-----se puede añadir algo de complejidad para manejar matrices impares(no es el caso actual)-----------------------
//---el resultado del proceso serial se obtiene en main serial()--usa directamente doolitle con AA como arg
//---el resultado del proceso paralelo se obtiene en main paralelo()--usa doolitle por bloques
public class PC2Preg3FactorizacionLU{
	private static String FILENAME = "DATAPC2Preg3FactorizacionLU.TXT"; 	//se crea la DATA y se guarda en FILENAME
    private static int    N =2500;              //se recomienda N>2000 pero N<<10000     
    private static int H =4 ;                   //para N=10000 la laptop colapsa    
    private static double CADENA;                       
    private static int BLOCK = 5;                       //tamaño de cada dato BLOCK-1     
    private static byte[] RECORD = new byte[BLOCK];     //para la lectura de cada dato 
    private static double [][] AA=new double[N][N];
    private static LinkedList<Thread> hilos = new LinkedList<Thread>();	
    //---------------------------------------------------------------------------------------
	public static void Generar(){
		 WriteData();
         AsignarDatosMatriz();
	}
    //-----------------------------------------------------------------------------
    private static double ObtenerElemento() {
        String CAD;
            CAD = "";
            for(int i=0;i<BLOCK-1;i++) {
                CAD = CAD + (char)(RECORD[i]);
            }
            return Double.parseDouble(CAD);
        }
    //----------------------------------------------------------------------------
	public static void WriteData() {       //crea la data 
    double X;
    long num;
        try {
            FileWriter FW = new FileWriter(FILENAME);
            for(int i=1;i<=N*N;i++) {
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
    //--------------------------------------------------------------------------------------
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
                 AA[i][j] =CADENA;                   //asignar elemento a matriz , un double de 4 cifras
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
    //----Factorizacion LU via doolitle -----------------------------------------------------------
	public static Result Doolitle(double [][]M){        //hecho en clases de un curso anterior
		int n =M.length;
		double[][] L = new double[n][n];
		double[][] U = new double[n][n];
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				if(i==j){
					L[i][j]=1;                       // para que el algoritmo sea consistente
				}
			}
		}
		for(int i=0;i<n;i++){
			for(int j=0;j<i;j++){
				double suma = 0;
				for(int k=0;k<j;k++){
					suma+=(L[i][k]*U[k][j]);
				}
				L[i][j]=(M[i][j]-suma)/U[j][j];
			}
			for(int j=i;j<n;j++){
				double suma = 0;
				for(int k=0;k<i;k++){
					suma+=(L[i][k]*U[k][j]);
				}
				U[i][j]=(M[i][j]-suma)/L[i][i];
			}
		}
		return new Result(L, U);                        // encapsulando en Result (Result.L  Result.U)    
	}
	//--------------------------------------------------------------------------------------------------------------  	  
	public static Result LuBloques(double [][]M){
		 //cuatro for's para la asignacion de los 4 cudrantes , 		A11 A12
		 // este proceso se paraleliza , justamente con H=4 hilos		A21	A22	 
	     final Matrices m = new Matrices();
		 Thread h1 = new Thread(new Runnable(){
			public void run(){
				int n = M.length;
		 		int g=(int)n/2;
		 		m.A11 =new double[g][g];
		 		for(int i=0;i<g;i++){
					for(int j=0;j<g;j++){
						m.A11[i][j] = M[i][j];
				}
		 	}
			}
		 });						// para luego hacer doolitle(A11) se obtiene L11 y U11			
		 Thread h2 = new Thread(new Runnable() {
			public void run(){
				int n = M.length;
		 		int g=(int)n/2;
				m.A12=new double[g][g];
		 		for(int i=0;i<g;i++){
					for(int j=0;j<g;j++){
						m.A12[i][j] = M[i][j+g];
					}
		 		}
			}
		 });
		 Thread h3 = new Thread(new Runnable() {
			public void run(){
				int n = M.length;
		 		int g=(int)n/2;
				m.A21=new double[g][g];
		 		for(int i=0;i<g;i++){
					for(int j=0;j<g;j++){
						m.A21[i][j] = M[i+g][j];
			}
		 }
			}	
		 });
		 Thread h4 = new Thread(new Runnable() {
			public void run(){
				int n = M.length;
		 		int g=(int)n/2;
				m.A22=new double[(int)n/2][(int)n/2];
		 		for(int i=0;i<g;i++){
					for(int j=0;j<g;j++){
						m.A22[i][j] = M[i+g][j+g];
			}
		 }
			}	
		 });
		 h1.start();
		 h2.start();
		 h3.start();
		 h4.start();
		 hilos.add(h1);
		 hilos.add(h2);
		 hilos.add(h3);
		 hilos.add(h4);
		 //L11 y U11 de doolitle(A11) -> L11 U11 
		 //L11 U12 = A12 resolver:U12 = L11^-1*A12 y hallar U12 
		 //L21 U11 = A21  resolver: L21 =A21*U11^-1  hallar L21
		 //actualizacion schur A22 = A22 - L21*U12
		 //todas las operaciones se haran despues de que los hilos terminen
		 for (Thread hil:hilos) {
				try{
					hil.join();
				}catch(InterruptedException e){
					e.printStackTrace();
				}						
		}
		/* para matrices pequeñas: Imprimir(m.A11);Imprimir(m.A12);Imprimir(m.A21);Imprimir(m.A22);*/
		Result lu=Doolitle(m.A11);  					/*Imprimir(lu.L);Imprimir(lu.U);*/
		double[][] U12 = resolverU(lu.L,m.A12);			//L11*U12 = A12 resolver y hallar U12	
		double[][] L21=resolverL(lu.U,m.A21);			//L21*U11 = A21 resolver y hallar L21
		Producto L21U12 =new Producto(L21, U12);
		m.A22=resta(m.A22 , L21U12.PROD); 				//actualizacion schur A22 = A22 - L21*U12
		Result lu22= Doolitle(m.A22);					//doolitle(A22) podria implementar una recursividad 
		double[][]L= Lglobal(lu.L, L21, lu22.L);		//ensamblar los L(L12=0) y ensamblar los U(U21=0)								
		double[][]U= Uglobal(lu.U, U12, lu22.U);		//				lu.L   0				 lu.U	U12	
		return new Result(L, U);						//		 L  =						U =
	}													// 				L21   lu.22.L			   0	lu.22.U			 
	//-----------------------------------------------------------------------------------------------------
	public static double[][] Lglobal(double[][] L11, double[][] L21, double[][] L22) {
    	int n = L11.length + L21.length;
    	double[][] L = new double[n][n];
    	// Copiar L11
    	for (int i = 0; i < L11.length; i++) {
        	for (int j = 0; j < L11.length; j++) {
            	L[i][j] = L11[i][j];
       	 }
    	}
    	// Copiar L21
    	for (int i = 0; i < L21.length; i++) {
        	for (int j = 0; j < L21.length; j++) {
            	L[i + L11.length][j] = L21[i][j];
        	}
    	}
    	// Copiar L22
    	for (int i = 0; i < L22.length; i++) {
        	for (int j = 0; j < L22.length; j++) {
            	L[i + L11.length][j + L11.length] = L22[i][j];
        	}
    	}
    	return L;
	}
	//-------------------------------------------------------------------------------------------------------------
	public static double[][] Uglobal(double[][] U11, double[][] U12, double[][] U22) {
    	int n = U11.length + U12.length;
    	double[][] U = new double[n][n];
    	// Copiar U11
    	for (int i = 0; i < U11.length; i++) {
        	for (int j = 0; j < U11.length; j++) {
            	U[i][j] = U11[i][j];
        	}
    	}
    	// Copiar U12
    	for (int i = 0; i < U12.length; i++) {
        	for (int j = 0; j < U12.length; j++) {
            	U[i][j + U11.length] = U12[i][j];
        	}
    	}
    	// Copiar U22
    	for (int i = 0; i < U22.length; i++) {
        	for (int j = 0; j < U22.length; j++) {
            	U[i + U11.length][j + U11.length] = U22[i][j];
        	}
    	}
    	return U;
	}
	//---------------------------------------------------------------------------------------------------------------
	public static double[][] resta(double [][]A,double[][]B){
		int n = A.length;
		double [][] resta = new double[n][n]; 
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				resta[i][j] = A[i][j]-B[i][j];
			}
		}
		return resta;
	}
	//--------------------------------------------------------------------------------------------------------------
	public static double[][] resolverU(double[][]A,double[][]B){
		//se tiene L11 y A12 : L11*U12 	= A12  			hallar U12
		//							U12 = L11^-1 A12
		double[][] inv=inversa(A);
		Producto prod=new Producto(inv,B);
		return prod.PROD;
	}
	//--------------------------------------------------------------------------------------------------------------
	public static double[][] resolverL(double[][]A,double[][]B){
		//se tiene U11 y A21: 	L21*U11 = A21    		hallar L21
		//						L21     = A21 U11^-1	
		Producto prod=new Producto(B,inversa(A));
		return prod.PROD;
	}
	//----------------------------------------------------------------------------------------------------------
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
	//----------------------------------------------------------------------------------------------------
	public static double[][]inversa(double [][] mm){
		int n=mm.length;			
		double [][] I = new double[n][n];
		double [][] m=Copiar(mm);
		for(int i=0;i<n;i++){
			I[i][i] = 1;
		}
		for(int i=0;i<n;i++){
			if(m[i][i]==0){
				boolean cambiado = false;
				for(int j=i+1;j<n;j++){
					if(m[j][i]!=0){
						double[]temp = m[i];
						m[i]=m[j];
						m[j]=temp;
						temp=I[i];
						I[i]=I[j];
						I[j]=temp;
						cambiado=true;
						break;
					}
				}
			
			if(!cambiado){
				throw new ArithmeticException("no existe inversa");
			}
		}
		double pivot =m[i][i];
		for(int j=0;j<n;j++){
			m[i][j] /=pivot;
			I[i][j] /=pivot;
		}
		for(int j=0;j<n;j++){		//operaciones por filas gauss
			if(i!=j){
				double factor=m[j][i];
				for(int k=0;k<n;k++){
					m[j][k]-=factor*m[i][k];
					I[j][k]-=factor*I[i][k];		
				}
			}
		}	
	}
		return I;
	}
	//----------------------------------------------------------------------------------------------------------
	public static void Imprimir(double[][] m){
		int n=m.length;
        System.out.println();
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				System.out.printf("%12.2f",m[i][j]);
			}
			System.out.println();
		}
		System.out.println();
	}	 
//--------------------------------------------------------------------------------------------------
    public static void main(String[] args){             //descomentar los "//\\" para corroborar, matrices pequeñas
		Generar();
        long start = System.nanoTime();
		Result serial=Doolitle(AA);
        long endSerial = (System.nanoTime()-start)/1000000000;
		System.out.println("\nel tiempo serial es "+endSerial+" segundos");
        //\\Imprimir(AA);
		//\\Imprimir(serial.L);
		//\\Imprimir(serial.U);   
        //\\Producto LUs=new Producto(serial.L, serial.U); 
		//\\Imprimir(LUs.PROD);  							//corroborando que L@U  = A , visible para matrices pequeñas
		//\\System.out.println("\n");
		start = System.nanoTime();
        Result paralelo=LuBloques(AA);
        long endParalelo = (System.nanoTime()-start)/1000000000;
        //\\Producto LUp=new Producto(paralelo.L,paralelo.U);
		//\\Imprimir(paralelo.L);
		//\\Imprimir(paralelo.U);
		//\\Imprimir(LUp.PROD);    						//corroborando L@U = A distingible para matrices pequeñas
        System.out.println("\nel tiempo paralelo es "+endParalelo+" segundos");
		System.out.println("la diferencia de tiempos es "+(endSerial-endParalelo)+"segundos");
    }
}
//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
class Result{
	double[][]L;
	double[][]U;
	public Result(double[][] L,double[][] U){
		this.L = L;
		this.U = U;
	}
}	
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Matrices{
	public double[][] A11;
	public double[][] A12;
	public double[][] A21;
	public double[][] A22;
}
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
class Producto{
	double [][]A;
	double[][]B;
	double[][]PROD;
	private static LinkedList<Thread> hilos = new LinkedList<Thread>();
	public Producto(double[][]A,double[][]B){
		this.A=A;
		this.B=B;
		productoParalelo();
	}
	//---------------------------------------------------------
	public  void Imprimir(double[][] m){
		int n=m.length;
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				System.out.printf("%12.2f",m[i][j]);
			}
			System.out.println();
		}
		System.out.println();
	}
	//--------------------------------------------------------
	public void productoParalelo(){
		int n =A.length;
		PROD = new double[n][n];
		Thread hil1 = new Thread(new Runnable(){
			public void run(){    
        	for(int i=0;i<(int)n/2;i++){
        		for(int j=0;j<(int)n/2;j++){
        			for(int k=0;k<n;k++){
        				PROD[i][j]+=A[i][k]*B[k][j];
        				}
        			}
        		}
				}
			});
			hilos.add(hil1);
			hil1.start();
			Thread hil2 = new Thread(new Runnable(){
			public void run(){
        		for(int i=0;i<(int)n/2;i++){
        			for(int j=(int)n/2;j<n;j++){
        				for(int k=0;k<n;k++){
        					PROD[i][j]+=A[i][k]*B[k][j];
        				}
        			}
        		}
        	}
			});
			hilos.add(hil2);
			hil2.start();
			Thread hil3 = new Thread(new Runnable(){
			public void run(){
        		for(int i=(int)n/2;i<n;i++){
        			for(int j=0;j<(int)n/2;j++){
        				for(int k=0;k<n;k++){
        					PROD[i][j]+=A[i][k]*B[k][j];
        				}
        			}
        		}
			}
			});
			hilos.add(hil3);
			hil3.start();
			Thread hil4 = new Thread(new Runnable(){
			public void run(){
        		for(int i=(int)n/2;i<n;i++){
        			for(int j=(int)n/2;j<n;j++){
        				for(int k=0;k<n;k++){
        					PROD[i][j]+=A[i][k]*B[k][j];
        				}
        			}
        		}
			}
			});
			hilos.add(hil4);
			hil4.start();
			for (Thread hil:hilos) {
				try{
					hil.join();
				}catch(InterruptedException e){
					e.printStackTrace();
				}						
			}
	}
}