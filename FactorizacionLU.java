//import java.io.*;
//import javax.swing.*;
import java.util.*;
public class FactorizacionLU{
	static int N = 20;//tiene que ser par.
	static double [][] a;
	static double [][] l_s;
	static double [][] u_s;
	static double [] b;
	static double [] x;
	static double [] y;
    private static LinkedList<Thread> hilos = new LinkedList<Thread>();	
	public static void Generar(){
		Random rnd =new Random();
		a= new double[N][N];
		b=new double[N];
		l_s=new double[N][N];
		u_s=new double[N][N];
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				a[i][j] = rnd.nextDouble() * 100;
			}
			b[i] = rnd.nextDouble() * 100;
		}
	}
	
	public static Result Doolitle(double [][]M){
		int n =M.length;
		double[][] L11 = new double[n][n];
		double[][] U11 = new double[n][n];
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				if(i==j){
					L11[i][j]=1;
				}
			}
		}
		for(int i=0;i<n;i++){
			for(int j=0;j<i;j++){
				double suma = 0;
				for(int k=0;k<j;k++){
					suma+=(L11[i][k]*U11[k][j]);
				}
				L11[i][j]=(M[i][j]-suma)/U11[j][j];
			}
			for(int j=i;j<n;j++){
				double suma = 0;
				for(int k=0;k<i;k++){
					suma+=(L11[i][k]*U11[k][j]);
				}
				U11[i][j]=(M[i][j]-suma)/L11[i][i];
			}
		}
		return new Result(L11, U11);
	}
	  	  
	public static Result LuBloques(double [][]M){
		 //cuatro for para la asignacion , ademas se asume que son matrices cuadradas,,, rectangulares ni d cñ 
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
		 });
		 hilos.add(h1);
		 h1.start();
		 // lu para A11  = L11 U11 se tiene L11 U 11
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
		 hilos.add(h2);
		 h2.start();
		 //L11 se tiene del punto anterior 
		 //L11 U12 = A12 resolver y hallar U12
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
		 hilos.add(h3);
		 h3.start();
		 //se tiene U11
		 //L21 U11 = A21    hallar L21
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
		 hilos.add(h4);
		 h4.start();
		 // actualizacion schur
		 //todas las operaciones se haran despues de que los hilos terminen
		 for (Thread hil:hilos) {
				try{
					hil.join();
				}catch(InterruptedException e){
					e.printStackTrace();
				}						
		}
		/*Imprimir(m.A11);
		Imprimir(m.A12);
		Imprimir(m.A21);
		Imprimir(m.A22);*/
		Result lu=Doolitle(m.A11); //L11 se tiene del punto anterior  
		/*Imprimir(lu.L);
		Imprimir(lu.U);*/
		Producto U12 = resolverU(lu.L,m.A12);			//L11 U12 = A12 resolver y hallar U12
		//Imprimir(U12.PROD);	
		Producto L21=resolverL(lu.U,m.A21);				//se tiene U11 	//L21 U11 = A21    hallar L21
		//Imprimir(L21.PROD);
		Producto L21U12 =new Producto(L21.PROD, U12.PROD);
		m.A22=resta(m.A22 , L21U12.PROD); //actualizacion schur
		//Imprimir(m.A22);
		Result lu22= Doolitle(m.A22);//llamar a luBloques para el A22
		double[][]L= Lglobal(lu.L, L21.PROD, lu22.L);//ensamblar los L    y ensamblar los U								
		double[][]U= Uglobal(lu.U, U12.PROD, lu22.U);
		return new Result(L, U);
	}
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

	public static Producto resolverU(double[][]A,double[][]B){
		//L11 U12 = A12 resolver y hallar U12
		double[][] inv=inversa(A);
		Producto prod=new Producto(inv,B);
		return prod;
	}

	public static Producto resolverL(double[][]A,double[][]B){
		//se tiene U11 	//L21 U11 = A21    hallar L21
		double[][] inv=inversa(A);
		Producto prod=new Producto(B,inv);
		return prod;
	}

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
		for(int j=0;j<n;j++){
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

	public static void Imprimir(double[][] m){
		int n=m.length;
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				System.out.printf("%12.2f",m[i][j]);
			}
			System.out.println();
		}
		System.out.println();
	}
	public static void Imprimir(double[] x){
		int n=x.length;
		for(int i=0;i<n;i++){
			System.out.printf("%12.2f",x[i]);
		}
		System.out.println();		
	}
	public static double[][] AxB(double[][]A,double[][]B) {
        int n=A.length;
		double[][] C=new double[n][n];
        for(int i=0;i<n;i++) {
            for(int j=0;j<n;j++) {
                C[i][j] = 0;
                for(int k=0;k<n;k++) {
                    C[i][j] = C[i][j] + A[i][k]*B[k][j];
                }
            }
        }
        return C;
    }
//--------------------------------------------------------------------------------------------------
    public static void main(String[] args){
		Generar();
		Result serial=Doolitle(a);
		Imprimir(a);
		Imprimir(serial.L);
		Imprimir(serial.U);
		Imprimir(AxB(serial.L, serial.U));
		System.out.println();
		System.out.println();
		Result paralelo=LuBloques(a);
		Imprimir(paralelo.L);
		Imprimir(paralelo.U);
		Imprimir(AxB(paralelo.L, paralelo.U));
	}
}
//------------------------------------------------------------------------------------------------------
class Result{
	double[][]L;
	double[][]U;
	public Result(double[][] L,double[][] U){
		this.L = L;
		this.U = U;
	}
}	
class Matrices{
	public double[][] A11;
	public double[][] A12;
	public double[][] A21;
	public double[][] A22;
}

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