//import java.io.*;
import javax.swing.*;
import java.util.*;
public class ProductoParalelo_4{
	private static double[][] a;
	private static double[][] b;
	private static double [][] prod;
	private static int n;
	private static boolean flag1 = false;
	private static boolean flag2 = false;
	private static boolean flag3 = false;
	private static boolean flag4 = false;
	private static LinkedList<Thread> hilos = new LinkedList<Thread>();
	public static void main(String[] args){
		System.out.println("el tiempo es  milisegundos");
		a = DataSet.ReadFile(DataSet.filas, DataSet.columnas);
        int n = a.length;
        int potencia;
        potencia=Integer.parseInt(JOptionPane.showInputDialog(null," ingresa la potencia"));
        prod= new double[n][n];
        b = new double[n][n];
        Asignar(b,a);
        //imprimir(a);
		long Start = System.nanoTime();
		for(int p=0;p<potencia;p++){
			if(p>0){
				Asignar(b,prod);
			}
			Thread hil1 = new Thread(new Runnable(){
			public void run(){    
        	for(int i=0;i<(int)n/2;i++){
        		for(int j=0;j<(int)n/2;j++){
        			for(int k=0;k<n;k++){
        				prod[i][j]+=a[i][k]*b[k][j];
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
        					prod[i][j]+=a[i][k]*b[k][j];
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
        					prod[i][j]+=a[i][k]*b[k][j];
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
        					prod[i][j]+=a[i][k]*b[k][j];
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
		long end = (System.nanoTime()-Start)/1000000000;
		System.out.println("el tiempo paralelo es "+ end + " segundos");
		//imprimir(prod);		
	}
	
	public static void Control(){
		if((flag1 && flag2 && flag3 && flag4)==true){
			Asignar(b,prod);
		}
	}
	public static void pro(int fil_begin, int fil_end, int col_begin,int col_end){
		for(int i=fil_begin;i<fil_end;i++){
			for(int j=col_begin;j<col_end;j++){
				for(int k=0;k<n;k++){
					prod[i][j] += a[i][k]*b[k][j]; 
				}
			}
		}
	}
	public static void imprimir(double[][]m){
		for(int i=0;i<m.length;i++){
			for(int j=0;j<m[0].length;j++){
				System.out.printf("%1.2f\t",m[i][j]);
			}
			System.out.println();
		}
		System.out.println();
	}

	public static void Asignar(double[][] mm,double m[][]){
		int n = m.length;
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				mm[i][j]=m[i][j];
			}
		}
	}

	public static void Asignar0(double m[][]){
		int n = m.length;
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				m[i][j]=0.0;
			}
		}
	}
}
 
//