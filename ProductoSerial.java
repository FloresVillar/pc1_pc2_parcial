import java.io.*;
import javax.swing.*;
public class ProductoSerial{
	public static void Metodo(){
		//long Time1,Time2;
        //double x;
        double[][] a1;
        a1 = DataSet.ReadFile(DataSet.filas, DataSet.columnas);
        Matrix a = new Matrix(a1);
        int potencia;
        potencia=Integer.parseInt(JOptionPane.showInputDialog(null," ingresa la potencia"));
		int n=a.getRows();
		double [][] prod=new double[n][n];
		double[][] b=new double[n][n];
		Asignar(b,a);
		long Start = System.nanoTime();
		for(int i=0;i<potencia;i++){
			if(i>0){
				Asignar(b,prod);
			}
			prod=ProdSer(a,b);
		}
		long End=(System.nanoTime()-Start)/1000000000;
		System.out.printf("el tiempo es %d segundos\n",End);
		//a.imprimir();
		//imprimir(prod);		
	}
	public static void main(String[]args){
		Metodo();	
	}
	public static void imprimir(double[][]m){
		for(int i=0;i<m.length;i++){
			for(int j=0;j<m[0].length;j++){
				System.out.printf("%1.2f\t",m[i][j]);
			}
			System.out.println();
		}
	}
	
	public static void Asignar(double[][] mm,Matrix m){
		int n = m.getRows();
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				mm[i][j]=m.GetCell(i,j);
			}
		}
	}
	public static void Asignar(double[][] mm,double[][] m){
		int n = m.length;
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				mm[i][j]=m[i][j];
			}
		}
	}
	public static double[][]ProdSer( Matrix m1,double[][] m2){
		int n = m1.getRows();
		double PROD[][]=new double[n][n];
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				PROD[i][j]=0;
				for(int k=0;k<n;k++){
					PROD[i][j]+=m1.GetCell(i,k)*m2[k][j];
				}
			}
		}
		return PROD;
	}
}