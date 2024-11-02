import java.util.LinkedList;
import javax.swing.*;
//------------------------ejecutar DataSet para crear primero los datos correspondientes-----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------
public class PC1Preg1PotenciaMatriz{
    private static LinkedList<Thread> hilos = new LinkedList<Thread>();
	private static double PROD [][];
	private static double[][] a1;
	private static int N;
	//--DataSet crea y almacena la matriz [n=1000][n=1000] en DATAPC1Preg1PotenciaMatriz.TXT
	//--Se realiza el calculo de la potencia p-esima de modo serial en PotenciaSerial()
	//--Se realiza el calculo de la potencia p-esima de modo paralelo en main() 
	//--la diferencia es distingible para matrices [1000][1000] por ejemplo para p=3 potencia
	//--los tiempos se muestran en segundos
	//--los resultados de las potencias se guardan en POTENCIASERIAL.TXT   y POTENCIAPARALELA.TXT
    //--se podria realizar la prueba con matrices pequeñas comparando los resultados con numpy en python por ejemplo(se hizo para una matriz pequeña)
	//--------------------------------------------------------------------------------------------
	public static void PotenciaSerial(Matrix a,int potencia){ 
		int n=a.getRows();
		double [][] prod=new double[n][n]; 
		long Start = System.nanoTime();
        Asignar(prod, a);
		for(int i=0;i<potencia-1;i++){
			prod=ProdSer(a,prod);
		}
		long End=(System.nanoTime()-Start)/1000000;
		System.out.printf("\nmatriz %dx%d\tpotencia %d-esima\n",n,n,potencia);
		System.out.printf("el tiempo serial es %d milisegundos\n",End);
		Matrix M=new Matrix(prod);
        DataSet.WriteFile(M, "POTENCIASERIAL.TXT"); 		
	}
	//-------------------------------------------------------------------------------
    public static void imprimir(double[][]m){
		for(int i=0;i<m.length;i++){
			for(int j=0;j<m[0].length;j++){
				System.out.printf("%1.2f\t",m[i][j]);
			}
			System.out.println();
		}
	}
	//---------------------------------------------------------------------------------
    public static void Asignar(double[][] mm,Matrix m){
		int n = m.getRows();
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				mm[i][j]=m.GetCell(i,j);
			}
		}
	}
	//----------------------------------------------------------------------------------
	public static void Asignar(double[][] mm,double[][] m){
		int n = m.length;
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				mm[i][j]=m[i][j];
			}
		}
	}
	//-----------------------------------------------------------------------------------
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
	//----------------------------------------------------------------------------------
    public static void main(String[]args){
        a1 = DataSet.ReadFile(DataSet.filas, DataSet.columnas);
        Matrix a = new Matrix(a1);
        int potencia;
        potencia=Integer.parseInt(JOptionPane.showInputDialog(null," ingresa la potencia"));
        PotenciaSerial(a,potencia);		//llama a proceso serial
        N = a1.length;
        PROD= new double[N][N];
        double[][] b = new double[N][N];
        Asignar(PROD,a1);
		long Start = System.nanoTime();	//proceso paralelo
		for(int p=0;p<potencia-1;p++){
			hilos.clear();
			Asignar(b,PROD);			
			double[][] temp1;
			double[][] temp2;
			double[][] temp3;			  // este analisis maneja si N es par o impar
			double[][] temp4;			  // mitad = (int)n/2 + r={0,1}   y esto en los indices pero era menos legible
			int mitad = (int)N/2 + (N%2); //para n par n/2    para  n impar n/2 + 1
			if(N%2==0){					  // mitad =2  si n=4    mitad=3 si n=5
				temp1 = new double[mitad][mitad];
				temp2 = new double[mitad][mitad];
				temp3 = new double[mitad][mitad];
				temp4 = new double[mitad][mitad];
			}
			else{									//ejemplo n=5
				temp1 = new double[mitad][mitad];     // [3][3] superior izq
				temp2 = new double[mitad][mitad-1];	  // [3][2] superior der	
				temp3 = new double[mitad-1][mitad];	  // [2][3] inferior izq
				temp4 = new double[mitad-1][mitad-1]; // [2][2] inferior der
			}
			Thread hil1 = new Thread(new Runnable(){
			public void run(){    
        	for(int i=0;i<mitad;i++){
        		for(int j=0;j<mitad;j++){
        			double suma=0;
					for(int k=0;k<N;k++){
        				suma+=a1[i][k]*b[k][j];
        			}
					temp1[i][j]=suma;
        			}
        		}
				}
			});
			hilos.add(hil1);
			Thread hil2 = new Thread(new Runnable(){
			public void run(){
        		for(int i=0;i<mitad;i++){
        			for(int j=mitad;j<N;j++){
        				double suma=0;
						for(int k=0;k<N;k++){
        					suma+=a1[i][k]*b[k][j];
        				}
						temp2[i][j-mitad]=suma;
        			}
        		}
        	}
			});
			hilos.add(hil2);
			Thread hil3 = new Thread(new Runnable(){
			public void run(){
        		for(int i=mitad;i<N;i++){
        			for(int j=0;j<mitad;j++){
        				double suma=0;
						for(int k=0;k<N;k++){
        					suma+=a1[i][k]*b[k][j];
        				}
						temp3[i-mitad][j]=suma;
        			}
        		}
			}
			});
			hilos.add(hil3);
			Thread hil4 = new Thread(new Runnable(){
			public void run(){
        		for(int i=mitad;i<N;i++){
        			for(int j=mitad;j<N;j++){
        				double suma=0;
						for(int k=0;k<N;k++){
        					suma+=a1[i][k]*b[k][j];
        				}
						temp4[i-mitad][j-mitad]=suma;
        			}
        		}
			}
			});
			hilos.add(hil4);
			hil1.start();
			hil2.start();
			hil3.start();
			hil4.start();
			for (Thread hil:hilos) {
				try{
					hil.join();
				}catch(InterruptedException e){
					e.printStackTrace();
				}						
			}			//ahora reemzamblando a PROD, java proporciona este metodo corroborado
			for (int i = 0; i < mitad; i++) {										//ejemplo 	
				System.arraycopy(temp1[i], 0, PROD[i], 0, mitad); //N par n=4 mitad=2
				System.arraycopy(temp2[i], 0, PROD[i], mitad, N/2);		//N=5 impar mitad=3
			}																	//[3][3]    [3][2]
			for (int i = 0; i < N / 2; i++) {									//[2][3]	[2][2]	
				System.arraycopy(temp3[i], 0, PROD[i + mitad], 0, mitad);
				System.arraycopy(temp4[i], 0, PROD[i + mitad], mitad, N / 2);
			}
		}
		long end = (System.nanoTime()-Start)/1000000;
		System.out.println("el tiempo paralelo es "+ end + " milisegundos");
        Matrix M =new Matrix(PROD);
        DataSet.WriteFile(M, "POTENCIAPARALELA.TXT");	
	}
}
