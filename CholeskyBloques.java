
import java.io.FileWriter;
import java.io.RandomAccessFile;
import java.io.IOException;
import java.util.LinkedList;
//==============================================================================================
public class CholeskyBloques {
    private static int BLOCK = 4; //tamaño de datos BLOCK -1=3, el ultimo 'digito' es " ",al generar la simetrica el tamaño se duplica
    private static byte [] RECORD = new byte[BLOCK];
    private static String FILENAMEMATRIZ = "DATACholeskyBloques.TXT";
    private static String FILENAMEX = "DATACholeskyX.TXT";
    private static int N = 12;
    private static String CADENA;
    private static double [][] A =new double[N][N];
    private static LinkedList<Thread> hilos = new LinkedList<>();
    private static int NB = 4;
    private static int k=N/NB;
    private static double[][][][] ABloques = new double[k][k][NB][NB];
    private static double[][][][] L = new double[k][k][NB][NB];
    //------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------
    
    //------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------
    public static void main(String[]args){
        //mientras no sea simetrica definida positiva
        /*boolean flag = false;
        while(!flag){
            GenerarMatriz(FILENAMEMATRIZ);
            //flag = esPositiva(A);
            flag= true;
        } */
        //GenerarMatriz(FILENAMEMATRIZ);
        //esPositiva();
        //Imprimir(A);
        //ClassQR qr = QR(A);
        //Imprimir(qr.Q);
        //Imprimir(qr.R);
        //Imprimir(ProdParalelo(qr.Q, qr.R));
        //double t = NormaFrobenius(A, qr.Q);
        //System.out.println(t);
        //System.out.println(BuscarPositiva(A));
        /* para buscar positiva dado que ya es simetrica:boolean flag =false;
        ClassPositiva pos=new ClassPositiva(A,true);
        while(!flag){
            GenerarMatriz(FILENAMEMATRIZ);
            pos = BuscarPositiva(A);
            flag = pos.flag;
        }
        A =pos.M;*/
        GenerarMatrizSimetrica(FILENAMEMATRIZ);
        Imprimir(A);
        double [][] G = CholeskySerial(A);
        Imprimir(G);
        //Imprimir(ProdParalelo(G, Transponer(G)));
        ObtenerABloques();
        for(int i=0;i<ABloques.length;i++){
            for(int j=0;j<ABloques[0].length;j++){
                Imprimir(ABloques[i][j]); //
            }
        }
        //Imprimir(A);
        //double[][] X=GenerarX(FILENAMEX);
        //double nTest=ProductoTriple(X, A);
        //System.out.println(nTest); se estaba testeando estos metodos
       // Choleskybloques(); falta revisar
    }
    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------
    public static void Choleskybloques(){
        ObtenerABloques();
        for(int j=0;j<k;j++){
            L[j][j] =CholeskySerial(ABloques[j][j]);
            Imprimir(L[j][j]);
            //trsm A21  A31   y   syrk A22=A22 -   A33,paralelizar
            for(int i=j+1;i<k;i++){
                System.out.println("abajo: " + i);
                L[i][j] =TRSM(ABloques[i][j],L[j][j]);
            }
            for(int i=j+1;i<k;i++){
                ABloques[i][i]=SYRK(ABloques[i][i],L[i][j]);
            }
        }
    }
    //-------------------------------------------------------------------------
    //--------------------------------------------------------------------------
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
    //-------------------------------------------------------------------------
    //------------------------------------------------------------------------
    public static double[][]MatrizInversa(double [][] mm){
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
    //-----------------------------------------------------------------------
    //-------------------------------------------------------------------------
    public static double[][] TRSM(double[][]M,double[][]P){
        double[][]INVP =MatrizInversa(P);
        //double [][] SOL =new double[M.length][INVP[0].length];
        return ProdParalelo(M, INVP); 
    }
    //---------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    public static double[][]SYRK(double[][]M,double[][]P){
        return Resta(M,ProdParalelo(P, Transponer(P)));
    }
    //------------------------------------------------------------------------
    //-------------------------------------------------------------------------
    public static double[][] Resta(double [][]A,double[][]B){
		int n = A.length;
		double [][] resta = new double[n][n]; 
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				resta[i][j] = A[i][j]-B[i][j];
			}
		}
		return resta;
	}
    //-----------------------------------------------------------------------
    //--------------------------------------------------------------------
    public static void ObtenerABloques(){
        int k =N/NB;
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < k; j++) {
                ABloques[i][j] = new double[NB][NB];  // Asignando submatrices de tamaño NB x NB
            }
        }
        System.out.println("dentro luego de ABloques forasin k ="+ k );
        //dividir en k bloques de tamaño NB
        /*for(int columna=0;columna<k;columna++){
            System.out.println("columna: "+columna);
            for(int abajo=columna;abajo<k;abajo++){
                System.out.println("abajo: "+abajo);
                for(int i=abajo*NB;i<(abajo+1)*NB;i++){
                    System.out.println("i: "+i);
                    System.arraycopy(A[i], columna*NB, ABloques[abajo][columna][i-abajo*NB],0, NB);
                }
                Imprimir(ABloques[abajo][columna]);
            }
        }*/
            
        Thread [] columnas =new Thread[k];
        for (int hil=0;hil<k;hil++) {
            final int  hilo = hil;
            columnas[hilo]=new Thread(new Runnable(){
                public void run(){
                    for(int abajo=hilo;abajo<k;abajo++){
                        for(int i=abajo*NB;i<(abajo+1)*NB;i++){
                            System.arraycopy(A[i], hilo*NB, ABloques[abajo][hilo][i-abajo*NB], 0, NB);
                        }
                        //Imprimir(ABloques[abajo][hilo]);
                    }
                }
            });
        }
        for (Thread thread : columnas) {
            thread.start();
        }
        for (Thread thread : columnas) {
            try{
                thread.join();
            }catch(InterruptedException e){
                Thread.currentThread().interrupt();
            }
        }
        /*double [][] A11=new double[NB][NB];
            for(int i=(1-1)*NB;i<1*NB;i++){
                System.arraycopy(A[i], 0, A11[i], 0, NB);
            }
            double [][] A21=new double[NB][NB];
            for(int i=(2-1)*NB;i<2*NB;i++){
                System.arraycopy(A[i], 0, A21[i], 0, NB);
            }
            double [][] A31=new double[NB][NB];
            for(int i=(3-1)*NB;i<3*NB;i++){
                System.arraycopy(A[i], 0, A31[i], 0, NB);
        }*/
        
    }
    //-----------------------------------------------------------------------------
    //---------------------------------------------------------------------------------
    public static ClassPositiva BuscarPositiva(double[][]M){
        ClassQR qr;
        double[][] M1 = new double[M.length][M[0].length];
        double epsilon=10000;
        int iter=0;
        while(iter<100&&epsilon>0.000001){
            qr = QR(M);
            M1 = ProdParalelo(qr.R, qr.Q);
            //norma de probenius o infinita M1 - M
            epsilon = NormaFrobenius(M,M1);
            //Imprimir(M1);
            System.out.println(epsilon);
            M = M1;
            iter++;
        }
        System.out.println("dentro de BuscarPositiva");
        Imprimir(M);
        //los elemetos de la diagonal son positivos?
        boolean flag=true;
        for(int i=0;i<M.length;i++){
            if(M[i][i]<0){
                flag =false;
            }
        }  
        return new ClassPositiva(M, flag);
    }
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    public static double NormaFrobenius(double[][]M,double[][]P){
        double norma = 0;
        for(int i=0;i<M.length;i++){
            double suma =0;
            for(int j=0;j<M[0].length;j++){
                if(i!=j){
                    suma = Math.pow(P[i][j] -M[i][j],2);
                }
            }
            norma +=suma;
        }
        return Math.sqrt(norma);
    }
    //---------------------------------------------------------------------------
    //----------------------------------------------------------------------------
    public static double[][] ProdParalelo(double[][]M1,double[][]M2){
        double [][]PROD=new double[M1.length][M2[0].length];
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
                    suma+=M1[i][k]*M2[k][j];
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
                        suma+=M1[i][k]*M2[k][j];
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
                        suma+=M1[i][k]*M2[k][j];
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
                        suma+=M1[i][k]*M2[k][j];
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
    return PROD;
    }
    //------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------
    public static ClassQR QR(double[][]M){
        int filas=M.length;
        int columnas = M[0].length;
        double x;
        double [][]Q =new double[filas][columnas];
        double[][]R=new double[filas][columnas];
        for (int i = 0; i < columnas; i++) {
            double res=0;
            for(int k=0;k<filas;k++){    //x =Math.sqrt(A.prodEsc(i, i));
                res+=(M[k][i]* M[k][i]);
            }
            x = Math.sqrt(res);
            R[i][i] =x;
            for (int k = 0; k < filas; k++) {
                    x = M[k][i] / R[i][i];
                    Q[k][i] = x;
            }
            for (int j = i + 1; j < columnas; j++) {//filas?gpt
                    //primero completar los R[i][j]
                    double aux=0;
                    for(int m=0;m<columnas;m++){
                            aux+=(Q[m][i]*M[m][j]);
                    }
                    R[i][j]=aux;
                    for (int k = 0; k < filas; k++) {
                            x = M[k][j] - R[i][j] * Q[k][i];
                            M[k] [j]= x;
                    }
            }
        }
        return new ClassQR(Q, R);
    }
    //-------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------
    public static double[][] CholeskySerial(double[][] M){
        int n = M.length;
        double [][]G =new double[M.length][M.length];
        //Imprimir(G);
        for(int j =0;j<n;j++){
            double suma = 0;
            for(int k =0;k<j;k++){
                suma+=(G[j][k]*G[j][k]);
            }
            //System.out.printf("\nM[j][j]-suma:%1.2f\n",M[j][j]-suma);
            G[j][j] = Math.pow(M[j][j]-suma,0.5);
            for(int i =j+1;i<n;i++){
                double sumai = 0;
                for(int m =0;m<j;m++){
                    sumai+=(G[i][m]*G[j][m]);       
                }
                G[i][j] = (A[i][j]-sumai)/G[j][j];
            }
        }
        return G;
    }   
    //-------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------
   /* def cholesky(a,g):
    n = len(a)
    for j in range(n):
        suma = 0
        for k in range(j):
            suma = suma + g[j,k]**2
        g[j,j] = mat.pow(a[j,j]-suma,0.5)
        for i in range(j+1,n):
            suma = 0
            for k in range(j):
                suma = suma + g[i,k]*g[j,k]
            g[i,j] = (a[i,j]-suma)/g[j,j]*/
    //-------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------
    public static void WriteDataMatriz(String FILENAME){
        double num;
        try{
            FileWriter fw = new FileWriter(FILENAME);
            //generando el aletario
            double max= Math.pow(10,BLOCK-1);
            double min = Math.pow(10,BLOCK-2);
            //int l = (int)(N*(N+1)/2);   //para una matriz simetrica
            for (int i=0;i<N*N;i++){
                num = Math.random()*(max - min ) + min;
                fw.write((long)num + " ");
            }
            fw.close();
        }catch(IOException e){
            e.printStackTrace();
        }
    }
    //----------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------
    public static void GenerarMatriz(String FILENAME){
        WriteDataMatriz(FILENAME);
        int n,k,i=0,j=0;
        double num;
        try{
            RandomAccessFile RAF = new RandomAccessFile(FILENAME,"r");
            n = (int)RAF.length()/BLOCK;
            k=0;
            //System.out.println("n: "+ n +" k"+ k);
            while(k<n){
                RAF.seek(k*BLOCK);      //posicionando el puntero
                RAF.read(RECORD);       //almacenando den RECORD "123 ", "432 "....
                num = convertir();      //llamando a convertir
                //almacenar en matriz
                A[i][j]=A[j][i]= num;      //pues sera simetrica
                if(j == i){           // j hasta j==i (la diagonal)
                    j = 0;
                    i++;
                    if(i ==N){
                        break;  //fila N=4 , no posible
                    }
                }else{
                    j++;
                } 
                k++;
            }
            RAF.close();
        }catch(IOException e){
            e.printStackTrace();
        }
       //Imprimir(A);
    }
    //-----------------------------------------------------------------------
    //----------------------------------------------------------------------
    public static void GenerarMatrizSimetrica(String FILENAME){
        WriteDataMatriz(FILENAME);
        int n,k,i=0,j=0;
        double num;
        try{
            RandomAccessFile RAF = new RandomAccessFile(FILENAME,"r");
            n = (int)RAF.length()/BLOCK;
            k=0;
            //System.out.println("n: "+ n +" k"+ k);
            while(k<n){
                RAF.seek(k*BLOCK);      //posicionando el puntero
                RAF.read(RECORD);       //almacenando den RECORD "123 ", "432 "....
                num = convertir();      //llamando a convertir
                //almacenar en matriz
                A[i][j]= num;      //pues sera simetrica
                if(j == N-1){           // j hasta j==i (la diagonal)
                    j = 0;
                    i++;
                }else{
                    j++;
                } 
                k++;
            }
            RAF.close();
        }catch(IOException e){
            e.printStackTrace();
        }
       A=ProdParalelo(A, Transponer(A));
    }
    //-----------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------
    public static double convertir(){
        String CAD = " ";
        for(int i=0;i<BLOCK-1;i++){
            CAD = CAD + (char)RECORD[i];
        }
        return Double.parseDouble(CAD);
    }
    //-----------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------
    public static synchronized void Imprimir(double[][] M){
        System.out.println();
        for(int i=0;i<M.length;i++){
            for(int j=0;j<M[0].length;j++){
                System.out.printf("%1.1f ",M[i][j]);
            }
            System.out.println();
        }
        System.out.println();
    }
    //------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------
    public static void WriteDataX(String FILENAME){
        double num;
        try{
            FileWriter fw = new FileWriter(FILENAME);
            for(int i=0;i<N;i++){
                double max = Math.pow(10, BLOCK-1);
                double min = Math.pow(10,BLOCK-2);
                num = Math.random()*(max - min)+min;
                fw.write((long)num + " ");
            }
            fw.close();
        }catch(IOException e){
            System.out.println(e.getMessage());
        }
    }
    
    public static double[][] GenerarX(String FILENAME){
        WriteDataX(FILENAME);
        double num;
        int k,n;
        double [][]X =new double[1][N];
        try{
            RandomAccessFile RAF = new RandomAccessFile(FILENAME, "r");
            n=(int)RAF.length()/BLOCK;
            k = 0;
            while(k<n){
                RAF.seek(k*BLOCK);
                RAF.read(RECORD);
                num = convertir();
                X[0][k] =num;
                k++;                
            }
            RAF.close();
        }catch(IOException e){
            System.out.println(e.getMessage());
        }
        //Imprimir(X);
        return X;
    }
    //--------------------------------------------------------------------
    //-------------------------------------------------------------------------
    public static double [][] Transponer(double[][]M){
        double[][]T=new double  [M[0].length][M.length];
        for(int i=0;i<M.length;i++){
            for(int j=0;j<M[0].length;j++){
                T[j][i]=M[i][j];
            }
        }
        return T;
    }
    //---------------------------------------------------------------------------
    //-----------------------------------------------------------------------
    public static double ProductoTriple(double[][]X,double[][]M){
        double [][]Xt = Transponer(X);
        double[][]PROD=new double[X.length][M[0].length];
        for(int i=0;i<X.length;i++){
            for(int j=0;j<M[0].length;j++){
                double sum=0;
                for(int k=0;k<M.length;k++){
                    sum += X[i][k]*M[k][j];
                }
                PROD[i][j] = sum;
            }
        }
        double XAXt=0;
        for(int i=0;i<PROD.length;i++){
            for(int j=0;j<Xt[0].length;j++){
                double sum=0;
                for(int k=0;k<Xt.length;k++){
                    sum += X[i][k]*M[k][j];
                }
                XAXt = sum;
            }
        }
        return XAXt;
    }
    //----------------------------------------------------------------------
    //--------------------------------------------------------------------------
    
    //-------------------------------------------------------------------------
    //----------------------------------------------------------------------
    public static void esPositiva(){
        boolean flag=false;
        while(!flag){
            int n=50;
            GenerarMatriz(FILENAMEMATRIZ);
            for(int  i=0;i<n;i++){
                double [][]x = GenerarX(FILENAMEX);
                double p = ProductoTriple(x, A);
                if(p<0){
                    flag = false;
                }else{
                    flag = true;
                }
            }
        }
    }
    //--------------------------------------------------------------------------------
    //------------------------------------------------------------------
}
class ClassQR{
    double[][]Q;
    double[][]R;
    ClassQR(double[][]q,double[][]r){
        this.Q=q;
        this.R=r;
    }
}
//------------------------------------------------------------------------------
class ClassPositiva{
    double[][]M;
    boolean flag;
    ClassPositiva(double[][]m,boolean f){
        this.M=m;
        this.flag=f;
    }
}