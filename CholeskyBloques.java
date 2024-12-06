
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.RandomAccessFile;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.LinkedList;
import java.util.Random;
import java.util.Scanner;
import java.io.IOException;
import java.io.FileNotFoundException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.io.File;
import java.io.BufferedReader;
//==============================================================================================
public class CholeskyBloques {
    private static int BLOCK = 5; //tamaño de datos BLOCK -1=3, el ultimo 'digito' es " ",al generar la simetrica el tamaño se duplica
    private static byte [] RECORD = new byte[BLOCK];
    private static String FILENAMEMATRIZ = "DATACholeskyBloques.TXT";
    private static String FILENAMEX = "DATACholeskyX.TXT";
    private static int N_global = 12;
    private static String CADENA;
    private static double [][] A_Global =new double[N_global][N_global];
    private static LinkedList<Thread> hilos = new LinkedList<>();
    private static int NB = 4;
    private static int k_global=N_global/NB;
    private static double[][][][] ABloques = new double[k_global][k_global][NB][NB];
    private static double[][][][] LBloques = new double[k_global][k_global][NB][NB];
    //------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------
    public static void main(String[]args){
        //GenerarMatriz(FILENAMEMATRIZ);
        //esPositiva();
        //Imprimir(A);   
         //TESTEANDO Clase QR para determinar si una matriz es definida positiva via A'=RQ
        /*
        GenerarMatriz(FILENAMEMATRIZ);
        ClassQR qr = QR(A);
        Imprimir(qr.Q);
        Imprimir(qr.R);
        Imprimir(ProdParalelo(qr.Q, qr.R));
        double t = NormaFrobenius(A, qr.Q);
        System.out.println(t);
        System.out.println(BuscarPositiva(A));
        //para buscar positiva dado que ya es simetrica:boolean flag =false;
        ClassPositiva positiva=new ClassPositiva(A,true);
        while(!flag){
            GenerarMatriz(FILENAMEMATRIZ);
            positiva = BuscarPositiva(A);
            flag = pos.flag;
        }
        A =positiva.M;
        */
        /*  //En lugar de usar QR se hizo A= AA^t y se tiene una simetrica
        GenerarMatrizSimetrica(FILENAMEMATRIZ);
        Imprimir(A);
        double [][] G = CholeskySerial(A);
        Imprimir(G);
            //TESTEANDO si cholesky serial realemente descompuso A en G (triangular inferior)
        /*
        Imprimir(ProdParalelo(G, Transponer(G)));
        */
        // TESTEANDO la obtencion de los bloques para el tratamiento paralelo, se otendran solo desde la diagonal hacia↓
        // es decir (k)(k+1)/2  BLOQUES
        /*
        ObtenerABloques();
        for(int i=0;i<ABloques.length;i++){
            for(int j=0;j<ABloques[0].length;j++){
                
                Imprimir(ABloques[i][j]); //
            }
        } 
        */
        //TESTEANDO DataSetCholesky
        /* 
        DataSetCholesky data = new DataSetCholesky(FILENAMEMATRIZ, N,N,BLOCK,"FileWriter");//2modos de CrearData PrintWriter(Scanner) o FileWriter(RandomAccessFile)
        A = DataSetCholesky.ReadDataRAF(FILENAMEMATRIZ);
        DataSetCholesky.WriteDataRAF("DATACholeskyBloquesCOPIA.TXT");
        Imprimir(DataSetCholesky.ReadDataRAF("DATACholeskyBloquesCOPIA.TXT"));
        Imprimir(A);
        */
        //Imprimir(A);
        //double[][] X=GenerarX(FILENAMEX);
        //double nTest=ProductoTriple(X, A);
        //System.out.println(nTest); se estaba testeando estos metodos
       // TESTEANDO CholeskyBloques        falta revisar
       
       /*
       for(int i=0;i<LBloques.length;i++){
            for(int j=0;j<LBloques[0].length;j++){
                Imprimir(LBloques[i][j]); //
            }
        }  
        */
            // TESTEANDO las triangulares 
        /* 
        ObtenerABloques();
        double [][] L11= CholeskySerial(ABloques[0][0]);
        Imprimir(L11);
        ImprimirInversa(InversaTriangularInferior(L11));
        Imprimir(ProdParalelo(L11, InversaTriangularInferior(L11)));
        Imprimir(Transponer(L11));
        ImprimirInversa(InversaTriangularSuperior(Transponer(L11)));
        Imprimir(ProdParalelo(Transponer(L11), InversaTriangularSuperior(Transponer(L11))));
        */
        /*     FUNCIONA OBTENIENDO BLOQUES AFUERA
        GenerarMatrizSimetrica(FILENAMEMATRIZ);
        Imprimir("A",A_Global);
        ObtenerABloques();
       for(int i=0;i<ABloques.length;i++){
         for(int j=0;j<ABloques[0].length;j++){
                
            Imprimir("ABloques",ABloques[i][j]); 
        }
        }
    */
    /* 
        DataSetCholesky dataGlobal = new DataSetCholesky("DATASETCHOLESKY.TXT", 12, 12, 3, "FileWriter");
        A_Global = dataGlobal.ReadDataRAF();
        Imprimir("A con DataSetCholseky", A_Global);
        ObtenerABloques();
        for(int i=0;i<ABloques.length;i++){
            for(int j=0;j<ABloques[0].length;j++){
                Imprimir("ABloques",ABloques[i][j]); 
            }
        }
        Imprimir("ABloques[0][0]", ABloques[0][0]);
        double [][] L11= CholeskySerial(ABloques[0][0]);
        Imprimir("L11",L11);
        double[][]L21=TRSM(ABloques[1][0], L11);
        Imprimir("L21 = TRSM(ABloques[2][1] L11)",L21);
        Imprimir("ABloques[3][1]", ABloques[2][0]);
        double[][]L31=TRSM(ABloques[2][0], L11);
        Imprimir("L31 =TRSM ABloques[3][1]", L31);
        Imprimir("ABloques[2][2]",ABloques[1][1]);
        double[][]A22p=SYRK(ABloques[1][1], L21,L21);
        Imprimir("SYRK A22=A22 -L21L21^t",A22p);
        double[][]A32p=SYRK(ABloques[2][1], L31, L21);
        Imprimir("SYRK A32=A32 -L31 L21^t", A32p);
        double[][]A33p =SYRK(ABloques[2][2],L31,L31);
        Imprimir("SYRK A33 =A33 -L31 L31^t", A33p);
        double[][]L22=CholeskySerial(A22p);
        Imprimir("L22 ", L22);
        double[][]L32=TRSM(A32p, L22);
        Imprimir("L32", L32);
        double[][]A33pp =SYRK(A33p, L32,L32);
        Imprimir("A33pp", A33pp);
        double[][]L33 = CholeskySerial(A33pp);
        Imprimir("L33", L33);
        double [][]LBloquesF =new double[N_global][N_global];
        CopiarBloque(LBloquesF, L11, 0, 0);  // L11 -> posición (0, 0)
        CopiarBloque(LBloquesF, L21, 4, 0);  // L21 -> posición (4, 0)
        CopiarBloque(LBloquesF, L31, 8, 0);  // L31 -> posición (8, 0)
        CopiarBloque(LBloquesF, L22, 4, 4);  // L22 -> posición (4, 4)
        CopiarBloque(LBloquesF, L32, 8, 4);  // L32 -> posición (8, 4)
        CopiarBloque(LBloquesF, L33, 8, 8);
        Imprimir("LBloques reemzamblado", LBloquesF);
        //corroborando con choleskySerial
        double[][]ASerial = Copiar(A_Global);
        Imprimir("factorizacion cholesky Serial",CholeskySerial(ASerial));
    */        
        //LOS RESULTADOS INDICAN QUE TODO ES CORRECTO LO CUAL ES GENIAL GENIAL GENIAL !!!
            //EL TESTEO ANTERIOR SE HIZO DE FORMA NO ITERADA PASO A PASO y RESULTA EN LA OBTENCION DE LOS LBLOQUES
            //AHORA SE INTENTA LA FORMA ITERADA y RESULTA UN ERROR en el manejo de los indices seguramente FALTA CORREGIR ELLO
            //UNICAMENTE SE REQUIERE CORREGIR EL MANEJO DE LOS INDICES DENTRO DE LOS FOR DE Choleskybloques()
    
        DataSetCholesky dataGlobal = new DataSetCholesky("DATASETCHOLESKY.TXT", 12, 12, 3, "FileWriter");
        A_Global = dataGlobal.ReadDataRAF();
        Imprimir("A con DataSetCholseky", A_Global);
        ObtenerABloques();
       for(int i=0;i<ABloques.length;i++){
            for(int j=0;j<ABloques[0].length;j++){
                Imprimir("ABloques",ABloques[i][j]); 
            }
        } 
        Choleskybloques();
        for(int i=0;i<ABloques.length;i++){
            for(int j=0;j<ABloques[0].length;j++){
                Imprimir("LBloques",LBloques[i][j]); 
            }
        }

    }
    //-------------------------------------------------------------------------------
    public static void ObtenerLGlobal(){
        for(int k=0;k<k_global;k++){
           for(int abajo=k*NB;abajo<k_global;abajo++){
                for(int filas=abajo;filas<(abajo+1)*NB;filas++){
                    
                }
           } 
        }
    }
    //-------------------------------------------------------------------------------
    public static void Choleskybloques(){
        //ObtenerABloques();
        /*
        for(int i=0;i<ABloques.length;i++){
            for(int j=0;j<ABloques[0].length;j++){
                Imprimir("ABloques",ABloques[i][j]); 
            }
        }
        */
        for(int j=0;j<k_global;j++){
            System.out.println("k = "+(j+1));
            LBloques[j][j] =CholeskySerial(ABloques[j][j]);
            Imprimir("LBLOQUE[j][j]",LBloques[j][j]);
            //trsm L21=A21 L11^-t       y   
            //trsm L31 =A31 L11^-t      y  
            for(int i=j+1;i<k_global;i++){
                System.out.printf("hacia abajo ↓ %d\t",i); 
                LBloques[i][j] = TRSM(ABloques[i][j],LBloques[j][j]); //
            }
            //syrk A22=A22 -  L21L21^t,    paralelizar
            //syrk A32=A32 - L31 L21^t
            //syrk A33=A33 -  L31L31^t
            for(int i=j+1;i<k_global;i++){
                System.out.printf("en las diagonales y ↓ %d\t",i);
                for(int ii=i;ii<k_global;ii++){
                    ABloques[ii][i] = SYRK(ABloques[ii][i],LBloques[ii][j],LBloques[i][j]);
                }
            }
        }    
    }
    //------------------------------------------------------------------------------
    public static void CopiarBloque(double[][] destino, double[][] bloque, int filaInicio, int columnaInicio) {
        int rows = bloque.length;
        int cols = bloque[0].length;

        for (int i = 0; i < rows; i++) {
            System.arraycopy(bloque[i], 0, destino[filaInicio + i], columnaInicio, cols);
        }
    }
    //------------------------------------------------------------------------------
    public static double[][] TRSM(double[][]Ab,double[][]L){
        double[][]U=Transponer(L);
        Imprimir("dentro de TRSM transpuesta L", L);
        Imprimir("L^t", U);
        double[][]INV =InversaTriangularSuperior(U);
        Imprimir("L^-t",INV);
        Imprimir("L^t L^-t", ProdParalelo(U, INV));
        //double [][] SOL =new double[M.length][INVP[0].length];
        //double[][] PROD = ProdSerial(Ab, INV);
        //Imprimir("actualizacion", PROD);
        return ProdParalelo(Ab, INV); 
    }
    //---------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    public static double[][]SYRK(double[][]M,double[][]L,double[][]Lt){
        //syrk A22=A22 -  L21L21^t,    paralelizar
        //syrk A32=A32 - L31 L21^t
        //syrk A33=A33 -  L31L31^t
        double[][]LT = Transponer(Lt);
        double[][] LLT = ProdParalelo(L, LT); 
        return Resta(M,LLT);
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
        
        for (int i = 0; i < k_global; i++) {
            for (int j = 0; j < k_global; j++) {
                ABloques[i][j] = new double[NB][NB];  // Asignando submatrices de tamaño NB x NB
            }
        } 
        //dividir en k bloques de tamaño NB(modo SERIAL)
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
        Thread [] columnas =new Thread[k_global];
        for (int hil=0;hil<k_global;hil++) {       //un hilo por columna se obtienen k(k+1)/2 bloques ,pues A es simetrica
            final int  hilo = hil;
            columnas[hilo]=new Thread(new Runnable(){
                public void run(){
                    for(int abajo=hilo;abajo<k_global;abajo++){
                        for(int i=abajo*NB;i<(abajo+1)*NB;i++){
                            System.arraycopy(A_Global[i], hilo*NB, ABloques[abajo][hilo][i-abajo*NB], 0, NB);
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
    //-------------------------------------------------------------------------------------------
    public static double[][] CholeskySerial(double[][] M){
        int n = M.length;
        double [][]G =new double[M.length][M.length];
        double[][]MM=Copiar(M);
        //Imprimir(G);
        for(int j =0;j<n;j++){
            double suma = 0;
            for(int k =0;k<j;k++){
                suma+=(G[j][k]*G[j][k]);
            }
            //System.out.printf("\nM[j][j]-suma:%1.2f\n",M[j][j]-suma);
            G[j][j] = Math.pow(MM[j][j]-suma,0.5);
            for(int i =j+1;i<n;i++){
                double sumai = 0;
                for(int m =0;m<j;m++){
                    sumai+=(G[i][m]*G[j][m]);       
                }
                G[i][j] = (MM[i][j]-sumai)/G[j][j];
            }
        }
        return G;
    }   
    //-------------------------------------------------------------------------------------------
    public static double[][] ProdSerial(double[][]M1,double[][]M2){
        double[][]PROD = new double[M1.length][M2[0].length];
        for(int i=0;i<M1.length;i++){
            for(int j=0;j<M1[0].length;j++){
                double suma = 0;
                for(int k=0;k<M2.length;k++){
                    suma +=M1[i][k]*M2[k][j];
                }
                PROD[i][j] = suma;
            }
        }
        return PROD;
    }
 //----------------------------------------------------------------------------
 public static double[][] ProdParalelo(double[][]M1,double[][]M2){
    double [][]PROD=new double[M1.length][M2[0].length];
    int N=M1.length;
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
    //------------------------------------------------------------------------------
    public static double[][] InversaTriangularSuperior(double[][]U){
        double [][] INV = new double[U.length][U.length];
        double factor = Math.pow(10,15);
        for(int j=U.length-1;j>=0;j--){
            if(U[j][j]==0){
                throw new ArithmeticException("no existe inversa");
            }
            INV[j][j] =1/U[j][j];
            for(int i=j-1;i>=0;i--){
                double suma = 0;
                for(int k=i+1;k<=j;k++){
                    suma +=(U[i][k]*INV[k][j]);    
                }
                INV[i][j] = -Math.round(suma*factor/U[i][i])/factor;
            }
        }
        return INV;
    }
    //------------------------------------------------------------------------------
    public static double[][] InversaTriangularInferior(double[][]L){
        double factor = Math.pow(10,15);
        double[][]INV = new double[L.length][L.length];
        for(int j=0;j<L.length;j++){
            if(L[j][j]==0){
                throw new ArithmeticException("La matriz no es invertible");
            }
            INV[j][j]=Math.round(factor/L[j][j])/factor;
            //System.out.println(" "+INV[j][j]);
            for(int i=j+1;i<L.length;i++){
                double suma = 0;
                for(int k=j;k<=i-1;k++){
                    suma += (L[i][k] * INV[k][j]);  
                }
                INV[i][j] = -Math.round(suma*factor/L[i][i])/factor; 
            }
        }
        return INV;
    }
    //------------------------------------------------------------------------
    public static void ImprimirInversa(double[][]M){
        for(int i=0;i<M.length;i++){
            for(int j=0;j<M.length;j++){
                System.out.printf("%1.15f ",M[i][j]);
            }
            System.out.println();
        }
    }
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
		for(int j=0;j<n;j++){		        //operaciones por filas gauss
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
        Imprimir("M",M);
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
            for (int i=0;i<N_global*N_global;i++){
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
                A_Global[i][j]=A_Global[j][i]= num;      //pues sera simetrica
                if(j == i){           // j hasta j==i (la diagonal)
                    j = 0;
                    i++;
                    if(i ==N_global){
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
                A_Global[i][j]= num;      //pues sera simetrica
                if(j == N_global-1){           // j hasta j==i (la diagonal)
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
       A_Global=ProdParalelo(A_Global, Transponer(A_Global));
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
    public static synchronized void Imprimir(String cad,double[][] M){
        System.out.println(cad);
        for(int i=0;i<M.length;i++){
            for(int j=0;j<M[0].length;j++){
                System.out.printf("%1.5f ",M[i][j]);
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
            for(int i=0;i<N_global;i++){
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
        double [][]X =new double[1][N_global];
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
                double p = ProductoTriple(x, A_Global);
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
//-----------------------------------------------------------------------------------------------------------------------------------
class ClassPositiva{
    double[][]M;
    boolean flag;
    ClassPositiva(double[][]m,boolean f){
        this.M=m;
        this.flag=f;
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------
class DataSetCholesky{
    private static int filas;
    private static int columnas;
    private static String name;
    private static int BLOCK;
    private static String mode;
    private static LinkedList<Thread> hilos = new LinkedList<>();
    DataSetCholesky(String nam,int fil,int col,int block,String mod){
        this.filas = fil;
        this.columnas = col;
        this.name =nam;
        this.BLOCK=block;
        this.mode =mod;
        if(mod.equals("PrintWriter")){
            CreateDataPrintWriter();
        }else{
            CreateDataFileWriter();
        }
        
    }
    //---------------------------------------------------------------------------
    public static void CreateDataPrintWriter(){
        double num;
        Random ran = new Random();
        PrintWriter pw = null;
        try{
            pw = new PrintWriter(name);
            for(int i = 0;i<filas;i++){
                for(int j=0;j<columnas;j++){
                    double min = Math.pow(10,BLOCK-2);
                    double max= Math.pow(10,BLOCK-1);
                    System.out.printf("max= %f , min =%f\n",max,min);
                    num = min + ran.nextDouble()*(max - min);
                    pw.printf("%1.1f",num);
                }
                pw.println("");
            }
            pw.close();
        }catch(IOException e){
            System.out.println(e.getMessage());
        }
    }
    //---------------------------------------------------------------------------------
    public static void CreateDataFileWriter(){
        double num;
        Random ran = new Random();
        FileWriter pw = null;
        try{
            pw = new FileWriter(name);
            for(int i = 0;i<filas;i++){
                for(int j=0;j<columnas;j++){
                    double min = Math.pow(10,BLOCK-2);
                    double max= Math.pow(10,BLOCK-1);
                    //System.out.printf("max= %f , min =%f\n",max,min);
                    num = min + ran.nextDouble()*(max - min);
                    pw.write((long)num+" ");
                }
            }
            pw.close();
        }catch(IOException e){
            System.out.println(e.getMessage());
        }
    }
    //-------------------------------------------------------------------------------
    public static double[][] ReadDataScanner(){
        if(mode.equals("PrintWriter")){
            double [][]M=new double[filas][columnas];
        Scanner scanner = null;
        double num;
        Path ruta = Paths.get(name);
        try{
            scanner = new Scanner(Files.newInputStream(ruta));
            scanner.useDelimiter("[;\\s]+");
            int i =0;
            int j =0;
            while(scanner.hasNext()){
                if(scanner.hasNextDouble()){
                    M[i][j] = scanner.nextDouble();
                    System.out.printf("%1.1f\t",M[i][j]);
                    j++;
                    if(j == columnas){
                        j = 0;
                        i++;
                        System.out.println();
                    }
                }else{
                    scanner.next();
                }
            }
            scanner.close();
        }catch(IOException e){
            System.out.println(e.getMessage());
        }
        return M;
        }else return null;
    }
    //---------------------------------------------------------------------------------------------------------------
    public static void WriteData(String FILE){
        if(mode.equals("PrintWriter")){
            PrintWriter pw = null;
        double[][] A= ReadDataScanner();
        try{    
            pw = new PrintWriter(FILE);
            for(int i =0;i<A.length;i++){
                for(int j=0;j<A[0].length;j++){
                    //System.out.printf("write%1.1f\t",A[i][j]);
                    pw.printf("%1.1f;",A[i][j]);
                }
                System.out.println();
                pw.println("");
            }
            pw.close();
        }catch(FileNotFoundException e){
            System.out.println(e.getMessage());
        }
        }
    }
    //----------------------------------------------------------------------------------------------------------
    public static void WriteData(String FILE,double[][]M){
        if(mode.equals("PrintWriter")){
            PrintWriter pw = null;
        try{
            pw = new PrintWriter(FILE);
            for(int i=0;i<M.length;i++){
                for(int j=0;j<M[0].length;j++){
                    //System.out.printf("write%1.1f\t",M[i][j]);
                    pw.printf("%1.1f;",M[i][j]);
                }
                //System.out.println();
                pw.println("");
            }
            pw.close();
        }catch(FileNotFoundException e){
            System.out.println(e.getMessage());
        }
        finally{
            pw.close();
        }
        }
    }
    //--------------------------------------------------------------------------------------------------------
    public static String BufferToString(byte[]BUFFER){
        String cad = "";
        for(int i=0;i<BUFFER.length;i++){
            cad = cad + (char)BUFFER[i];
        }
        return cad;
    }
    //-----------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------
    public static double[][] ReadDataRAF(String FILE){
        if(mode.equals("FileWriter")){
            double [][]M=new double[filas][columnas];
        RandomAccessFile RAF =null;
        try{
            RAF = new RandomAccessFile(FILE,"r");
            int W = (int)RAF.length()/(filas*columnas);
           // System.out.println("dentro de readRAF"+W);
            byte []BUFFER = new byte[W];
            for(int i =0;i<filas;i++){
                for(int j=0;j<columnas;j++){
                    RAF.seek(i*columnas*W + W*j);       //buena forma de direccionar
                    RAF.read(BUFFER);
                    String cad = BufferToString(BUFFER).trim();
                    M[i][j] = Double.parseDouble(cad);
                    System.out.printf("%1.1f\t", M[i][j]);
                    if(j==columnas - 1){
                        System.out.println();
                    }
                }
            }
            RAF.close();
        }catch(IOException e){
            System.out.println(e.getMessage());
        }
        return M;
        }
        else{return null;}
    }
    //---------------------------------------------------------------------------------------------------
    public static double[][] ReadDataRAF(){
        if(mode.equals("FileWriter")){
            double [][]M=new double[filas][columnas];
        RandomAccessFile RAF =null;
        try{
            RAF = new RandomAccessFile(name,"r");
            int W = (int)RAF.length()/(filas*columnas);
           // System.out.println("dentro de readRAF"+W);
            byte []BUFFER = new byte[W];
            for(int i =0;i<filas;i++){
                for(int j=0;j<columnas;j++){
                    RAF.seek(i*columnas*W + W*j);       //buena forma de direccionar
                    RAF.read(BUFFER);
                    String cad = BufferToString(BUFFER).trim();
                    M[i][j] = Double.parseDouble(cad);
                    System.out.printf("%1.1f\t", M[i][j]);
                    if(j==columnas - 1){
                        System.out.println();
                    }
                }
            }
            RAF.close();
        }catch(IOException e){
            System.out.println(e.getMessage());
        }
        return M = ProdParalelo(M, Transponer(M));
        }
        else{return null;}
    }
    //---------------------------------------------------------------------------------------------------------
    public static void WriteDataRAF(String FILE){
        if(mode.equals("FileWriter")){
            double [][]M=new double[filas][columnas];
        RandomAccessFile RAF =null;
        FileWriter fw = null;
        try{
            RAF = new RandomAccessFile(name,"r");
            fw = new FileWriter(FILE);
            int W = (int)RAF.length()/(filas*columnas);
           // System.out.println("dentro de readRAF"+W);
            byte []BUFFER = new byte[W];
            for(int i =0;i<filas;i++){
                for(int j=0;j<columnas;j++){
                    RAF.seek(i*columnas*W + W*j);
                    RAF.read(BUFFER);
                    String cad = BufferToString(BUFFER).trim();
                    M[i][j] = Double.parseDouble(cad);
                    fw.write((long)M[i][j]+" ");
                }
            }
            RAF.close();
            fw.close();
        }catch(IOException e){
            System.out.println(e.getMessage());
            }
         
        }
        
    }
    //---------------------------------------------------------------------------------------------
    public static double[][] ProdParalelo(double[][]M1,double[][]M2){
        double [][]PROD=new double[M1.length][M2[0].length];
        int N=M1.length;
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
    //--------------------------------------------------------------------------------------------------------
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
    //---------------------------------------------------------------------------------------------------------
}
