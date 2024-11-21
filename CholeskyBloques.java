import java.io.File;
import java.io.FileWriter;
import java.io.RandomAccessFile;
import java.io.IOException;
import java.util.random.*;
//==============================================================================================
public class CholeskyBloques {
    private static int BLOCK = 4; //tamaño de datos BLOCK -1=3, el ultimo 'digito' es " "
    private static byte [] RECORD = new byte[BLOCK];
    private static String FILENAME = "DATACholeskyBloques.TXT";
    private static int N = 8;
    private static String CADENA;
    private static double [][] A =new double[N][N];
    //------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------
    
    //------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------
    public static void main(String[]args){
        //mientras no sea simetrica definida positiva
        boolean flag = false;
        while(!flag){
            Generar();
            //flag = esPositiva(A);
            flag= true;
        }
    }
    //-------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------
    public static double[][] CholeskySerial(double[][] M){
        int n = M.length;
        double [][]G =new double[M.length][M.length];
        for(int j =0;j<n;j++){
            double suma = 0;
            for(int k =0;k<j;k++){
                suma+=G[j][k]*G[j][k];
            }
            G[j][j] = Math.pow(M[j][j]-suma,0.5);
            for(int i =j+1;i<n;i++){
                double sumai = 0;
                for(int k =0;k<j;k++){
                    sumai+=G[i][k]*G[j][k];       
                }
                G[i][j] = (A[i][j]-sumai)/G[j][j];
            }
        }
        return G;
    }   
    //-------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------
    public static void WriteData(String FILENAME){
        double num;
        try{
            FileWriter fw = new FileWriter(FILENAME);
            //generando el aletario
            double max= Math.pow(10,BLOCK-1);
            double min = Math.pow(10,BLOCK-2);
            int l = (int)(N*(N+1)/2);   //para una matriz simetrica
            for (int i=0;i<l;i++){
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
    public static void Generar(){
        WriteData(FILENAME);
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
                A[i][j]=A[j][i]= num;
                if(j == i){           // j hasta j==i (la diagonal)
                    j = 0;
                    i++;
                    if(i ==N){
                        break;  //fila 4 , no posible
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
        Imprimir(A);
        
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
    public static void Imprimir(double[][] M){
        System.out.println();
        for(int i=0;i<M.length;i++){
            for(int j=0;j<M[0].length;j++){
                System.out.printf("%1.2f\t",M[i][j]);
            }
            System.out.println();
        }
        System.out.println();
    }
    //------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------
    public static boolean esPositiva(double [][]M){
        
        return true;
    }
    
}
