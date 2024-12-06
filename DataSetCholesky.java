import java.util.Scanner;
import java.util.Random;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
//file para la data de cholesky
//---------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------
public class DataSetCholesky{
    private static int filas;
    private static int columnas;
    private static String name;
    DataSetCholesky(String nam,int fil, int col){
        this.name = nam;
        this.filas= fil;
        this.columnas = col;
        CreateData();
    }
    //--------------------------------------------------------------------------------------------------------------------
    public static void CreateData(){
        Random rnd = new Random();
        PrintWriter pw = null;
        double num;
        try{
            pw = new PrintWriter(name);
            for(int i=0;i<filas;i++){
                for(int j=0;j<columnas;j++){
                    num = rnd.nextDouble()*100;
                    pw.printf("%1.1f;",num);  
                    //System.out.printf("%1.1f\t",num);     
                }
                pw.println("");
                //System.out.println();
            }
        }catch(FileNotFoundException e){
            System.out.println(e.getMessage());
        }   
        finally{
            pw.close();
        }
    }
    //---------------------------------------------------------------------------------------------------------------------------------------
    public static double[][] ReadData(String FILENAME){
        System.out.println("dentro de ReadData");
        double[][] M=new double[filas][columnas];
        Scanner scanner = null;
        Path path = Paths.get(FILENAME);
        try{    
            scanner = new Scanner(path);//gpt new Scanner(File.newBufferedReader(path));
            scanner.useDelimiter("[;\\s]+");
        }   catch(IOException e){
            System.out.println(e.getMessage());
        } 
        int i=0,j=0;
        while(scanner.hasNext()){
            //System.out.println("dentro de while ReadData");
            if(scanner.hasNextDouble()){
                //System.out.println("dentro de if");
                M[i][j] = scanner.nextDouble();
                System.out.printf("%1.1f\t",M[i][j]);
                j++;
                if(j==columnas){
                    j = 0;
                    i++;
                    System.out.println();
                }
            }
            else{
                scanner.next(); // aqui salta el puntero a la siguiente linea?; pues " " no es entero?
            } 
        } 
        return M;
    }
    //----------------------------------------------------------------------------------------------------------
     public static void WriteData(String FILE){
        PrintWriter pw = null;
        double [][]M=ReadData(name);
        try{
            pw = new PrintWriter(FILE);
            for(int i=0;i<M.length;i++){
                for(int j =0;j<M[0].length;j++){
                    pw.printf("%1.1f ", M[i][j]);
                }
                pw.println("");
            }
            
        }catch(FileNotFoundException e){
            System.out.println(e.getMessage());
        }
        finally{
            pw.close();
        }
     }
    //------------------------------------------------------------------------------------------------------------
    public static void WriteData(String FILE,double[][]M){
        PrintWriter pw = null;
        try{
            pw = new PrintWriter(FILE);
            for(int i=0;i<M.length;i++){
                for(int j =0;j<M[0].length;j++){
                    pw.printf("%1.1f ", M[i][j]);
                }
                pw.println("");
            }
            
        }catch(FileNotFoundException e){
            System.out.println(e.getMessage());
        }
        finally{
            pw.close();
        }
     }
    //------------------------------------------------------------------------------------------------------------------------
    public static void main(String[] args) {
        DataSetCholesky data = new DataSetCholesky("DATACHOLESKY.TXT", 10, 10);
        DataSetCholesky.ReadData(DataSetCholesky.name);
        DataSetCholesky.WriteData("FILECOPY.TXT");
    }
}