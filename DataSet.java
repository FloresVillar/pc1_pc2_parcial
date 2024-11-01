/* Realizar un programa Java serial y paralelo 
(usando hilos explícitos creados e instanciados en el 
programa) para solucionar el problema de la 
descomposición matricial QR de una matriz cuadrada 
NxN (). Verificar que los resultados sean equivalentes
 y comparar los tiempos de ejecución.
PROGRAMA QUE GENERA LOS DATOS*/
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Random;
import java.util.Scanner;
import java.io.FileNotFoundException;
public class DataSet {
static int filas = 1000;
static int columnas = 1000;
static String DATAFILE = "DATAPC1Preg1PotenciaMatriz.TXT";
//------------------------------------------------------
public static void CreateFile(int M, int N) {
    Random random = new Random();
    PrintWriter pw = null;
    try {   pw = new PrintWriter(DATAFILE);
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < N; j++) {
                    double valor = random.nextDouble() * 100;
                    pw.printf("%10.0f", valor);
                }
            pw.println("");
            }
    } catch (FileNotFoundException ex) {
        System.out.println(ex.getMessage());
    } finally {
        pw.close();
    }
}
//------------------------------------------------------------------
public static double[][] ReadFile(int M, int N) {
    double[][] MTX = new double[M][N];
    Path filePath = Paths.get(DATAFILE);//Paths es una clase muy simple con un unico meodo static get()
    Scanner scanner = null;               //devuelve un objeto de tipo cadena o URI
        try {
                scanner = new Scanner(filePath);
        }
        catch (IOException ex) {
                System.out.println(ex.getMessage());
        }
        int i = 0, j = 0;       
        while (scanner.hasNext()) {
                if (scanner.hasNextDouble()) {
                        MTX[i][j] = scanner.nextDouble();
                        j++;
                        if (j == N) {
                                j = 0;
                                i++;
                        }
                }
                else {
                        scanner.next();
                }
        }
        return MTX;
}
//------------------------------------------------------
public static void WriteFile(Matrix M, String archivo){
    int mm = M.getRows();
    int nn = M.getCols();
        PrintWriter pw = null;
        try {
                pw = new PrintWriter(archivo);
                for (int i = 0; i < mm; i++) {
                        for (int j = 0; j < nn; j++) {
                                double valor = M.GetCell(i, j);
                                pw.printf("%10.0f ", valor);
                        }
                        pw.println("");
                }
        }
        catch (FileNotFoundException ex) {
                System.out.println(ex.getMessage());
        }
        finally {
                pw.close();
        }
}
//------------------------------------------------------
public static void main(String[] args) {
   Scanner scanner = new Scanner(System.in);
   CreateFile(filas, columnas);
   scanner.close(); 
   //double [][]A=ReadFile();   
}
}
//fin clase DataSet	
/*Al procesar los programas anteriores 
(Serial y Paralelo)  se crearon los siguientes 
archivos
Comparando los archivos generados 
en modo serial y paralelo
Se verifica que son correctos*/ 
class Matrix {
   private int ROWS;
   private int COLS;
   private double[][] M;

   public Matrix(double[][] var1) {
      this.M = var1;
      this.ROWS = var1.length;
      this.COLS = var1[0].length;
   }

   public double GetCell(int var1, int var2) {
      return this.M[var1][var2];
   }

   public void SetCell(int var1, int var2, double var3) {
      this.M[var1][var2] = var3;
   }

   public int getRows() {
      return this.ROWS;
   }

   public int getCols() {
      return this.COLS;
   }

   public void imprimir() {
      for(int var1 = 0; var1 < this.ROWS; ++var1) {
         for(int var2 = 0; var2 < this.COLS; ++var2) {
            System.out.printf("%12.2f", this.M[var1][var2]);
         }

         System.out.println();
      }

   }

   public synchronized void incrementar(int var1, int var2, double var3) {
      double[] var10000 = this.M[var1];
      var10000[var2] += var3;
   }

   public double prodEsc(int var1, int var2, int var3, int var4) {
      double var5 = 0.0;

      for(int var7 = var3; var7 <= var4; ++var7) {
         var5 += this.M[var7][var1] * this.M[var7][var2];
      }

      return var5;
   }

   public double prodEsc(int var1, int var2) {
      double var3 = 0.0;

      for(int var5 = 0; var5 < this.ROWS; ++var5) {
         var3 += this.M[var5][var1] * this.M[var5][var2];
      }

      return var3;
   }

   Matrix prod(Matrix var1) {
      int var2 = var1.getCols();
      Matrix var3 = new Matrix(new double[this.ROWS][var2]);

      for(int var4 = 0; var4 < this.ROWS; ++var4) {
         for(int var5 = 0; var5 < var2; ++var5) {
            for(int var6 = 0; var6 < this.COLS; ++var6) {
               double var7 = this.GetCell(var4, var6) * var1.GetCell(var6, var5);
               var3.SetCell(var4, var5, var3.GetCell(var4, var5) + var7);
            }
         }
      }

      return var3;
   }
}
