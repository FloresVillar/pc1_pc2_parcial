import java.io.*;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;
//Buenas tardes profesor. El proceso Paralelo se ubica dentro del metodo con el mismo nombre ProcesoParalelo()
public class BusquedaSecuencialFloresVillar {
private static String FILENAME = "DATOS.TXT";
private static int    N = 1000000;
private static String KEY = "01234567890";
private static String CADENA;
private static int BLOCK = 11;
private static byte[] RECORD = new byte[BLOCK];
private static LinkedList<Thread> hilos = new LinkedList<Thread>();
//

    //------------------------------------------------
    private static String GetString() {
    String CAD;
        CAD = "";
        for(int i=0;i<=10;i++) {
            CAD = CAD + (char)(RECORD[i]);
        }
        return CAD;
    }
    //------------------------------------------------
    /*private static void PrintRecord() {
        CADENA = "";
        for(int i=0;i<=10;i++) {
            CADENA = CADENA + (char)(RECORD[i]);
        }
        System.out.println(CADENA);
    }*/
    //------------------------------------------------
    private static void ProcesoSerial() {
    long n,P,T,Time1,Time2;
    int i;
        try {
             RandomAccessFile RAF = new RandomAccessFile(FILENAME,"r");
             Time1 = System.currentTimeMillis();
             T = RAF.length();
             n = T/BLOCK;
             P = -1;
             i = 0;
             while((i<=n-1)&&(P==-1)) { //El método java.io.RandomAccessFile.seek(long pos) establece el desplazamiento del puntero de archivo
                 RAF.seek(i*BLOCK); //, medido desde el comienzo de este archivo, en el que se produce la siguiente lectura o escritura.
                 RAF.read(RECORD);  //Lee hasta b.lengthbytes de datos de este archivo en una matriz de bytes.
                 CADENA = GetString();
                 if(CADENA.equals(KEY)==true) {
                    P = i;
                 }
                 i++;
             }
             RAF.close();
             if(P>=0) {
                System.out.println("Elemento " + KEY + " Existente en la Posicion " + P);
             }
             else {
                System.out.println("Elemento " + KEY + " No Existe");
             }
             Time2 = System.currentTimeMillis();
             System.out.printf("Tiempo de Procesamiento Serial: %d%n", (Time2-Time1));
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
//Profesor ,tambien expuse del proceso paralelo para la potencia de una matriz, ello por puntos el viernes, dia que hubo clases virtuales debido al paro......
    private static void ProcesoParalelo(){
        Thread hil1 = new Thread(new Runnable(){
			public void run(){    
                long n,P,T,Time1,Time2;
                int i;
                try {
                    RandomAccessFile RAF = new RandomAccessFile(FILENAME,"r");
                    Time1 = System.currentTimeMillis();
                    T = RAF.length();
                    n = T/BLOCK;
                    P = -1;
                    i = 0;
                    while((i<n/4)&&(P==-1)) { //El método java.io.RandomAccessFile.seek(long pos) establece el desplazamiento del puntero de archivo
                        RAF.seek(i*BLOCK); //, medido desde el comienzo de este archivo, en el que se produce la siguiente lectura o escritura.
                        RAF.read(RECORD);  //Lee hasta b.lengthbytes de datos de este archivo en una matriz de bytes.
                        CADENA = GetString();
                        if(CADENA.equals(KEY)==true) {
                           P = i;
                        }
                        i++;
                    }
                    RAF.close();
                    if(P>=0) {
                       System.out.println("Elemento " + KEY + " Existente en la Posicion " + P);
                    }
                    else {
                       System.out.println("Elemento " + KEY + " No Existe");
                    }
                    Time2 = System.currentTimeMillis();
                    System.out.printf("Tiempo de Procesamiento Serial: %d%n", (Time2-Time1));
               } catch (IOException ex) {
                   ex.printStackTrace();
               }
				}
			});
			hilos.add(hil1);
			hil1.start();
			Thread hil2 = new Thread(new Runnable(){
			public void run(){
        		long n,P,T,Time1,Time2;
                int i;
                try {
                    RandomAccessFile RAF = new RandomAccessFile(FILENAME,"r");
                    Time1 = System.currentTimeMillis();
                    T = RAF.length();
                    n = T/BLOCK;
                    P = -1;
                    i = (int)n/4;
                    while((i<(int)(2*n)/4)&&(P==-1)) { //El método java.io.RandomAccessFile.seek(long pos) establece el desplazamiento del puntero de archivo
                        RAF.seek(i*BLOCK); //, medido desde el comienzo de este archivo, en el que se produce la siguiente lectura o escritura.
                        RAF.read(RECORD);  //Lee hasta b.lengthbytes de datos de este archivo en una matriz de bytes.
                        CADENA = GetString();
                        if(CADENA.equals(KEY)==true) {
                           P = i;
                        }
                        i++;
                    }
                    RAF.close();
                    if(P>=0) {
                       System.out.println("Elemento " + KEY + " Existente en la Posicion " + P);
                    }
                    else {
                       System.out.println("Elemento " + KEY + " No Existe");
                    }
                    Time2 = System.currentTimeMillis();
                    System.out.printf("Tiempo de Procesamiento Serial: %d%n", (Time2-Time1));
               } catch (IOException ex) {
                   ex.printStackTrace();
               } 
        	}
			});
			hilos.add(hil2);
			hil2.start();
			Thread hil3 = new Thread(new Runnable(){
			public void run(){
                long n,P,T,Time1,Time2;
                int i;
                try {
                    RandomAccessFile RAF = new RandomAccessFile(FILENAME,"r");
                    Time1 = System.currentTimeMillis();
                    T = RAF.length();
                    n = T/BLOCK;
                    P = -1;
                    i = (int)(2*n)/4;
                    while((i<(int)(3*n)/4)&&(P==-1)) { //El método java.io.RandomAccessFile.seek(long pos) establece el desplazamiento del puntero de archivo
                        RAF.seek(i*BLOCK); //, medido desde el comienzo de este archivo, en el que se produce la siguiente lectura o escritura.
                        RAF.read(RECORD);  //Lee hasta b.lengthbytes de datos de este archivo en una matriz de bytes.
                        CADENA = GetString();
                        if(CADENA.equals(KEY)==true) {
                           P = i;
                        }
                        i++;
                    }
                    RAF.close();
                    if(P>=0) {
                       System.out.println("Elemento " + KEY + " Existente en la Posicion " + P);
                    }
                    else {
                       System.out.println("Elemento " + KEY + " No Existe");
                    }
                    Time2 = System.currentTimeMillis();
                    System.out.printf("Tiempo de Procesamiento Serial: %d%n", (Time2-Time1));
               } catch (IOException ex) {
                   ex.printStackTrace();
               }
        		 
			}
			});
			hilos.add(hil3);
			hil3.start();
			Thread hil4 = new Thread(new Runnable(){
			public void run(){
                long n,P,T,Time1,Time2;
                int i;
                try {
                    RandomAccessFile RAF = new RandomAccessFile(FILENAME,"r");
                    Time1 = System.currentTimeMillis();
                    T = RAF.length();
                    n = T/BLOCK;
                    P = -1;
                    i = (int)(3*n)/4;
                    while((i<n)&&(P==-1)) { //El método java.io.RandomAccessFile.seek(long pos) establece el desplazamiento del puntero de archivo
                        RAF.seek(i*BLOCK); //, medido desde el comienzo de este archivo, en el que se produce la siguiente lectura o escritura.
                        RAF.read(RECORD);  //Lee hasta b.lengthbytes de datos de este archivo en una matriz de bytes.
                        CADENA = GetString();
                        if(CADENA.equals(KEY)==true) {
                           P = i;
                        }
                        i++;
                    }
                    RAF.close();
                    if(P>=0) {
                       System.out.println("Elemento " + KEY + " Existente en la Posicion " + P);
                    }
                    else {
                       System.out.println("Elemento " + KEY + " No Existe");
                    }
                    Time2 = System.currentTimeMillis();
                    System.out.printf("Tiempo de Procesamiento Serial: %d%n", (Time2-Time1));
               } catch (IOException ex) {
                   ex.printStackTrace();
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
//Los resultados indican que el ultimo hilo encentra lo buscado(aunque no se imprime en orden hilo1, hilo2, hilo3,hilo4) 
//con un tiempo  que aproximadamente la cuarta parte para cada hilo, lo cual es genial 
/*
 * private static void ProcesoParalelo_Idea_Base() {
        long n,P,T,Time1,Time2;
        int i;
            try {
                 RandomAccessFile RAF = new RandomAccessFile(FILENAME,"r");
                 Time1 = System.currentTimeMillis();
                 T = RAF.length();
                 n = T/BLOCK;
                 P = -1;
                 i = 0; //la idea es paralelizar esos 4 bloques de codigo
                 while((i<=n-1)&&(P==-1)) {  
                     RAF.seek(i*BLOCK);  
                     RAF.read(RECORD);   
                     CADENA = GetString();
                     if(CADENA.equals(KEY)==true) {
                        P = i;
                     }
                     if(i+1>n-1){
                        break;
                     }
                     RAF.seek((i+1)*BLOCK); //, medido desde el comienzo de este archivo, en el que se produce la siguiente lectura o escritura.
                     RAF.read(RECORD);  //Lee hasta b.lengthbytes de datos de este archivo en una matriz de bytes.
                     CADENA = GetString();
                     if(CADENA.equals(KEY)==true) {
                        P = i+1;
                     }
                     if(i+2>n-1){
                        break;
                     }
                     RAF.seek((i+2)*BLOCK); //, medido desde el comienzo de este archivo, en el que se produce la siguiente lectura o escritura.
                     RAF.read(RECORD);  //Lee hasta b.lengthbytes de datos de este archivo en una matriz de bytes.
                     CADENA = GetString();
                     if(CADENA.equals(KEY)==true) {
                        P = i+2;
                     }
                     if(i+3>n-1){
                        break;
                     }
                     RAF.seek((i+3)*BLOCK); //, medido desde el comienzo de este archivo, en el que se produce la siguiente lectura o escritura.
                     RAF.read(RECORD);  //Lee hasta b.lengthbytes de datos de este archivo en una matriz de bytes.
                     CADENA = GetString();
                     if(CADENA.equals(KEY)==true) {
                        P = i+3;
                     }
                     i=i+4;
                 }
                 RAF.close();
                 if(P>=0) {
                    System.out.println("Elemento " + KEY + " Existente en la Posicion " + P);
                 }
                 else {
                    System.out.println("Elemento " + KEY + " No Existe");
                 }
                 Time2 = System.currentTimeMillis();
                 System.out.printf("Tiempo de Procesamiento Serial: %d%n", (Time2-Time1));
            } catch (IOException ex) {
                ex.printStackTrace();
            }
    }
 * 
 */
    

    public static String ReadFile() { //obtiene la info , supongo ,,,,,,investigar despues
    String LINE="",CADENA="";
       try {
         File FILE = new File(FILENAME);
         BufferedReader BR = new BufferedReader(new FileReader(FILE));
         while ((LINE = BR.readLine()) != null) {
           CADENA = CADENA + LINE;
         }
         BR.close();
       }
       catch(IOException e) {
       }
       return CADENA;
    }

    public static void WriteData(int N) {       //crea la data
    double X;
    long num;
        try {
            FileWriter FW = new FileWriter(FILENAME);
            for(int i=1;i<=N-1;i++) {
                X = Math.random()*9000000000.0;
                num = 1000000000 + (long)X;
                FW.write(num + " ");
            }
            FW.write("01234567890");
            FW.close();
        }
        catch (IOException E) {
            System.out.print(E.getMessage());
        }
    }


    //--------------------------------------------------
    public static void main(String[] args) {
      WriteData(N);
      ProcesoSerial();
      System.out.println();
      ProcesoParalelo();
    }

}