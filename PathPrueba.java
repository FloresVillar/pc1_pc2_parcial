import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.Files;
public class PathPrueba {
    public static void main(){
    Path archivo = Paths.get("C:\\Users\\FLORES VILLAR\\Desktop\\file_prueba.txt");
    System.out.println(archivo.getFileName());
    System.out.println(archivo.getParent());
    System.out.println(archivo.getRoot());
    System.out.println(archivo.startsWith("la"));
    System.out.println(archivo.isAbsolute());
    //Path file = Files.createFile(Paths.get("C:\\Users\\FLORES VILLAR\\Desktop\\file_prueba.txt"));
    System.out.println(Files.exists(Paths.get("C:\\Users\\FLORES VILLAR\\Desktop\\file_prueba.txt"), null));
    }
}
