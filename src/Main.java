import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.Authenticator;
import java.net.Socket;


public class Main {
     
    public static boolean ComprobarConexion(){
        String dirWeb = "www.ebi.ac.uk";       
        int puerto = 80;

                try{
                        Socket s = new Socket(dirWeb, puerto);
                        if(s.isConnected()){                            
                            return true;
                        }
                        return false;
        }catch(Exception e){
                return false;
        }        
    }    
    public static void AlinearSecuencias() throws IOException{
        System.out.print("¿Desea alinear las secuencias? (y/n)");
        BufferedReader bufferRead = new BufferedReader(new InputStreamReader(System.in));
        String resp = bufferRead.readLine();
        if(resp.equals("y"))
        {
            try {
                String cmd = "clustalw2 secuencia.fasta"; 
                Runtime.getRuntime().exec(cmd);
                System.out.print("Secuencia alineada con exito");
                } catch (IOException ioe) {
                    System.out.println (ioe);}
        }
    }
    public static void  Iniciar() throws IOException{
            String nombre_especie;
            String gen_id;
            String DBversion;
            ADN objAdn;                                                                 
            
            System.out.println("-------------Genes-------------");

            System.out.print("Ingrese el nombre de la especie: ");
            BufferedReader bufferRead = new BufferedReader(new InputStreamReader(System.in));
            nombre_especie = bufferRead.readLine();

            System.out.print("Ingrese el ID del gen: ");
            gen_id = bufferRead.readLine();

            System.out.print("Ingrese la version de la base de datos: ");
            DBversion = bufferRead.readLine();

            
            try{
                //objAdn = new ADN("cow", "ENSBTAG00000021527", "68");
                objAdn = new ADN(nombre_especie, gen_id, DBversion);            
                objAdn.getInformacionEspecie();        
                System.out.println("");
                //System.out.println("---------------------------------------");
                System.out.println("");
                objAdn.getInformacionGen();

                System.out.println("");
                //System.out.println("---------------------------------------");
                System.out.println("");
                objAdn.getInformacionHomologo();
                
                System.out.println("---------------------------------------");
                AlinearSecuencias();
            }catch(Exception e){
                System.out.println("No se encontro informacion");
            }
    }    
    public static void main(String[] args) throws Exception {
                
        if(ComprobarConexion())
        {            
            Iniciar();            
        }
        else{
            String ans,proxy,port,username,password;
            BufferedReader bufferRead = new BufferedReader(new InputStreamReader(System.in));
            System.out.println("Esta navegando detras de un proxy? (y/n)");
            ans = bufferRead.readLine();
            
            if(ans.equalsIgnoreCase("y"))                
            {
                System.out.println("Direccion del proxy: ");
                proxy = bufferRead.readLine();
                System.out.println("puerto: ");
                port = bufferRead.readLine();
                System.out.println("Usuario: ");
                username = bufferRead.readLine();
                System.out.println("Pass: ");
                password = bufferRead.readLine();                                               
                
                Authenticator.setDefault(new ProxyAuthenticator(username, password));                
                System.setProperty("http.proxyHost", proxy);
                System.setProperty("http.proxyPort", port);                
                if(ComprobarConexion())
                   Iniciar();
                else
                    System.out.println("Error al conectar el proxy");
            }else{
                System.out.println("Verifique su conexion a internet");
            }
           
        }
    }        
}
// Para ejecutar ingresar
// cow
// ENSBTAG00000021527
// 68