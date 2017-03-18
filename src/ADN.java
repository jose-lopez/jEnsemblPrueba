import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import uk.ac.roslin.ensembl.config.DBConnection;
import uk.ac.roslin.ensembl.config.EnsemblCoordSystemType;
import uk.ac.roslin.ensembl.config.FeatureType;
import uk.ac.roslin.ensembl.dao.database.DBCollectionSpecies;
import uk.ac.roslin.ensembl.dao.database.DBRegistry;
import uk.ac.roslin.ensembl.dao.database.DBSpecies;
import uk.ac.roslin.ensembl.datasourceaware.DAXRef;
import uk.ac.roslin.ensembl.datasourceaware.compara.DAHomologyPairRelationship;
import uk.ac.roslin.ensembl.datasourceaware.core.DAChromosome;
import uk.ac.roslin.ensembl.datasourceaware.core.DACoordinateSystem;
import uk.ac.roslin.ensembl.datasourceaware.core.DAExon;
import uk.ac.roslin.ensembl.datasourceaware.core.DAGene;
import uk.ac.roslin.ensembl.datasourceaware.core.DATranscript;
import uk.ac.roslin.ensembl.exception.ConfigurationException;
import uk.ac.roslin.ensembl.exception.DAOException;
import uk.ac.roslin.ensembl.exception.NonUniqueException;
import uk.ac.roslin.ensembl.model.Mapping;
import uk.ac.roslin.ensembl.datasourceaware.core.DADNASequence;
import uk.ac.roslin.ensembl.model.core.Chromosome;

public class ADN {
    public String nombre_especie;
    public String DBversion;
    public String gen_id;
    DBRegistry eReg;
    public DBSpecies especie;
    public DAGene Gen;
    public List<DAHomologyPairRelationship> homologos;
    public List<Homologo> region_promotora;
    Mapping mapping;
    public ArrayList secuencias;
    
    
    public ADN(String nombre_especie, String gen_id, String DBversion) throws ConfigurationException, DAOException, NonUniqueException {
        this.nombre_especie = nombre_especie;
        this.DBversion = DBversion;
        try{           
            eReg = new DBRegistry(DBConnection.DataSource.ENSEMBLDB);
            
            
            try{
            especie = eReg.getSpeciesByAlias(this.nombre_especie);
            Gen = especie.getGeneByStableID(gen_id, DBversion);
            mapping = Gen.getChromosomeMapping();
            homologos = Gen.getHomologies();
            region_promotora = new ArrayList<Homologo>();
            secuencias = new ArrayList();
            }catch(Exception e){
                System.out.println("No se encontro la especie");
            };
        }catch(Exception e){
            
            System.out.println("Error al conectar la base de datos");
        };
    }
    
    
    public void getInformacionEspecie() {
        System.out.println("-----------------Información de la especie -------------------------");
        System.out.println("Nombre: " + Gen.getSpecies().getCommonName());
        System.out.println("Taxionomia: " + especie.getTaxonomyID());
    }
    
    
    public void getInformacionGen() throws DAOException, NonUniqueException {
        System.out.println("-----------------Información del gen "+Gen+"-------------------------");
        System.out.println("Gene: " + Gen.getStableID());
        System.out.println(homologos.size() + " homólogos encontrados");
        System.out.println("\tchr start: " + mapping.getTargetCoordinates().getStart());
        System.out.println("\tchr stop: " + mapping.getTargetCoordinates().getEnd());
        System.out.println("\tassembly: " + Gen.getAssembly());
        System.out.println("\tdescription: " + Gen.getDescription());
        System.out.println("\tsymbol: " + Gen.getDisplayName());
        System.out.println("\tstrand: " + mapping.getTargetCoordinates().getStrand());
        System.out.println("\ttaxonID: " + Gen.getSpecies().getTaxonomyID());
        System.out.println("\tstatus: " + Gen.getStatus());
        System.out.println("\ttype: " + Gen.getBiotype());
            
        //System.out.println("Secuencia: " + Gen.getSequenceAsString().length());
    }
    
    String _id;
    String _CadGenerada="";
    public void getInformacionHomologo() throws DAOException, NonUniqueException, IOException {
        System.out.println("-----------------------Homólogos------------------------------");
        System.out.println(homologos.size() + " homólogos encontrados");
        
        BufferedReader bufferRead = new BufferedReader(new InputStreamReader(System.in));
        String ans;
        
        int i = 1;
        boolean exit = false;
        System.out.print("¿Buscar region promotora de todos los homologos? y/n");
        ans = bufferRead.readLine();
        if(ans.equalsIgnoreCase("y")){
            System.out.print("Ingrese el valor inferior: ");
            String va = bufferRead.readLine();
            System.out.print("Ingrese el valor superior: ");
            String vb = bufferRead.readLine();
            int valorA = Integer.parseInt(va);
            int valorB = Integer.parseInt(vb);
            
            for (DAHomologyPairRelationship h : homologos) {
                System.out.print(i + " ");
                System.out.print(h.getTargetProperties().getSpeciesName());
                System.out.print(" Gen: " + h.getTarget().getStableID());
                _id=">"+h.getTarget().getStableID();
                System.out.print(" [" + h.getType().toString() + "] (Ancestro común: " + h.getLastCommonAncestor() + ")");
                System.out.println("["+h.getTargetProperties().getCoords().toString()+"]");                
                regionPromotora(valorA, valorB, i - 1);                
                i++;                
            }       
            generarArchivo(secuencias, "secuencias");
            GenerarFasta(_CadGenerada);
        }        
        else
        {
            for (DAHomologyPairRelationship h : homologos) {
                System.out.print(i + " ");
                System.out.print(h.getTargetProperties().getSpeciesName());
                System.out.print(" Gen: " + h.getTarget().getStableID());
                _id=">"+h.getTarget().getStableID();
                System.out.print(" [" + h.getType().toString() + "] (Ancestro común: " + h.getLastCommonAncestor() + ")");
                System.out.println("["+h.getTargetProperties().getCoords().toString()+"]");

                do {
                    System.out.print("¿Desea buscar region promotora para este homologo? (y/n/f)");
                    ans = bufferRead.readLine();
                    if (ans.equalsIgnoreCase("y")) {
                        System.out.print("Ingrese el valor inferior: ");
                        String va = bufferRead.readLine();
                        System.out.print("Ingrese el valor superior: ");
                        String vb = bufferRead.readLine();

                        int valorA = Integer.parseInt(va);
                        int valorB = Integer.parseInt(vb);
                        regionPromotora(valorA, valorB, i - 1);
                    }
                    else if (ans.equalsIgnoreCase("n")) {
                        //no se ingresa el homologo a las regiones promotoras escogidas 
                    }
                    else if (ans.equalsIgnoreCase("f")) {
                        exit = true;
                        break;
                    }
                } while(ans.equalsIgnoreCase("y") && ans.equalsIgnoreCase("n") && ans.equalsIgnoreCase("f"));
                i++;
                if (exit){
                    generarArchivo(secuencias, "secuencias");
                    GenerarFasta(_CadGenerada);
                    break;
                }
            }
        }
    }
    
    
    public void regionPromotora(int valorA, int valorB, int index) {
        String[] splitted = homologos.get(index).getTargetProperties().getCoords().toString().split(" ");
        int inicio = Integer.parseInt(splitted[0]);
        int inicio_promotor = inicio - valorA;
        if (inicio_promotor < 1) {
            System.out.println("Region Promotora: ");
            return ;
        }
        int fin_promotor = inicio - valorB;
        
        DAGene gen_cr = null;
        Chromosome cromosoma_homo = null;
        try {
            gen_cr = especie.getGeneByStableID(homologos.get(index).getTarget().getStableID(), DBversion);
            try {
                cromosoma_homo = (Chromosome) gen_cr.getChromosomeMapping().getTarget();
            } catch (NonUniqueException ex) {
                //Logger.getLogger(ADN.class.getName()).log(Level.SEVERE, null, ex);
            } catch (NullPointerException nex) {
                
            }
        } catch (DAOException ex) {
            //Logger.getLogger(ADN.class.getName()).log(Level.SEVERE, null, ex);
        }
        if (cromosoma_homo != null) {
            Homologo h = new Homologo(homologos.get(index), cromosoma_homo.getSequenceAsString(inicio_promotor, fin_promotor));
            region_promotora.add(h);
            secuencias.add(region_promotora.get(region_promotora.size() - 1).region_promotora);
            System.out.println("Region Promotora : " + region_promotora.get(region_promotora.size() - 1).region_promotora); //ultima posicion
            _CadGenerada+=_id+"\n"+region_promotora.get(region_promotora.size() - 1).region_promotora+"\n";                        
        }
        else {
            Homologo h = new Homologo(homologos.get(index), "null");
            region_promotora.add(h);
            System.out.println("Region Promotora: ");
        }
        System.out.println("");
    }
    
    
    public void getInformacionExones() {
        System.out.println("-----------------Información de los Exones -------------------------");
        System.out.println("\tTranscript Count: " + Gen.getTranscripts().size());
        for (DATranscript t : Gen.getTranscripts()) {
            System.out.println("\t\tTranscript: " + t.getStableID());
            System.out.println("\t\t " + t.getDisplayName());
            System.out.println("\t\t " + t.getStatus());
            System.out.println("\t\t " + t.getBiotype());
            System.out.println("\t\t " + t.getDescription());
            System.out.println("\t\t " + t.getGene().getStableID());
            System.out.println("\t\tCanonical ?  " + t.isCanonical());
            System.out.println("\t\tXREF: " + t.getDisplayXRef().getDBDisplayName());
            System.out.println("\t\tXREF: " + t.getDisplayXRef().getDisplayID());
            System.out.println("\t\tXREF: " + t.getDisplayXRef().getInfoType());
            System.out.println("\t\tXREF: " + t.getDisplayXRef().getInfo());
            for (Mapping m : t.getLoadedMappings(EnsemblCoordSystemType.chromosome)) {
                System.out.println("\t\tMapping: " + m.getTargetHashID());
                System.out.println("\t\t\tCoords: " + m.getTargetCoordinates().toString());
            }

            // look at all the exons of the transcript 
            System.out.println("EXONS");
            for (DAExon e : t.getExons()) {
                System.out.println("\t\tRank: " + e.getRank());
                System.out.println("\t\tStableID: " + e.getStableID());
                System.out.println("\t\tID: " + e.getId());
                System.out.println("\t\tstart phase: " + e.getPhase());
                System.out.println("\t\tend phase: " + e.getEndPhase());
                System.out.println("\t\tcurrent: " + e.isCurrent());
                System.out.println("\t\tconstitutive: " + e.isConstitutive());

                //get locations of exon
                for (Mapping m : e.getLoadedMappings(EnsemblCoordSystemType.chromosome)) {
                    System.out.println("\t\tMapping: " + m.getTargetHashID());
                    System.out.println("\t\t\tCoords: " + m.getTargetCoordinates().toString());
                }
            }
        }
    }
    public void GenerarFasta(String cadenas){
        FileWriter fichero = null;
        PrintWriter pw = null;
        try
        {
            fichero = new FileWriter("secuencia.fasta");
            pw = new PrintWriter(fichero);            
            pw.println(cadenas.toString());
 
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
           try {
           if (null != fichero)
              fichero.close();
           } catch (Exception e2) {
              e2.printStackTrace();
           }
        }
    }
    public void generarArchivo(ArrayList cadenas,String nomarch){
        FileWriter fichero = null;
        PrintWriter pw = null;
        try
        {
            fichero = new FileWriter(nomarch+".stam");
            pw = new PrintWriter(fichero);
 
            for (int i = 0; i < cadenas.size(); i++)
                pw.println(cadenas.get(i).toString());
 
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
           try {
           if (null != fichero)
              fichero.close();
           } catch (Exception e2) {
              e2.printStackTrace();
           }
        }
System.out.println("Se ha generado el archivo "+nomarch+".stam con las secuencias de ADN solicitadas");        
    }    
}
