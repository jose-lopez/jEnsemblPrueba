
import uk.ac.roslin.ensembl.datasourceaware.compara.DAHomologyPairRelationship;

public class Homologo {
    public DAHomologyPairRelationship homologo;
    public String region_promotora;
    
    public Homologo(DAHomologyPairRelationship homologo, String region_promotora)
    {
        this.homologo = homologo;
        this.region_promotora = region_promotora;
    }
}
