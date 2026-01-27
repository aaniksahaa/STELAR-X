package taxon;

public class Taxon {

    public static int count = 0;

    public int id;
    public String label;
    
    
    public Taxon(int i, String lb){
        id = i;
        label = lb;
    }

    public Taxon(String lb){
        id = count++;
        label = lb;
    }

    @Override
    public String toString(){
        return label;
    }

    
} 