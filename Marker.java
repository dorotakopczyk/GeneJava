
public class Marker {
	public String Name; 
    public String Chromosome;
    public int Position;
    public double Pvalue;
    
    public void setName(String name) {
        this.Name = name;
    }

    public String getName() {
        return this.Name;
    }
    
    public void setChromosome(String name) {
        this.Chromosome = name;
    }

    public String getChromosome() {
        return this.Chromosome;
    }
    
    public void setPosition(int name) {
        this.Position = name;
    }

    public int getPosition() {
        return this.Position;
    }
    
    public void setPvalue(double name) {
        this.Pvalue = name;
    }

    public double getPvalue() {
        return this.Pvalue;
    }
}
