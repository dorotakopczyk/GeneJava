public class Region {

    public int RegionIndex;
    public String MarkerName;
    public int Chr;
    public double Pvalue;
    public int Position;
    public int RegionStart;
    public int RegionStop;
    public int NumSigMarkers;
    public int NumSuggestiveMarkers;
    public int NumTotalMarkers;

    
    public int SizeOfRegion;

    public void setSizeOfRegion() {
    	this.SizeOfRegion = this.RegionStop - this.RegionStart; 
    }
    
    public int getSizeOfRegion() {
    	return this.SizeOfRegion;
    }
    
    public void setRegionIndex(int name) {
        this.RegionIndex = name;
    }

    public int getRegionIndex() {
        return this.RegionIndex;
    }


    public void setMarkerName(String name) {
        this.MarkerName = name;
    }

    public String getMarkerName() {
        return this.MarkerName;
    }
    
    public void setChr(int name) {
        this.Chr = name;
    }

    public int getChr() {
        return this.Chr;
    }

    public void setPvalue(double name) {
        this.Pvalue = name;
    }

    public double getPvalue() {
        return this.Pvalue;
    }
    
    public void setPosition(int name) {
        this.Position = name;
    }

    public int getPosition() {
        return this.Position;
    }
    
    public void setRegionStart(int name) {
        this.RegionStart = name;
    }

    public int getRegionStart() {
        return this.RegionStart;
    }

    public int getRegionStop() {
        return this.RegionStop;
    }
    
    public void setRegionStop(int name) {
        this.RegionStop = name;
    }

    public int getNumSigMarkers() {
        return this.NumSigMarkers;
    }
    
    public void setNumSigMarkers(int name) {
        this.NumSigMarkers = name;
    }

    public int getNumSuggestiveMarkers() {
        return this.NumSuggestiveMarkers;
    }
    
    public void setNumSuggestiveMarkers(int name) {
        this.NumSuggestiveMarkers = name;
    }

    public int getNumTotalMarkers() {
        return this.NumTotalMarkers;
    }
    
    public void setNumTotalMarkers(int name) {
        this.NumTotalMarkers = name;
    }
}