
public class UniqueRegions {

	// java UniqueRegions input.txt output.txt 0.00001 0.0001 500000
	public static void main(String[] args) {
		var inputFileLocation = args[0];
        var outputFileLocation = args[1];
        var indexPvalueThreshold = args[2];
        var suggestivePvalueThreshold = args[3];
        var searchSpace = args[4];
        
        var geneAnalyzer = new RecursiveGeneAnalyzer(Double.parseDouble(indexPvalueThreshold), 
        		Double.parseDouble(suggestivePvalueThreshold), 
        		inputFileLocation, 
        		Integer.parseInt(searchSpace), 
        		outputFileLocation);
        
        //geneAnalyzer.GetMyRegions();
        geneAnalyzer.RunThroughDataset();
        
        System.out.println("Completed. Please find output at " + outputFileLocation);
	}

}
