import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class RecursiveGeneAnalyzer {
	double _indexPvalueThreshold;
    double _suggestivePvalueThreshold;
    Path _inputFileLocation;
    int _searchSpace;
    String _outputFileLocation;

    List<Region> _resultSet;
    List<Marker> _totalDataSet;
    
    public RecursiveGeneAnalyzer(double indexPvalueThreshold, double suggestivePvalueThreshold, String inputFileLocation, int searchSpace, String outputFileLocation) {
        _indexPvalueThreshold = indexPvalueThreshold;
        _suggestivePvalueThreshold = suggestivePvalueThreshold;
        _inputFileLocation = Paths.get(inputFileLocation);
        _searchSpace = searchSpace;
        _outputFileLocation = outputFileLocation;
        
        _resultSet = new ArrayList<Region>();
    }
    
    public void RunThroughDataset() {
    	
    	UploadDataset();
    	final String[] chromosomes = new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"};
    	
    	for(int x = 0; x < chromosomes.length; x++) {
    		String currentChromosomeUnderAnalysis = chromosomes[x]; 
    		
    		List<Marker> workingChromosomeSet = _totalDataSet.stream()
    											.filter(y -> y.Chromosome == currentChromosomeUnderAnalysis)
    											.collect(Collectors.toList());
    		if(!ShouldRunFurtherAnalysis(workingChromosomeSet)) {
    			System.out.println("Could not find any markers under the p-value threshold from chromosome " + currentChromosomeUnderAnalysis);
    		}
    	}
    	
    }
    

	private boolean ShouldRunFurtherAnalysis(List<Marker> workingChromosomeSet) {
		boolean hasIndexValueThresholdResult = workingChromosomeSet.stream()
											.filter(x -> x.getPvalue() <= _indexPvalueThreshold)
											.collect(Collectors.toList()).size() >= 1; 
											
		return hasIndexValueThresholdResult;
	}

	public void UploadDataset() {
    	
    	List<String> dataset = null;
    	try (Stream<String> lines = Files.lines(_inputFileLocation)) 
    	{
    		System.out.println("Found file and loading...");
    		dataset  = lines.collect(Collectors.toList());
    	} 
    	catch (Exception e) {
    		System.out.println(e);
    	}
    	
    	dataset.remove(0); 
    	
    	_totalDataSet = TransformInputFileToListOfObjects(dataset);
    }

	private List<Marker> TransformInputFileToListOfObjects(List<String> dataset) {
		List<Marker> filedata = new ArrayList<Marker>();
		for(String line : dataset) {
			var wordsArray = line.split("\t");
			
			var record = new Marker();
            record.Name = wordsArray[0];
            record.Chromosome = wordsArray[1];
            record.Position = Integer.parseInt(wordsArray[2]);
            
            try {
                /* We suspect that this block of statement can throw 
                 * exception when p-value is NA
                 */
                Integer.parseInt(record.Chromosome);
            	record.setPvalue(Double.parseDouble(wordsArray[wordsArray.length-1]));
                filedata.add(record);

             }
             catch (NumberFormatException e) { 
                System.out.println("Could not parse pvalue for " + wordsArray[3]);
             }
		}
		
    	var cleaned = filedata.stream().filter(x -> !x.Chromosome.isEmpty()).toArray(Marker[]::new);
    	var noBlanks = Arrays.asList(cleaned);
    	
		return noBlanks; 
	}
}
