import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
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
    List<Marker> _defInRegion; 
    
    public RecursiveGeneAnalyzer(double indexPvalueThreshold, double suggestivePvalueThreshold, String inputFileLocation, int searchSpace, String outputFileLocation) {
        _indexPvalueThreshold = indexPvalueThreshold;
        _suggestivePvalueThreshold = suggestivePvalueThreshold;
        _inputFileLocation = Paths.get(inputFileLocation);
        _searchSpace = searchSpace;
        _outputFileLocation = outputFileLocation;
        
        _resultSet = new ArrayList<Region>();
        _defInRegion = new ArrayList<Marker>();
    }
    
    public void RunThroughDataset() {
    	
    	UploadDataset();
    	final String[] chromosomes = new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"};
    	
    	for(int x = 0; x < chromosomes.length; x++) {
    		String currentChromosomeUnderAnalysis = chromosomes[x]; 
    		List<Marker> workingChromosomeSet = _totalDataSet.stream()
    											.filter(y -> y.getChromosome().equals(currentChromosomeUnderAnalysis))
    											.collect(Collectors.toList());
    		
    		if(!ShouldRunFurtherAnalysis(workingChromosomeSet)) {
    			System.out.println("Could not find any markers under the p-value threshold from chromosome " + currentChromosomeUnderAnalysis);
    		}
    		else {
    			List<Marker> IndexLevelMarkers = GetPositionsOfIndexMarkers(workingChromosomeSet);

    			for(Marker candidate : IndexLevelMarkers){
    				System.out.println(candidate.Chromosome);
    				List<Marker> entireRegionSet = TraverseSearch(workingChromosomeSet, candidate.Position);
    				if(!_defInRegion.containsAll(entireRegionSet)) {
    					_defInRegion.addAll(entireRegionSet);
    					
    					
    					ConstructRegionAndAppendToResultSet(entireRegionSet);
    				}
    				
    				System.out.println("---- ");
    			}
    			
    		}
    	}
    	
    	try {
			BuildResultFile();
		} catch (IOException e) {
			e.printStackTrace();
		}
    	
    }

	private void ConstructRegionAndAppendToResultSet(List<Marker> entireRegionSet) {
		Collections.sort(entireRegionSet, new Comparator<Marker>(){
	   	     public int compare(Marker o1, Marker o2){
	   	         if(o1.Position == o2.Position)
	   	             return 0;
	   	         return o1.Position < o2.Position ? -1 : 1;
	   	     }
	   	});
    	var sug = entireRegionSet.stream()
		.filter(p -> p.getPvalue() <= _suggestivePvalueThreshold)
		.collect(Collectors.toList());
		
    	sug.forEach(x -> System.out.println("Position: " + x.Position + " Sugg P Val: " + x.Pvalue));
    	
		Region newRegion = new Region(){};
		newRegion.setRegionStart(entireRegionSet.get(0).Position);
       	
		Marker smallestMarker = GetMarkerWithLowestPValue(entireRegionSet); 
		newRegion.setChr(Integer.parseInt(smallestMarker.getChromosome()));
    	newRegion.setMarkerName(smallestMarker.getName());
    	newRegion.setPosition(smallestMarker.getPosition());
    	newRegion.setPvalue(smallestMarker.getPvalue());
    	
    	newRegion.setNumSigMarkers(entireRegionSet.stream()
					.filter(p -> p.getPvalue() <= _indexPvalueThreshold)
					.collect(Collectors.toList()).size());
    	

		newRegion.setNumSuggestiveMarkers(entireRegionSet.stream()
							.filter(p -> p.getPvalue() <= _suggestivePvalueThreshold)
							.collect(Collectors.toList()).size());
		
		newRegion.setNumTotalMarkers(entireRegionSet.size());
		
		Collections.sort(entireRegionSet, new Comparator<Marker>(){
	   	     public int compare(Marker o1, Marker o2){
	   	         if(o1.Position == o2.Position)
	   	             return 0;
	   	         return o1.Position > o2.Position ? -1 : 1;
	   	     }
	   	});
		
		newRegion.setRegionStop(entireRegionSet.get(0).Position);
		newRegion.setSizeOfRegion();
		newRegion.setRegionIndex(_resultSet.size() + 1);
		
		_resultSet.add(newRegion);
	}

	private Marker GetMarkerWithLowestPValue(List<Marker> entireRegionSet) {
		Collections.sort(entireRegionSet, new Comparator<Marker>(){
   	     public int compare(Marker o1, Marker o2){
   	         if(o1.Pvalue == o2.Pvalue)
   	             return 0;
   	         return o1.Pvalue < o2.Pvalue ? -1 : 1;
   	     }
   	});
		return entireRegionSet.get(0);
	}

	private List<Marker> GetPositionsOfIndexMarkers(List<Marker> workingChromosomeSet) {
		 return workingChromosomeSet.stream()
				.filter(x -> x.getPvalue() <= _indexPvalueThreshold)
				.collect(Collectors.toList());
	}

	private List<Marker> TraverseSearch(List<Marker> workingChromosomeSet, int indexPosition) {
		//While there are results within +/- the searchSpace under the suggestive value threshold, keep expanding your search
		
		System.out.println("Traverse search: ");
		// Minus
		int startSearchStart = SeekStartPosition(workingChromosomeSet, indexPosition); 
		// Plus
		int stopSearchStop = SeekStopPosition(workingChromosomeSet, indexPosition); 
		System.out.println("Start: " + startSearchStart); 
		System.out.println("Stop: " + stopSearchStop); 
		
		return workingChromosomeSet.stream()
				.filter(x -> x.getPosition() >= startSearchStart && x.getPosition() <= stopSearchStop)
				.collect(Collectors.toList());
	}
	
	private int SeekStopPosition(List<Marker> workingChromosomeSet, int indexPosition) {
		
		int stopSearchStop = indexPosition; 
		int previousSearchStop = 0;
		
		
		boolean keepGoingDown = true; 
		while(keepGoingDown) {
			var stopping = ExpandPlusDirection(workingChromosomeSet, stopSearchStop); 
			
			var suggestive = stopping.stream()
					.filter(y -> y.getPvalue() <= _suggestivePvalueThreshold)
					.collect(Collectors.toList());
			
			
			Collections.sort(suggestive , new Comparator<Marker>(){
		   	     public int compare(Marker o1, Marker o2){
		   	    	return o1.Position > o2.Position ? -1 :(o1.Position < o2.Position ? 1 : 0);
		   	     }
		       	});
			
			
			stopSearchStop = suggestive.get(0).Position;

			if(suggestive.size() > 0 && previousSearchStop != stopSearchStop) {
				stopSearchStop = suggestive.get(0).Position;
			}
			else {
				keepGoingDown = false; 
			}
			
			previousSearchStop = stopSearchStop; 
		}
		
		return stopSearchStop; 	
	}
	
	private int SeekStartPosition(List<Marker> workingChromosomeSet, int indexPosition) {
		
		int startSearchStart = indexPosition; 
		int previousSearchStart = 0;
		
		boolean keepGoingUp = true; 
		while(keepGoingUp) {
			var starting =  ExpandMinusDirection(workingChromosomeSet, startSearchStart);
			
			var suggestive = starting.stream()
					.filter(y -> y.getPvalue() <= _suggestivePvalueThreshold)
					.collect(Collectors.toList());
		
			if(suggestive.size() > 0 && previousSearchStart != startSearchStart) {
				startSearchStart = suggestive.get(0).Position;
				
			}
			else {
				keepGoingUp = false;
			}
			
			previousSearchStart = startSearchStart; 
		}
		
		return startSearchStart; 
	}
	
	private List<Marker>ExpandPlusDirection(List<Marker> workingChromosomeSet, int position){
		return workingChromosomeSet
		.stream()
		.filter(y -> y.getPosition() >= position && 
				y.getPosition() <= position + _searchSpace )
		.collect(Collectors.toList());
	}
	
	private List<Marker>ExpandMinusDirection(List<Marker> workingChromosomeSet, int position){
		return workingChromosomeSet
		.stream()
		.filter(y -> y.getPosition() <= position && 
				y.getPosition() >= position - _searchSpace )
		.collect(Collectors.toList());
	}

	private boolean ShouldRunFurtherAnalysis(List<Marker> workingChromosomeSet) {
		return workingChromosomeSet.stream()
											.filter(x -> x.getPvalue() <= _indexPvalueThreshold)
											.collect(Collectors.toList()).size() > 0; 
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
	
	private void BuildResultFile()throws IOException {
    	File fout = new File(_outputFileLocation);
    	FileOutputStream fos = new FileOutputStream(fout);
     
    	BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fos));
     
    	var header = new String[] { "Region", "MarkerName", "Chr", "Position", "P-value", "RegionStart",
                "RegionStop", "NumSigMarkers", "NumSuggestiveMarkers",  "NumTotalMarkers", "SizeOfRegion" };
    	
    	bw.write(String.join("\t", header));
    	bw.newLine();
    	
        for(Region result : _resultSet)
        {
            var arrayResult = new String[]
            {
                String.valueOf(result.getRegionIndex()),
                result.getMarkerName(),
                String.valueOf(result.getChr()),
                String.valueOf(result.getPosition()),
                String.valueOf(result.getPvalue()),
                String.valueOf(result.getRegionStart()),
                String.valueOf(result.getRegionStop()),
                String.valueOf(result.getNumSigMarkers()),
                String.valueOf(result.getNumSuggestiveMarkers()),
                String.valueOf(result.getNumTotalMarkers()),
                String.valueOf(result.getSizeOfRegion()),
            };
            
            bw.write(String.join("\t", arrayResult));
            bw.newLine();
        }
     
    	bw.close();
    }
}
