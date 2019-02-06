import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;
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


public class GeneAnalyzer {

    double _indexPvalueThreshold;
    double _suggestivePvalueThreshold;
    Path _inputFileLocation;
    int _searchSpace;
    String _outputFileLocation;

    List<Region> _resultSet;

    public GeneAnalyzer(double indexPvalueThreshold, double suggestivePvalueThreshold, String inputFileLocation, int searchSpace, String outputFileLocation) {
        _indexPvalueThreshold = indexPvalueThreshold;
        _suggestivePvalueThreshold = suggestivePvalueThreshold;
        _inputFileLocation = Paths.get(inputFileLocation);
        _searchSpace = searchSpace;
        _outputFileLocation = outputFileLocation;
        
        _resultSet = new ArrayList<Region>();
    }

    public List<Region> GetMyRegions() {
    	
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
 
    	
    	var markers = TransformInputFileToListOfObjects(dataset);
    	
    	//Split everything up by chromosome
    	List<List<Marker>> chromosomeSets = markers.stream()
                .collect(Collectors.groupingBy(Marker :: getChromosome)) 
                .entrySet()
                .stream()
                .map(Map.Entry::getValue)
                .collect(Collectors.toList());
    	
    	Comparator<List<Marker>> comparator = (list1, list2) -> {
    	       
    		    var chrom1 = Integer.valueOf(list1.get(0).Chromosome);
    		    var chrom2 = Integer.valueOf(list2.get(0).Chromosome);
    		    return chrom1.compareTo(chrom2);     
    	};
    	
    	chromosomeSets.sort(comparator);
    	
    	// All the work is done for a set of chromosomes,
        // where first the chromosomes are in order and then so are positions. (genomic order)
    	for(List<Marker> chromosomeSet : chromosomeSets) {
    		
    		// ChromosomeSet cant be assigned to bc its a foreach iteration variable
            // thus it's immutable. It needs to be assignable. This is important in case we need to 
            // manipulate our list 
            var workingChromosome = chromosomeSet;
            
            //Verify position
            var isSorted = IsSorted(workingChromosome);
            if (!isSorted)
            {
            	System.out.println("Markers for chromosome" + workingChromosome.get(0).Chromosome + " are not in order. Reshuffling...");
            	throw new ArithmeticException("Not in genomic order");
            }
            
            var stepOneCandidates = GetRecordsExceedingIndexThreshold(workingChromosome);
            
            if (stepOneCandidates.size() == 0)
            {
            	System.out.println("For chromosome" + workingChromosome.get(0).Chromosome + 
            			", no markers found with an index p-value threshold exceeding " + _indexPvalueThreshold );
            }
            else
            {
            	// We will then search 500,000 base pairs in both directions 
                // We can now begin defining a region. Expand the search +/- 500k (position) 
            	for(Marker candidate : stepOneCandidates) {
            		var stepTwoCandidates = GetExpandedSearchSpace(workingChromosome, candidate.getPosition());
            	
            		// We can now build our regions and add them to the result set.
                    var regionCandidates = GetRecordsExceedingSuggestiveThreshold(stepTwoCandidates);
                    if(regionCandidates.size() > 0) {
                    	for(Marker regionCandidate : regionCandidates) {
                    		
                    		// Then we will extend the window for another 500,000 base pairs beyond that and continue searching.
                            // First, define the region by expanding the search results +/- 500k again
                            var expandedResults = GetExpandedSearchSpace(workingChromosome, regionCandidate.getPosition());
                            
                            // We will define the start and stop positions of the region as the positions of the first and last marker
                            // in the region that meet the SUGGESTIVE THRESHOLD.
                            var newRegion = BuildRegion(_resultSet, expandedResults, regionCandidate);
                
                            if (IsNewMarker(newRegion.MarkerName) && !OverlapsWithPreviousRegion(newRegion))
                            {
                                newRegion.setRegionIndex(_resultSet.size() + 1);
                                _resultSet.add(newRegion);
                            }
                            else {
                            	FixUpRegion(newRegion, expandedResults);
                            }
                    	}
                    }
            	}
            }
            
    	}
    	
    	   // However, there are several regions with many markers below this threshold that are not unique. These markers are correlated due to the underlying
        // structure of the genome (linkage disequilibrium), and we want to identify and summarize all of the UNIQUE regions.
        try {
			BuildResultFile(_outputFileLocation);
		} catch (IOException e) {
			e.printStackTrace();
		}
        
        return _resultSet;
        
    }
    
    private void FixUpRegion(Region newRegion, List<Marker> chromosomeSet) {
    	 List<Region> regionsNeedingFixup = _resultSet.stream().filter(x -> (newRegion.RegionStart - x.RegionStart) < _searchSpace &&
                 (newRegion.RegionStart - x.RegionStart) >= 0 &&
                 x.Chr == newRegion.Chr).collect(Collectors.toList());

			if (regionsNeedingFixup.size() == 0)
			{
			return;
			}
			
			Region regionNeedingFixup = regionsNeedingFixup.get(0);
			
			if (regionNeedingFixup.RegionStart > newRegion.RegionStart)
			{
			regionNeedingFixup.RegionStart = newRegion.RegionStart;
			}
			if (regionNeedingFixup.RegionStop < newRegion.RegionStop)
			{
			regionNeedingFixup.RegionStop = newRegion.RegionStop;
			}
			
			var region = chromosomeSet.stream().filter(x -> x.Pvalue < _suggestivePvalueThreshold).collect(Collectors.toList());
			
			regionNeedingFixup.NumSigMarkers = region.stream()
						.filter(p -> p.getPvalue() <= _indexPvalueThreshold)
						.collect(Collectors.toList()).size();
			
			regionNeedingFixup.NumSuggestiveMarkers = region.stream()
						.filter(p -> p.getPvalue() <= _suggestivePvalueThreshold)
						.collect(Collectors.toList()).size();
			
			regionNeedingFixup.NumTotalMarkers = chromosomeSet.stream().filter(x -> x.Position >= regionNeedingFixup.RegionStart 
					&& x.Position <= regionNeedingFixup.RegionStop).collect(Collectors.toList()).size();
			regionNeedingFixup.SizeOfRegion = regionNeedingFixup.RegionStop - regionNeedingFixup.RegionStart + 1;
		
	}

	private boolean IsNewMarker(String newRegionMarkerName)
    {
        if (_resultSet.stream().filter(x -> x.MarkerName == newRegionMarkerName).collect(Collectors.toList()).size() > 0)
        {
            return false;
        }

        return true;
    }

    
    private boolean OverlapsWithPreviousRegion(Region newRegion)
    {
        if (_resultSet.size() < 1)
        {
            return false;
        }
        
        var index = _resultSet.size();
        var previousRegion = _resultSet.get(index - 1);
        if ((newRegion.RegionStart - previousRegion.getRegionStart() < 500000) 
        		&& (newRegion.RegionStart - previousRegion.RegionStart > 0))
        {
            return true;
        }
        return false; 
    }
    
    private void BuildResultFile(String fileLocation)throws IOException {
    	File fout = new File(fileLocation);
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



	private Region BuildRegion( List<Region> resultSet, List<Marker> chromosomeSet, Marker regionCandidate) {
		
    	// We will define the start and stop positions of the region as the positions of the first and last marker
        // in the region that meet the SUGGESTIVE THRESHOLD. 
    	var region = chromosomeSet.stream().filter(x -> x.getPvalue() < _suggestivePvalueThreshold)
    			.collect(Collectors.toList());
    	
    	Collections.sort(region, new Comparator<Marker>(){
    	     public int compare(Marker o1, Marker o2){
    	         if(o1.Pvalue == o2.Pvalue)
    	             return 0;
    	         return o1.Pvalue < o2.Pvalue ? -1 : 1;
    	     }
    	});
    	
    	var newRegion = new Region();
    	newRegion.setChr(Integer.parseInt(region.get(0).getChromosome()));
    	newRegion.setMarkerName(region.get(0).getName());
    	newRegion.setPosition(region.get(0).getPosition());
    	newRegion.setPvalue(region.get(0).getPvalue());
    	
       	Collections.sort(region, new Comparator<Marker>(){
   	     public int compare(Marker o1, Marker o2){
   	         if(o1.Position == o2.Position)
   	             return 0;
   	         return o1.Position < o2.Position ? -1 : 1;
   	     }
       	});
       	
       	
       	newRegion.setRegionStart(region.get(0).Position);
       	newRegion.setRegionStop(region.get(region.size() - 1).Position);
       	newRegion.setNumSigMarkers(region.stream()
       								.filter(p -> p.getPvalue() <= _indexPvalueThreshold)
       								.collect(Collectors.toList()).size());
       	newRegion.setNumSuggestiveMarkers(region.stream()
       								.filter(p -> p.getPvalue() <= _suggestivePvalueThreshold)
       								.collect(Collectors.toList()).size());
       	
    	newRegion.setNumTotalMarkers(region.size());

    	newRegion.setSizeOfRegion();
    	
		return newRegion;
	}

	private List<Marker> GetRecordsExceedingSuggestiveThreshold(List<Marker> workingSet) {
    	List<Marker> results = new ArrayList<Marker>(); 
    	
    	int len=workingSet.size();
    	for(int i=0; i<len; i++) {
    	    if (workingSet.get(i).getPvalue() < _suggestivePvalueThreshold) {
    	        results.add(workingSet.get(i));
    	    }
    	}
    	
        
        return results;
	}

	private List<Marker> GetExpandedSearchSpace(List<Marker> workingChromosome, int position) {
    	var startingPosition = position - _searchSpace;
        var endingPosition = position + _searchSpace;
        var stepTwoCandidates = workingChromosome.stream()
        	    .filter(p -> p.getPosition() >= startingPosition && 
        	    		p.getPosition() <= endingPosition)
        	    .collect(Collectors.toList());  // On that chromosome.

        return stepTwoCandidates;
	}

	public List<Marker> GetRecordsExceedingIndexThreshold( List<Marker> workingSet)
    {
    	// First we search for an index SNP exceeding the index SNP threshold (p<0.00001)
    	List<Marker> results = new ArrayList<Marker>(); 
    	
    	int len=workingSet.size();
    	for(int i=0; i<len; i++) {
    	    if (workingSet.get(i).getPvalue() < _indexPvalueThreshold) {
    	        results.add(workingSet.get(i));
    	    }
    	}
    	
        
        return results;
    }

	public static boolean IsSorted(List<Marker> listMarkers)
    {
    	Marker[] arrMarkers = new Marker[listMarkers.size()];
    	arrMarkers = listMarkers.toArray(arrMarkers);
    	
        for (int i = 1; i < listMarkers.size(); i++)
        {
            var previousPosition = arrMarkers[i - 1].Position;
            var currentPosition = arrMarkers[i].Position;
            if (previousPosition > currentPosition)
            {
                return false;
            }
        }
        return true;
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