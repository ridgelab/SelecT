package envi;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

import errors.FileParsingException;
import tools.Log;
import tools.GeneticMap;

public class MapParser {
	
	private static int TEST_LINES = 10; //for testing the first 10 lines of a file to ensure it is the right format
	
	private boolean skip_first_line;
	
	private String map_path;
	private Scanner map_scan;
	private Log log;
	
	public MapParser(String map_file, Log log) throws FileParsingException {
		
		this.map_path = map_file;
		this.log = log;
		
		skip_first_line = false;
		
		try {
			File file = new File(map_path);
			map_scan = new Scanner(file);
			
		} catch (FileNotFoundException e) {
			String msg = "Error: Genetic Map file not found in the map directory";
			throw new FileParsingException(log, msg);
		}
	}
	
	/**
	 * For use on a tab separated Hapmap II genetic Map. Recombination rate is 
	 * defined as the the cM/Mb between the current line's position the position
	 * in the next line (in physical position, not genetic map position)
	 * 
	 * @return		Returns a GeneticMap object with a key defining a position 
	 * 				range on the chromosome and a value equal to combined rate
	 * 				of recombination (cM/Mb)
	 */
	@SuppressWarnings("unused")
	public GeneticMap parseGeneMap() throws FileParsingException {
		log.addLine("Importing genetic map from " + map_path);
		
		checkMapFile();
		
		GeneticMap gm = new GeneticMap();
		
		String[] cur_ln = new String[0];
		String[] nxt_ln = new String[0];
		
		if(skip_first_line)
			map_scan.nextLine();
		
		nxt_ln = map_scan.nextLine().split("\\s+");
		
		if(!nxt_ln[0].matches("[0-9]+")) {
			String msg = "Error: Genetic Map has invalid number format";
			throw new FileParsingException(log, msg);
		}
		if(nxt_ln.length > 3 && !nxt_ln[3].matches("\\s+")) {
			String msg = "Error: Map file " + map_path +" has invalid formatting";
			throw new FileParsingException(log, msg);
		}
		
		int line = 0;
		int cur_pos = 0;
		int nxt_pos = Integer.parseInt(nxt_ln[0]);
		double rate = Double.parseDouble(nxt_ln[1]);
		do {
			
			line++;
			try {
				
				gm.put(cur_pos, nxt_pos, rate);
				
				cur_ln = nxt_ln;
				cur_pos = nxt_pos;
				
				nxt_ln = map_scan.nextLine().split("\\s");
				nxt_pos = Integer.parseInt(nxt_ln[0]);
				
				rate = Double.parseDouble(nxt_ln[1]);
				
			} catch (NumberFormatException e) {
				log.addLine("\tError with line " + line);
				cur_ln = map_scan.nextLine().split("\\s");
				nxt_ln = map_scan.nextLine().split("\\s");
				line+=2;
			}
			
		} while(map_scan.hasNextLine());
		
//		for(Range r : gm.getRangeSet()) {
//			System.out.println(r + " => Double [" + gm.get(r) + "]");
//			log.addLine(r.toString() + " => Double [" + gm.get(r) + "]");
//		}
		
		return gm;
	}
	
	private void checkMapFile() throws FileParsingException {
		
		Scanner temp_scan = null;
		
		try {
			temp_scan = new Scanner(new File(map_path));
			
			String first_line = temp_scan.nextLine();
			String[] first_line_arr = first_line.split("\\s+");
			
			if(!first_line_arr[0].matches("[0-9]+"))
				skip_first_line = true;
			
			for(int i = 0; i < TEST_LINES; i++) {
				String line = temp_scan.nextLine();
				String [] line_arr = line.split("\\s+");
				
				if(line_arr.length != 3) {
					String msg = "Error: Map file " + map_path 
							+" has invalid formatting";
					throw new FileParsingException(log, msg);
				}
				
				Integer.parseInt(line_arr[0]);
				Double.parseDouble(line_arr[1]);
				Double.parseDouble(line_arr[2]);
			}
			
		} catch (FileNotFoundException e) {
			String msg = "Error: Genetic Map file not found in the /Map directory";
			throw new FileParsingException(log, msg);
		} catch (NumberFormatException e) {
			String msg = "Error: Genetic Map has invalid number format";
			throw new FileParsingException(log, msg);
		}	
	}
}
