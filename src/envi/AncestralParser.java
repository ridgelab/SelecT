package envi;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import errors.FileParsingException;
import tools.SNP;
import tools.Window;
import tools.Log;

public class AncestralParser {
	
	private static int TEST_LINES = 10; //for testing the first 10 lines of a file to ensure it is the right format
	private static String LEGEND_TYPE = "legend";
	private static String EMF_TYPE = "emf";
	private static String EMF_SEQ = "SEQ";
	private static String EMF_DATA = "DATA";
	private static String EMF_H_SAPIEN = "homo_sapien";
	
	private boolean skip_first_line;
	private int chr;
	
	private String anc_path;
	private Scanner anc_scan;
	private Log log;

	public AncestralParser(String anc_path, int chr, Log log) throws FileParsingException {
		
		this.chr = chr;
		this.anc_path = anc_path;
		this.log = log;
		
		skip_first_line = false;
		
		try {
			anc_scan = new Scanner(new File(anc_path));
		} catch (FileNotFoundException e) {
			String msg = "Error: Legend File not found in the /Ancestral directory";
			throw new FileParsingException(log, msg);
		}
	}
	
	/**
	 * For use on ancestral legend files; only the first four columns are parsed
	 * Ancestral type is defined as a0 allele
	 * 
	 * @return		Returns a list of SNP objects to be referenced later on when calculating population genetic statistics
	 */
	public List<Window> parseAncestralTypes() throws FileParsingException {//out file here
		log.addLine("Importing Ancestral data from " + anc_path);
		
		List<Window> anc_types = new ArrayList<Window>();
		
		if(anc_path.contains(LEGEND_TYPE))
			anc_types = parseLegendFile();
		else if(anc_path.contains(EMF_TYPE))
			anc_types = parseEmfFile();
		
		return anc_types;
	}
	
	private List<Window> parseEmfFile() throws FileParsingException {
		
		List<Window> anc_types = new ArrayList<Window>();
		
		int num_seq = 0;
		int st_pos = -1;
		int end_pos = -1;
		
		while(anc_scan.hasNextLine()) {
			
			String line = anc_scan.nextLine();
			String[] ln_arr = line.split("\\s+");
			
			if(line.length() >= 1 && ln_arr[0].charAt(0) != '#') {
				
				if(ln_arr[0].equals(EMF_SEQ)) {
					num_seq++;
					
					if(ln_arr[1].contains(EMF_H_SAPIEN) && num_seq == 1) {
						st_pos = Integer.parseInt(ln_arr[3]);
						end_pos = Integer.parseInt(ln_arr[4]);
					}
				}
				if(ln_arr[0].equals(EMF_DATA))  {
					
					int num_anc_seq = (num_seq - 1) / 2;
					
					anc_types = collectData(anc_types, anc_scan, st_pos, end_pos, num_anc_seq);
					
					num_seq = 0;
				}
					
			}
		}
		
		//************for printing to a file****************
		System.out.println("Starting Print\t" + anc_types.size());
//		try {
//			PrintWriter pw = new PrintWriter(out_file);
//		
//			for(int i = 0; i < anc_types.size(); i++) {
//				pw.write(anc_types.get(i) + "\n");
//			}
//			
//			pw.close();
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		}
//		System.out.println("Ending Print");
		//****************************************************
		
		return anc_types;
	}
	
	private List<Window> collectData(List<Window> anc_types, 
								Scanner anc_scan, 
								int st_pos, 
								int end_pos,
								int num_anc_seq) {
		
		int line_num = 0;
		String line = anc_scan.nextLine();
		while(!line.equals("//")) {
			if(line.charAt(0) != '-') {
				
				int cntA = 0;
				int cntT = 0;
				int cntG = 0;
				int cntC = 0;
				
				for(int i = 0; i < num_anc_seq; i++) {
					
					int seq_indx = 1 + i*2;
					if(seq_indx > line.length())
						break;
					else if(line.charAt(seq_indx) == 'a' || line.charAt(seq_indx) == 'A')
						cntA++;
					else if(line.charAt(seq_indx) == 't' || line.charAt(seq_indx) == 'T')
						cntT++;
					else if(line.charAt(seq_indx) == 'g' || line.charAt(seq_indx) == 'G')
						cntG++;
					else if(line.charAt(seq_indx) == 'c' || line.charAt(seq_indx) == 'C')
						cntC++;
				}
				
				if(line.charAt(1) != '-') {
					
					String a0 = getAncestralAllele(cntA, cntT, cntG, cntC, line);
					int pos = line_num + st_pos;
					String snp_id = chr + ":" + pos;
					
					if(hasValidWindow(anc_types, pos)) {
						
						for(int i = 0; i < anc_types.size(); i++) {
							Window w = anc_types.get(i);
							if(w.getStPos() < pos && w.getEndPos() >= pos) {
								w.addSNP(new SNP(pos, a0, "-", snp_id));
								anc_types.set(i, w);
							}
						}
					}
					else {
						
						Window w = new Window(st_pos, end_pos);
						w.addSNP(new SNP(pos, a0, "-", snp_id));
						anc_types.add(w);
					}
				}
				
				//increase only when you have a human allele
				line_num++;
			}
			line = anc_scan.nextLine();
		}
		
		return anc_types;
	}
	
	private String getAncestralAllele(int cntA, int cntT, int cntG, int cntC, String line) {
		
		if(cntA > cntT && cntA > cntG && cntA > cntC)
			return "A";
		if(cntT > cntA && cntT > cntG && cntT > cntC)
			return "T";
		if(cntG > cntT && cntG > cntA && cntG > cntC)
			return "G";
		if(cntC > cntT && cntC > cntG && cntC > cntA)
			return "C";
		
		return Character.toString(line.charAt(1));//the default is return the first ancestral value (if there is a tie)
	}
	
	//makes 1Mb windows
	private List<Window> parseLegendFile() throws FileParsingException {
		
		List<Window> anc_types = new ArrayList<Window>();
		
		checkLegendFile();
		
		if(skip_first_line)
			anc_scan.nextLine(); //skips the first line header
		
		final int WIN_SIZE = 1000000;
		int st_pos = 0;
		int end_pos = WIN_SIZE;
		while(anc_scan.hasNextLine()) {
			
			String[] ln_arr = anc_scan.nextLine().split("\\s+");
			
			int pos = Integer.parseInt(ln_arr[1]);
			
			if(hasValidWindow(anc_types, pos)) {
				
				for(int i = 0; i < anc_types.size(); i++) {
					Window w = anc_types.get(i);
					if(w.getStPos() < pos && w.getEndPos() >= pos) {
						w.addSNP(new SNP(pos, ln_arr[2], ln_arr[3], ln_arr[0]));
						anc_types.set(i, w);
					}
				}
			}
			else {
				
				while(st_pos < pos && end_pos < pos) {
					st_pos += WIN_SIZE;
					end_pos += WIN_SIZE;
				}
				
				Window w = new Window(st_pos, end_pos);
				w.addSNP(new SNP(pos, ln_arr[2], ln_arr[3], ln_arr[0]));
				anc_types.add(w);
				
				st_pos += WIN_SIZE;
				end_pos += WIN_SIZE;
			}
		}
		
//		for(Window w : anc_types)
//			System.out.println(w);
		
		return anc_types;
		
	}
	
	private boolean hasValidWindow(List<Window> wins, int pos) {
		
		for(Window w : wins) {
			if(w.getStPos() < pos && w.getEndPos() >= pos) 
				return true;
		}
		
		return false;
	}
	
	private void checkLegendFile() throws FileParsingException {
		
		Scanner temp_scan = null;
		try {
			temp_scan = new Scanner(new File(anc_path));
			
			//to check the first line; to skip or not to skip, that is the question
			String first_line = temp_scan.nextLine();
			String[] first_line_arr = first_line.split("\\s+");
			if(first_line_arr[1].contains("pos"))
				skip_first_line = true;
			if(first_line_arr[0].contains("rs") && first_line_arr[0].length() < 3)
				skip_first_line = true;
			
			for(int i = 0; i < TEST_LINES; i++) {
				
				String line = temp_scan.nextLine();
				String[] line_arr = line.split("\\s+");
				
				if(line_arr.length < 4 || line_arr.length > 7) {
					String msg = "Error: Legend file " + anc_path + " has invalid number of columns";
					throw new FileParsingException(log, msg);
				}
				
				Integer.parseInt(line_arr[1]);//check for position in position column
			}
			
		} catch (FileNotFoundException e) {
			String msg = "Error: Legend File not found in the /Ancestral directory";
			throw new FileParsingException(log, msg);
		} catch (NumberFormatException e) {
			String msg = "Error: Legend File " + anc_path + " has incorrect colum formatting";
			throw new FileParsingException(log, msg);
		}
	}
}
