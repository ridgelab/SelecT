package envi;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import errors.FileParsingException;
import tools.Individual;
import tools.Log;
import tools.SNP;
import tools.Window;

/**
 * VcfParser is used to parser input VCF files with population data for SelecT run.
 * Data parsed includes individuals, windows, and ancestral types
 */
public class VcfParser {
	
	//Getting rid of the first 9 columns in the 
	//VCF file to only include individuals
	private final int DEFAULT_COL = 9;
	
	private String file_path;
	private int chr; 
	private List<Window> windows;
	private List<Window> ancestral;
	private Individual[] individuals;
	private Log log;
	
	/**
	 * Simple constructor
	 */
	public VcfParser()
	{
		file_path = null; 
		chr = 0; 
		windows = null; 
		ancestral = null; 
		individuals = null;
		log = null; 
	}
	
	/**
	 * Constructor with inputs
	 * 
	 * @param file_path	path to the VCF file that will be parsed
	 * @param chr		chromosome number for the VCF file
	 * @param log		a simple log file to record how the parsing goes
	 */
	public VcfParser(String file_path, int chr, Log log)
	{
		this.file_path = file_path; 
		this.chr = chr; 
		this.log = log;
		
		individuals = new Individual[0]; // Place holder
		windows = new ArrayList<Window>();
		ancestral = new ArrayList<Window>();
	}
	
	/**
	 * Parses the VCF file and stores the resulting data 
	 * (individuals, ancestral data, and windows) as members of the VcfParser
	 * 
	 * @param win_size	how many bases to include in each window
	 * @param anc_data	whether or not this file contains ancestral data
	 * @throws FileParsingException
	 */
	public void parseVCF(int win_size, boolean anc_data) throws FileParsingException {
		log.addLine("Importing VCF data from " +  file_path);
		
		try (Scanner scan = new Scanner(new File(file_path))) {
			String line = "";
			
			//Setup Individual[] that will be populated
			while (scan.hasNext()) {
				
				line = scan.nextLine();
				if (line.charAt(1) != '#') {
					
					String[] ln_arr = line.split("\\s+");
					
					individuals = new Individual[ln_arr.length - DEFAULT_COL];
					for (int i = 0; i < individuals.length; i++) {
						individuals[i] = new Individual(i, chr);
					}
					
					break;
				}
			}
			
			//Add data to Individual[] and create Windows
		    int start_pos = 0;
		    int end_pos = win_size - 1;
		    int index = 0;
		    int pos = 0;
		    Window cur_win = new Window();
		    Window anc_win = new Window();
		    
		    while (scan.hasNext()) {
		    	
		    	line = scan.nextLine();
				String[] ln = line.split("\\s+");
				
				pos = Integer.parseInt(ln[1]);
				
				if (!isPhased(ln)) {
					String msg = "Error: VCF file is not phased";
					throw new FileParsingException(log, msg);
				}
				
				if (pos >= (start_pos + win_size)) {
					if (!cur_win.equals(new Window())) {
						cur_win.setEndIndex(index - 1);
						windows.add(cur_win);
					}
					
					if (anc_data && !anc_win.equals(new Window())) {
						ancestral.add(anc_win);
					}
					
					while (pos >= (start_pos + win_size)) {
						start_pos += win_size;
						end_pos += win_size;
					}
					
		    		cur_win = new Window(start_pos, end_pos, index);
		    		
		    		if (anc_data) {
						anc_win = new Window(start_pos, end_pos, index);
					}
				}
				
				if (ln[4].contains(",")) {
					if (validAncestralData(ln[7])) {
						
						//Fill in Window and Ancestral data
						String a0 = ln[3].toUpperCase();
						String anc_allele = getAncestralAllele(ln);
						String[] alt_alleles = ln[4].split(",");
						
						String[] alt_ids = createAltIDs(pos, a0, alt_alleles);
						
						//adding alternate SNP data
						for (int i = 0; i < alt_alleles.length; i++) {
							cur_win.addSNP(pos, a0, alt_alleles[i], alt_ids[i]);
							if (anc_data) {
								anc_win.addSNP(pos, anc_allele, "-", alt_ids[i]);
							}
						}
						//adding ref SNP data
						cur_win.addSNP(pos, alt_alleles[0], a0, alt_ids[alt_ids.length - 1]);
						if (anc_data) {
							anc_win.addSNP(pos, anc_allele, "-", alt_ids[alt_ids.length - 1]);
						}
						
						//Fill in Individual data
						for (int i = DEFAULT_COL; i < ln.length; i++) {
							//adding alternate Individual data
							String[] alleles = ln[i].split("\\|");
							
							for (int j = 0; j < alt_alleles.length; j++) {
								
								if (Integer.parseInt(alleles[0]) == j + 1) {
									individuals[i-DEFAULT_COL].addAlleleToStrand(true, true);
								}
								else {
									individuals[i-DEFAULT_COL].addAlleleToStrand(false, true);
								}
								
								if (Integer.parseInt(alleles[1]) == j + 1) {
									individuals[i-DEFAULT_COL].addAlleleToStrand(true, false);
								}
								else {
									individuals[i-DEFAULT_COL].addAlleleToStrand(false, false);
								}
								
							}
							//adding ref Individual data
							if (Integer.parseInt(alleles[0]) == 0) {
								individuals[i-DEFAULT_COL].addAlleleToStrand(true, true);
							}
							else {
								individuals[i-DEFAULT_COL].addAlleleToStrand(false, true);
							}
							
							if ( Integer.parseInt(alleles[1]) == 0) {
								individuals[i-DEFAULT_COL].addAlleleToStrand(true, false);
							}
							else {
								individuals[i-DEFAULT_COL].addAlleleToStrand(false, false);
							}
						}
						
						index += alt_alleles.length + 1;
					}
				}
				else {
					//Fill in Window and Ancestral data
					SNP s = new SNP(pos, ln[3].toUpperCase(), ln[4].toUpperCase(), ln[2]);
					boolean success = addSnpToWin(cur_win, s);
					
					if (success) {
				    	if (anc_data && validAncestralData(ln[7])) {
				    		anc_win.addSNP(pos, getAncestralAllele(ln), "-", ln[2]);
				    	}
					
				    	//Fill in Individual data
				    	for (int i = DEFAULT_COL; i < ln.length; i++) {
				    		String[] alleles = ln[i].split("\\|");
				    		individuals[i-DEFAULT_COL].addAlleleToStrand(alleles[0].charAt(0), true); 
				    		individuals[i-DEFAULT_COL].addAlleleToStrand(alleles[1].charAt(0), false); 
				    	}
				    	
				    	index++;
					}
				}
			}  
			
			cur_win.setEndIndex(index - 1);
			windows.add(cur_win);
			
			if (anc_data) {
				ancestral.add(anc_win);
			}
			
		} catch (IOException e) {
			String msg = "Error: Could not correctly read in VCF file " + file_path;
			throw new FileParsingException(log, msg);
		}
	}
	
	/*
	 * The containsSNP() function looks at position, allele0, and allele1
	 * 		-In the special case when w contains a SNP where s matches
	 * 		 position, s.allele0 == allele1 or allele0, AND s.allele1 == allele0 
	 * 		 or allele1, s is not incorporated into the window. It is considered 
	 * 		 redundant data and not realistic for further calculation.
	 */
	private boolean addSnpToWin(Window w, SNP s) {
		
		if (!w.containsSNP(s)) {
			w.addSNP(s);
			return true;
		}
		
		return false;
	}
	
	private boolean isPhased(String[] line) {
		
		for (int i = DEFAULT_COL; i < line.length; i++) {
			if (line[i].contains("/")) {
				return false;
			}
		}
		
		return true;
	}
	
	private String[] createAltIDs(int pos, String a0, String[] alleles) {
		
		String[] ids = new String[alleles.length + 1];
		
		for (int i = 0; i < alleles.length; i++) {
			ids[i] = chr + ":" + pos + ":" + a0 + ":" + alleles[i];
		}
		
		ids[alleles.length] = chr + ":" + pos + ":" + alleles[0] + ":" + a0 + "-ref"; 
		
		return ids;
	}	
	
	private boolean validAncestralData(String info) {
		
		return info.contains("AA=");
	}
	
	
	private String getAncestralAllele(String[] ln) {
		
		String[] info_arr = ln[7].split(";");
		String aa = "";
		for (int i = 0; i < info_arr.length; i++) {
			if (info_arr[i].contains("AA=")) {
				aa = info_arr[i];
			}
		}
		
		aa = aa.substring(3);
		String[] aa_arr = aa.split("\\|");
		
		if (aa_arr[0].equals("?") || aa_arr[0].length() > 1) {
			//if the ancestral indel allele is shorter than derived allele
			if (aa_arr[1].length() < aa_arr[2].length() || aa_arr[1].equals("-")) {
				
				if (ln[3].length() < ln[4].length()) {
					return ln[3].toUpperCase();
				}
				else {
					return ln[4].toUpperCase();
				}
			}
			//if the ancestral indel allele is longer than the derived allele
			else {
				
				if (ln[3].length() > ln[4].length()) {
					return ln[3].toUpperCase();
				}
				else {
					return ln[4].toUpperCase();
				}
			}
		} 
		
		return aa_arr[0].toUpperCase();
	}
	
	/**
	 * Gets the windows generated from the parsed VCF file
	 */
	public List<Window> getWindows() {
		return windows; 
	}
	
	/**
	 * Gets the ancestral data generated from the parsed VCF file 
	 */
	public List<Window> getAncestralTypes() {
		return ancestral;
	}
	
	/**
	 * Gets the individuals generated from the parsed VCF file
	 */
	public Individual[] getIndividuals() {
		return individuals; 
	}
	
}

