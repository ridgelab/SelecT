package envi;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import errors.FileParsingException;
import tools.*;


public class VcfParser {
	
	private final int DEFAULT_COL = 9;
	
	private String file_path;
	private int chr; 
	private List<Window> windows;
	private List<Window> ancestral;
	private Individual[] individuals;
	private Log log;
	
	public VcfParser()
	{
		file_path = null; 
		chr = 0; 
		windows = null; 
		ancestral = null; 
		individuals = null;
		log = null; 
	}
	
	public VcfParser(String file_path, int chr, Log log)
	{
		this.file_path = file_path; 
		this.chr = chr; 
		this.log = log;
		
		individuals = new Individual[0]; // Place holder
		windows = new ArrayList<Window>();
		ancestral = new ArrayList<Window>();
	}
	
	public void parseVCF(int win_size, boolean anc_data) throws FileParsingException {
		log.addLine("Importing VCF data from " +  file_path);
		
		try {
			
			Scanner scan = new Scanner(new File(file_path));
			String line = "";
			
			//Setup Individual[] that will be populated
			while(scan.hasNext()) {
				
				line = scan.nextLine();
				if(line.charAt(1) != '#') {
					
					String[] ln_arr = line.split("\\s+");
					
					individuals = new Individual[ln_arr.length - DEFAULT_COL];
					for(int i = 0; i < individuals.length; i++)
						individuals[i] = new Individual(i, chr);
					
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
		    
		    while(scan.hasNext()) {
		    	
		    	line = scan.nextLine();
				String[] ln = line.split("\\s+");
				
				pos = Integer.parseInt(ln[1]);
				
				if(!isPhased(ln)) {
					String msg = "Error: VCF file is not phased";
					throw new FileParsingException(log, msg);
				}
				
				if(pos >= (start_pos + win_size)) {
					if(!cur_win.equals(new Window())) {
						cur_win.setEndIndex(index - 1);
						windows.add(cur_win);
					}
					
					if(anc_data && !anc_win.equals(new Window()))
						ancestral.add(anc_win);
					
					while(pos >= (start_pos + win_size)) {
						start_pos += win_size;
						end_pos += win_size;
					}
					
		    		cur_win = new Window(start_pos, end_pos, index);
		    		
		    		if (anc_data) 
			    		anc_win = new Window(start_pos, end_pos, index);
				}
				
				if(ln[4].contains(",")) {
					if(validAncestralData(ln[7])) {
						
						//Fill in Window and Ancestral data
						String a0 = ln[3].toUpperCase();
						String anc_allele = getAncestralAllele(ln);
						String[] alt_alleles = ln[4].split(",");
						
						String[] alt_ids = createAltIDs(pos, a0, alt_alleles);
						
						//adding alternate SNP data
						for(int i = 0; i < alt_alleles.length; i++) {
							cur_win.addSNP(pos, a0, alt_alleles[i], alt_ids[i]);
							if(anc_data)
								anc_win.addSNP(pos, anc_allele, "-", alt_ids[i]);
						}
						//adding ref SNP data
						cur_win.addSNP(pos, alt_alleles[0], a0, alt_ids[alt_ids.length - 1]);
						if(anc_data)
							anc_win.addSNP(pos, anc_allele, "-", alt_ids[alt_ids.length - 1]);
						
						//Fill in Individual data
						for(int i = DEFAULT_COL; i < ln.length; i++) {
							//adding alternate Individual data
							String[] alleles = ln[i].split("\\|");
							
							for(int j = 0; j < alt_alleles.length; j++) {
								
								if(Integer.parseInt(alleles[0]) == j + 1) 
									individuals[i-DEFAULT_COL].addAlleleToStrand1(true);
								else
									individuals[i-DEFAULT_COL].addAlleleToStrand1(false);
								
								if(Integer.parseInt(alleles[1]) == j + 1)
									individuals[i-DEFAULT_COL].addAlleleToStrand2(true);
								else
									individuals[i-DEFAULT_COL].addAlleleToStrand2(false);
								
							}
							//adding ref Individual data
							if(Integer.parseInt(alleles[0]) == 0)
								individuals[i-DEFAULT_COL].addAlleleToStrand1(true);
							else
								individuals[i-DEFAULT_COL].addAlleleToStrand1(false);
							
							if(Integer.parseInt(alleles[1]) == 0)
								individuals[i-DEFAULT_COL].addAlleleToStrand2(true);
							else
								individuals[i-DEFAULT_COL].addAlleleToStrand2(false);
						}
						
						index += alt_alleles.length + 1;
					}
				}
				else {
					//Fill in Window and Ancestral data
					SNP s = new SNP(pos, ln[3].toUpperCase(), ln[4].toUpperCase(), ln[2]);
					boolean success = addSnpToWin(cur_win, s);
					
					if(success) {
				    	if (anc_data && validAncestralData(ln[7]))
				    		anc_win.addSNP(pos, getAncestralAllele(ln), "-", ln[2]);
					
				    	//Fill in Individual data
				    	for (int i = DEFAULT_COL; i < ln.length; i++) {
				    		
				    		String[] alleles = ln[i].split("\\|");
				    		individuals[i-DEFAULT_COL].addAlleleToStrand1(alleles[0]); 
				    		individuals[i-DEFAULT_COL].addAlleleToStrand2(alleles[1]); 
				    	}
				    	
				    	index++;
					}
				}
			}  
			
			cur_win.setEndIndex(index - 1);
			windows.add(cur_win);
			
			if(anc_data) {
				ancestral.add(anc_win);
			}
			
			//***********Testing***********
//			if(anc_data) {
//				PrintWriter pw = new PrintWriter(new File("indv.txt"));
//				pw.println("Individuals: " + individuals.length);
//				for(int i = 0; i < individuals.length; i++)
//					pw.println(individuals[i]);
//				pw.close();
//				
//				System.out.println("\n\nWindows: " + windows.size());
//				for(int i = 0; i < windows.size(); i++) 
//					System.out.println(windows.get(i));
//				
//				System.out.println("\n\nWindows: " + windows.size());
//				for(int i = 0; i < 2; i++) 
//					System.out.println(windows.get(i));
//				
//				System.out.println("\n\nAncestral:" + ancestral.size());
//				for (int i = 0; i < ancestral.size(); i++)
//					System.out.println(ancestral.get(i));
//			}
			//******************************
			
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
		
		if(!w.containsSNP(s)) {
			w.addSNP(s);
			return true;
		}
		
		return false;
	}
	
	private boolean isPhased(String[] line) {
		
		for(int i = DEFAULT_COL; i < line.length; i++) {
			if(line[i].contains("/"))
				return false;
		}
		
		return true;
	}
	
	private String[] createAltIDs(int pos, String a0, String[] alleles) {
		
		String[] ids = new String[alleles.length + 1];
		
		for(int i = 0; i < alleles.length; i++) 
			ids[i] = chr + ":" + pos + ":" + a0 + ":" + alleles[i];
		
		ids[alleles.length] = chr + ":" + pos + ":" + alleles[0] + ":" + a0 + "-ref"; 
		
		return ids;
	}	
	
	private boolean validAncestralData(String info) {
		
		//return info.matches("");
		return info.contains("AA=");
	}
	
	
	private String getAncestralAllele(String[] ln) {
		
		String[] info_arr = ln[7].split(";");
		String aa = "";
		for(int i = 0; i < info_arr.length; i++) {
			if(info_arr[i].contains("AA="))
				aa = info_arr[i];
		}
		
		aa = aa.substring(3);
		String[] aa_arr = aa.split("\\|");
		
		if(aa_arr[0].equals("?") || aa_arr[0].length() > 1) {
			//if the ancestral indel allele is shorter than derived allele
			if(aa_arr[1].length() < aa_arr[2].length() || aa_arr[1].equals("-")) {
				
				if(ln[3].length() < ln[4].length())
					return ln[3].toUpperCase();
				else
					return ln[4].toUpperCase();
			}
			//if the ancestral indel allele is longer than the derived allele
			else {
				
				if(ln[3].length() > ln[4].length())
					return ln[3].toUpperCase();
				else
					return ln[4].toUpperCase();
			}
		} 
		
		return aa_arr[0].toUpperCase();
	}
	
	public List<Window> getWindows() {
		return windows; 
	}
	
	public List<Window> getAncestralTypes() {
		return ancestral;
	}
	
	public Individual[] getIndividuals() {
		return individuals; 
	}
	
}
