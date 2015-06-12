package analysis;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;

import errors.FileParsingException;
import tools.*;

public class WindowParser {
	
	private File win_file;
	
	private WindowStats ws;
	
	private Log log;
	
	public WindowParser(Log log, File win_file, int st_pos, int end_pos) {
		
		this.log = log;
		this.win_file = win_file;
		
		ws = new WindowStats(st_pos, end_pos);
	}
	
	public WindowStats parseWindow(boolean filter) throws FileParsingException {
		
		if(filter)
			return runAnalysisParse();
		else
			return runGenericParse();
	}
	
	private WindowStats runGenericParse() throws FileParsingException {
		
		List<SNP> ihs_snps = new LinkedList<SNP>();
		List<SNP> xpehh_snps = new LinkedList<SNP>();
		List<SNP> ihh_snps = new LinkedList<SNP>();
		List<SNP> ddaf_snps = new LinkedList<SNP>();
		List<SNP> daf_snps = new LinkedList<SNP>();
		List<SNP> fst_snps = new LinkedList<SNP>();
		//List<SNP> tajd_snps = new LinkedList<SNP>();
		//List<SNP> new_snps = new LinkedList<SNP>();
		
		List<Double> ihs_stats = new LinkedList<Double>();
		List<Double> xpehh_stats = new LinkedList<Double>();
		List<Double> ihh_stats = new LinkedList<Double>();
		List<Double> ddaf_stats = new LinkedList<Double>();
		List<Double> daf_stats = new LinkedList<Double>();
		List<Double> fst_stats = new LinkedList<Double>();
		//List<Double> fajd_stats = new LinkedList<Double>();
		//List<Double> new_stats = new LinkedList<Double>();
		
		try {
			
			Scanner scan = new Scanner(win_file);
			scan.nextLine();//skip header line
			
			while(scan.hasNext()) {
				
				String[] line = scan.nextLine().split("\\s+");
				
				if(line.length != 12) {
					String msg = "Error: Window file " + win_file.getName() + " has irregular number of scores";
					throw new FileParsingException(log, msg);
				}
				
				SNP s = new SNP(Integer.parseInt(line[1]), line[0]);
				Double i_dbl = Double.parseDouble(line[2]);
				Double x_dbl = Double.parseDouble(line[3]);
				Double h_dbl = Double.parseDouble(line[4]);
				Double dd_dbl = Double.parseDouble(line[5]);
				Double d_dbl = Double.parseDouble(line[6]);
				Double f_dbl = Double.parseDouble(line[7]);
				//Double t_dbl = Double.parseDouble(line[X]);
				//Double new_dbl = Double.parseDouble(line[X]);
				
				Double un_pop_dbl = Double.parseDouble(line[8]);
				Double un_mop_dbl = Double.parseDouble(line[9]);
				Double pop_dbl = Double.parseDouble(line[10]);
				Double mop_dbl = Double.parseDouble(line[11]);
				
				if(i_dbl != Double.NaN) {
					ihs_snps.add(s);
					ihs_stats.add(i_dbl);
				}
				if(x_dbl != Double.NaN) {
					xpehh_snps.add(s);
					xpehh_stats.add(x_dbl);
				}
				if(h_dbl != Double.NaN) {
					ihh_snps.add(s);
					ihh_stats.add(h_dbl);
				}
				if(dd_dbl != Double.NaN) {
					ddaf_snps.add(s);
					ddaf_stats.add(dd_dbl);
				}
				if(d_dbl != Double.NaN) {
					daf_snps.add(s);
					daf_stats.add(d_dbl);
				}
				if(f_dbl != Double.NaN) {
					fst_snps.add(s);
					fst_stats.add(f_dbl);
				}
				//if(t_dbl != Double.NaN) {
				//	tajd_snps.add(s);
				//	tajd_stats.add(t_dbl);
				//}
				//if(new_dbl != Double.NaN) {
				//	new_snps.add(s);
				//	new_stats.add(new_dbl);
				//}
				if(un_pop_dbl != Double.NaN) 
					ws.addUnstdPopScore(s, un_pop_dbl);
				if(un_mop_dbl != Double.NaN) 
					ws.addUnstdMopScore(s, un_mop_dbl);
				if(pop_dbl != Double.NaN) 
					ws.addStdPopScore(s, pop_dbl);
				if(mop_dbl != Double.NaN) 
					ws.addStdMopScore(s, mop_dbl);
			}
			
		} catch (FileNotFoundException e) {
			String msg = "Error: Could not open data file";
			throw new FileParsingException(log, msg);
		} catch (NumberFormatException e) {
			String msg = "Error: Irregular data types in window windo file " + win_file.getName();
			throw new FileParsingException(log, msg);
		}
		
		ws.setIHS(ihs_stats, ihs_snps);
		ws.setIHH(ihh_stats, ihh_snps);
		ws.setXPEHH(xpehh_stats, xpehh_snps);
		ws.setDDAF(ddaf_stats, ddaf_snps);
		ws.setDAF(daf_stats, daf_snps);
		ws.setFST(fst_stats, fst_snps);
		//ws.setTAJD(tajd_stats, tajd_snps);
		//ws.setNEW(new_stats, new_snps);
		
		return ws;
	}
	
	private WindowStats runAnalysisParse() throws FileParsingException {
		
		List<SNP> ihs_snps = new LinkedList<SNP>();
		List<SNP> xpehh_snps = new LinkedList<SNP>();
		List<SNP> ihh_snps = new LinkedList<SNP>();
		List<SNP> ddaf_snps = new LinkedList<SNP>();
		List<SNP> daf_snps = new LinkedList<SNP>();
		List<SNP> fst_snps = new LinkedList<SNP>();
		//List<SNP> tajd_snps = new LinkedList<SNP>();
		//List<SNP> new_snps = new LinkedList<SNP>();
		
		List<Double> ihs_stats = new LinkedList<Double>();
		List<Double> xpehh_stats = new LinkedList<Double>();
		List<Double> ihh_stats = new LinkedList<Double>();
		List<Double> ddaf_stats = new LinkedList<Double>();
		List<Double> daf_stats = new LinkedList<Double>();
		List<Double> fst_stats = new LinkedList<Double>();
		//List<Double> fajd_stats = new LinkedList<Double>();
		//List<Double> new_stats = new LinkedList<Double>();
		
		try {
			
			Scanner scan = new Scanner(win_file);
			scan.nextLine();//skip header line
			
			while(scan.hasNext()) {
				
				String[] line = scan.nextLine().split("\\s+");
				
				if(line.length != 12) {
					String msg = "Error: Window file " + win_file.getName() + " has irregular number of scores";
					throw new FileParsingException(log, msg);
				}
				
				SNP s = new SNP(Integer.parseInt(line[1]), line[0]);
				Double i_dbl = Double.parseDouble(line[2]);
				Double x_dbl = Double.parseDouble(line[3]);
				Double h_dbl = Double.parseDouble(line[4]);
				Double dd_dbl = Double.parseDouble(line[5]);
				Double d_dbl = Double.parseDouble(line[6]);
				Double f_dbl = Double.parseDouble(line[7]);
				//Double t_dbl = Double.parseDouble(line[X]);
				//Double new_dbl = Double.parseDouble(line[X]);
				
				Double un_pop_dbl = Double.parseDouble(line[8]);
				Double un_mop_dbl = Double.parseDouble(line[9]);
				Double pop_dbl = Double.parseDouble(line[10]);
				Double mop_dbl = Double.parseDouble(line[11]);
				
				if(!i_dbl.equals(Double.NaN)
						&& !x_dbl.equals(Double.NaN)
						&& !h_dbl.equals(Double.NaN)
						&& !dd_dbl.equals(Double.NaN)
						&& !d_dbl.equals(Double.NaN)
						&& !f_dbl.equals(Double.NaN)
						//&& !t_dbl.equals(Double.NaN)
						//&& !new_dbl.equals(Double.NaN)
						&& !un_pop_dbl.equals(Double.NaN)
						&& !un_mop_dbl.equals(Double.NaN)
						&& !pop_dbl.equals(Double.NaN)
						&& !mop_dbl.equals(Double.NaN)) {
					
					ihs_snps.add(s);
					ihs_stats.add(i_dbl);
					xpehh_snps.add(s);
					xpehh_stats.add(x_dbl);
					ihh_snps.add(s);
					ihh_stats.add(h_dbl);
					ddaf_snps.add(s);
					ddaf_stats.add(dd_dbl);
					daf_snps.add(s);
					daf_stats.add(d_dbl);
					fst_snps.add(s);
					fst_stats.add(f_dbl);
					//tajd_snps.add(s);
					//tajd_stats.add(t_dbl);
					//new_snps.add(s);
					//new_stats.add(new_dbl);
					
					ws.addUnstdPopScore(s, un_pop_dbl);
					ws.addUnstdMopScore(s, un_mop_dbl);
					ws.addStdPopScore(s, pop_dbl);
					ws.addStdMopScore(s, mop_dbl);
				}
			}
			
			
		} catch (FileNotFoundException e) {
			String msg = "Error: Could not open data file";
			throw new FileParsingException(log, msg);
		} catch (NumberFormatException e) {
			String msg = "Error: Irregular data types in window windo file " + win_file.getName();
			throw new FileParsingException(log, msg);
		}
		
		ws.setIHS(ihs_stats, ihs_snps);
		ws.setIHH(ihh_stats, ihh_snps);
		ws.setXPEHH(xpehh_stats, xpehh_snps);
		ws.setDDAF(ddaf_stats, ddaf_snps);
		ws.setDAF(daf_stats, daf_snps);
		ws.setFST(fst_stats, fst_snps);
		//ws.setTAJD(tajd_stats, tajd_snps);
		//ws.setNEW(new_stats, new_snps);
		
		return ws;
	}
}
