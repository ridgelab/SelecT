package analysis;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import errors.FileParsingException;
import errors.IllegalInputException;
import tools.Log;
import tools.SNP;
import tools.WindowStats;

/**
 * Runs combination of windows and analysis data from statistical calculation
 *
 */
public class Combiner {

	private static String[] DEFAULT = {"i", "x", "h", "dd", "d", "f", "up", "um", "p", "m"};
	
	private int chr;
	
	private File wrk_dir;
	private File stats_dir;
	private String[] stat_str;
	
	private List<WindowStats> all_ws;
	
	private Log log;
	
	/**
	 * Constructor, instantiates the Combiner class will the input arguments from 
	 * the command line (already parsed by the SignificanceAnalyzer class)
	 * 
	 * @param arg_map					map of the arguments. Run SignificanceAnalyzer with -h option for more details 
	 * @param log						universal log file
	 * @throws IllegalInputException
	 */
	public Combiner(HashMap<String, Object> arg_map, Log log) throws IllegalInputException {
		
		this.log = log;
		
		chr = (Integer) arg_map.get("chromosome_num");
		wrk_dir = (File) arg_map.get("wrk_dir");
		
		stats_dir = new File(wrk_dir.getAbsolutePath() + File.separator + "stats_files");
		if (!stats_dir.isDirectory()) {
			String msg = "Error: Stat files directory path does not exist";
			throw new IllegalInputException(log, msg);
		}
		
		all_ws = new LinkedList<WindowStats>();
		stat_str = new String[0];
	}
	
	/**
	 * Runs combination for the windows created in statistical calculation with a 
	 * combination filter used to generate the specific list of stats desired.
	 * 
	 * @param combine_filter
	 * @throws Exception
	 */
	public void combineWindows(String combine_filter) throws Exception {
		
		log.addLine("\tRunning Window Combiner");
		System.out.println("\tRunning Window Combiner");
		
		stat_str = createStatString(combine_filter);
		
		log.addLine("Reading in files");
		System.out.println("Reading in files");
		importData(false);
	}
	
	/**
	 * Runs combination for the windows created in statistical calculation with a 
	 * default filter used to generate the list of stats printed.
	 * 
	 * @throws FileParsingException
	 */
	public void combineWindows() throws FileParsingException {
		
		log.addLine("\nRunning Window Combiner");
		System.out.println("\nRunning Window Combiner");
		
		stat_str = DEFAULT;
		log.addLine("Reading in files");
		System.out.println("Reading in files");
		importData(false);
		
	}
	
	/**
	 * Runs combination for the windows created in statistical calculation with a 
	 * default filter used to generate the list of stats printed. Filters incomplete data.
	 * 
	 * @throws FileParsingException
	 */
	public void combineAnalysisData() throws FileParsingException {
		
		log.addLine("\nRunning Window Combiner");
		System.out.println("\nRunning Window Combiner");
		
		stat_str = DEFAULT;
		log.addLine("Reading in files");
		System.out.println("Reading in files");
		importData(true);
	}
	
	/**
	 * Writes out the statistical results of SelecT analysis
	 * 
	 * @throws FileParsingException
	 */
	public void writeStats() throws FileParsingException {
		
		System.out.println("Writing combined windows");
		
		try {
			wrk_dir = new File(wrk_dir.getAbsolutePath() + File.separator + "final_out");
			if (!wrk_dir.exists()) {
				wrk_dir.mkdir();
			}
			
			File out_file = new File(wrk_dir.getAbsoluteFile() + File.separator 
					+ "combined_windows.tsv");
			int num = 1;
			while (out_file.exists()) {
				out_file = new File(wrk_dir.getAbsoluteFile() + File.separator 
						+ "combined_windows" + num + ".tsv");
				num++;
			}
			out_file.createNewFile();
			
			if (stat_str.equals(DEFAULT)) {
				simplePrint(out_file);
			}
			else {
				specificPrint(out_file);
			}
			
		} catch(IOException e) {
			String msg = "Error: There was a problem with printing to the output file in " + wrk_dir.getAbsolutePath();
			throw new FileParsingException(log, msg);
		}
		
	}
	
	public List<WindowStats> getAllStats() {
		return all_ws;
	}
	
	private void specificPrint(File out_file) throws FileNotFoundException {
		
		boolean i,x,h,dd,d,f,up,um,p,m;
		i = x = h = dd = d = f = up = um = p = m = false;
		
		if (ssContains("i")) {
			i = true;
		}
		if (ssContains("h")) {
			h = true;
		}
		if (ssContains("x")) {
			x = true;
		}
		if (ssContains("dd")) {
			dd = true;
		}
		if (ssContains("d")) {
			d = true;
		}
		if (ssContains("f")) {
			f = true;
		}
		if (ssContains("up")) {
			up = true;
		}
		if (ssContains("um")) {
			um = true;
		}
		if (ssContains("p")) {
			p = true;
		}
		if (ssContains("m")) {
			m = true;
		}
		
		PrintWriter pw = new PrintWriter(out_file);
		pw.print("snp_id\tposition");
		if (i) {
			pw.print("\tiHS");
		}
		if (x) {
			pw.print("\tXPEHH");
		}
		if (h) {
			pw.print("\tiHH");
		}
		if (dd) {
			pw.print("\tdDAF");
		}
		if (d) {
			pw.print("\tDAF");
		}
		if (f) {
			pw.print("\tFst");
		}
		if (up) {
			pw.print("\tUnstd_PoP");
		}
		if (um) {
			pw.print("\tUnstd_MoP");
		}
		if (p) {
			pw.print("\tPoP");
		}
		if (m) {
			pw.print("\tMoP");
		}
		pw.print("\n");
		
		for (int j = 0; j < all_ws.size(); j++) {
			
			WindowStats ws = all_ws.get(j);
			List<SNP> all_snps = ws.getAllSNPs();
			for (int k = 0; k < all_snps.size(); k++) {
				
				SNP snp = all_snps.get(k);
				pw.write(snp.getSnpID() + "\t" + snp.getPosition());
				if (i) {
					pw.write("\t" + ws.getIhsScore(snp).toString());
				}
				if (x) {
					pw.write("\t" + ws.getXpehhScore(snp).toString());
				}
				if (h) {
					pw.write("\t" + ws.getIhhScore(snp).toString());
				}
				if (dd) {
					pw.write("\t" + ws.getDDafScore(snp).toString());
				}
				if (d) {
					pw.write("\t" + ws.getDafScore(snp).toString());
				}
				if (f) {
					pw.write("\t" + ws.getFstScore(snp).toString());
				}
				if (up) {
					pw.write("\t" + ws.getUnstdPopScore(snp).toString());
				}
				if (um) {
					pw.write("\t" + ws.getUnstdMopScore(snp).toString());
				}
				if (p) {
					pw.write("\t" + ws.getStdPopScore(snp).toString());
				}
				if (m) {
					pw.write("\t" + ws.getStdMopScore(snp).toString());
				}
				pw.write("\n");	
			}
		}		
		pw.close();
	}
	
	private void simplePrint(File out_file) throws FileNotFoundException {
		
		PrintWriter pw = new PrintWriter(out_file);
		pw.print("snp_id\tposition\tiHS\tXPEHH\tiHH\tdDAF\tDAF\tFst"//\tTajD\tNew
				+ "\tunstd_PoP\tunstd_MoP\twin_PoP\twin_MoP\n");
		
		for (int i = 0; i < all_ws.size(); i++) {
			pw.print(all_ws.get(i));
		}
		
		pw.close();
	}
	
	private void importData(boolean filter_incomplete_data) throws FileParsingException {
		
		String[] all_file_names = stats_dir.list();
		
		for (int i = 0; i < all_file_names.length; i++) {
			
			File win_file = new File(stats_dir.getAbsolutePath() 
					+ File.separator + all_file_names[i]);
			addWindowFile(win_file, filter_incomplete_data);
		}
	}
	
	private void addWindowFile(File win_file, boolean filter_incomplete_data) throws FileParsingException {
		
		if (win_file != null && win_file.exists() 
					&& win_file.getName().charAt(0) != '.'
				&& win_file.getName().contains("chr" + chr)
				&& win_file.getName().contains("_s")
				&& win_file.getName().contains("-e")) {
			
			System.out.println("Combining Window:\t" + win_file.getName());
			
			int st_pos = getStart(win_file.getName());
			int end_pos = getEnd(win_file.getName());
			
			WindowParser wp = new WindowParser(log, win_file, st_pos, end_pos);
			WindowStats ws = wp.parseWindow(filter_incomplete_data);
			all_ws.add(ws);
		}
		else {
			log.addLine("\nWARNING: Skipping invalid file " + win_file.getName() + 
					". Please ensure that correct arguments were input for file location and chromosome number.");
		}
	}
	
	private int getEnd(String name) {
		
		int st_indx = name.indexOf("-e");
		int end_indx = name.indexOf(".tsv");
		
		return Integer.parseInt(name.substring((st_indx + 2), end_indx));
	}
	
	private int getStart(String name) {
		
		int st_indx = name.indexOf("_s");
		int end_indx = name.indexOf("-e");
		
		return Integer.parseInt(name.substring((st_indx + 2), end_indx));
	}
	
	private boolean ssContains(String stat) {
		
		for (int i = 0; i < stat_str.length; i++) {
			if (stat_str[i].equals(stat)) {
				return true;
			}
		}
		
		return false;
	}
	
	private String[] createStatString(String str) throws IllegalInputException {
		
		String[] fin_str = str.split(":");
		
		for (int i = 0; i < fin_str.length; i++) {
			
			String s = fin_str[i];
			if (!s.equals("i") && !s.equals("h") && !s.equals("x") 
					&& !s.equals("d") && !s.equals("f") && !s.equals("dd")
					&& !s.equals("up") && !s.equals("um") && !s.equals("p") 
					&& !s.equals("m")) {
				String msg = "Error: Input stat string has invalid characters or values";
				throw new IllegalInputException(log, msg);
			}
		}
		
		return fin_str;
	}
}
