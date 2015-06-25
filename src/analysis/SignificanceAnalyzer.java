package analysis;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;

import org.apache.commons.math3.special.Erf;

import errors.IllegalInputException;
import tools.Log;
import tools.SNP;
import tools.WindowStats;

public class SignificanceAnalyzer {
	
	public static void main(String[] args) {
		
		System.out.println("Starting the final phase of SelecT");
		
		Log log = new Log(Log.type.analysis);
		HashMap<String, Object> arg_map = setupArgs(args);
		
		try {
			
			if((Boolean) arg_map.get("combine_only")) {
				
				Combiner c = new Combiner(arg_map, log);
				if(arg_map.get("combine_fltr").equals(".:default:."))
					c.combineWindows();
				else
					c.combineWindows((String) arg_map.get("combine_fltr"));
				
				c.writeStats();
			}
			else {
				
				Combiner c = new Combiner(arg_map, log);
				
				if((Boolean) arg_map.get("use_incomp"))
					c.combineWindows();
				else
					c.combineAnalysisData();
				
				
				if((Boolean) arg_map.get("write_combine"))
					c.writeStats();
				
				System.out.println("\nStarting Analysis");
				boolean run_normalization = (Boolean) arg_map.get("run_norm");
				boolean ignore_mop = (Boolean) arg_map.get("ignore_mop");
				boolean ignore_pop = (Boolean) arg_map.get("ignore_pop");
				
				SignificanceAnalyzer sa = new SignificanceAnalyzer(arg_map, c.getAllStats(), log);
				sa.findSignificantSNPs(run_normalization, ignore_mop, ignore_pop);
			}
			 
		} catch (Exception e) {
			
			System.out.println("CMS Died Prematurely." 
					+ " Check log output for troubleshooting.");
			
			log.addLine("\n\nCMS Died Prematurely. Error in computation.");
			
			e.printStackTrace();
		}	
		
		System.out.println("\n\nSelecT significance analysis finished!");
		log.addLine("Significance analysis finished!");
	}
	
	private static HashMap<String, Object> setupArgs(String[] args) {
		
		ArgumentParser parser = ArgumentParsers.newArgumentParser("SignificanceAnalyzer")
				.defaultHelp(true)
                .description("Run Analysis on evolutionary statistics");
		
		//Creating required arguments
		parser.addArgument("wrk_dir").type(Arguments.fileType().verifyIsDirectory()
                .verifyCanRead()).help("SelecT workspace directory (created in phase 1)");
		
		parser.addArgument("chr").type(Integer.class).choices(Arguments.range(1, 22))
				.help("Chromosome number");
		
		//Creating optional arguments
		parser.addArgument("-co", "--combine_only").action(Arguments.storeTrue())
				.help("When present this flag only runs first half of analysis where "
						+ "windows are combined into one file");
		
		parser.addArgument("-cf", "--combine_fltr").type(String.class).setDefault(".:default:.")
				.help("Uses a specific filter for printing specific combination of stats in output file "
						+ "and can only be used in conjunction with the --combine_only flag. "
						+ "See api for definition of filter string");
		
		parser.addArgument("-p", "--p_value").type(Double.class).setDefault(0.01)
				.choices(Arguments.range(0.0, 1.0))
				.help("Sets the p-value cutoff for significance check on composite scores. "
						+ "This value must be a float data-type; add a decimal point if not working. "
						+ "If not included, defaults to 0.01");
		
		parser.addArgument("-wc", "--write_combine").action(Arguments.storeTrue())
				.help("During full analysis this flag creates a separate file "
						+ "complete with all unfiltered data");
		
		parser.addArgument("-rn", "--run_norm").action(Arguments.storeTrue())
				.help("Normalizes MoP and PoP probabilities across whole dataset instead of "
						+ "one window before finding significance");
		
		parser.addArgument("-ui", "--use_incomp").action(Arguments.storeTrue())
				.help("Use incomplete data when analyzing MoP scores");
		
		parser.addArgument("-im", "--ignore_mop").action(Arguments.storeTrue())
				.help("Ignores all MoP scores and finds significance based upon PoP only; "
						+ "when both --ignore flags are present significance is determined "
						+ "by finding either MoP or PoP significance");
		
		parser.addArgument("-ip", "--ignore_pop").action(Arguments.storeTrue())
				.help("Ignores all PoP scores and finds significance based upon MoP only; "
						+ "when both --ignore flags are present significance is determined "
						+ "by finding either MoP or PoP significance");
		
		//Parsing user-inputed arguments
		HashMap<String, Object> parsedArgs = new HashMap<String, Object>();
		
		//Checking to make sure input is correct
		try {parser.parseArgs(args, parsedArgs);}
    	catch (ArgumentParserException e) {
    		System.out.println("Fatal error in argument parsing: see log");
            e.printStackTrace();
            
            Log err_log = new Log(Log.type.stat);
            err_log.addLine("Error: Failed to parse arguments"); 
            err_log.addLine("\t*" + e.getMessage());
            err_log.addLine("\t*Go to api for more information or run with -h as first parameter for help");
			err_log.addLine("\t*You will need to redo this entire step--all new data is invalid");
			
			System.exit(0);
        }
		return parsedArgs;
	}
	
	
	private double sig_score;
	private File out_file;
	private List<WindowStats> all_ws;
	
	private Log log;
	
	public SignificanceAnalyzer(HashMap<String, Object> arg_map, List<WindowStats> all_ws, Log log) {
		
		this.log = log;
		this.all_ws = all_ws;
		
		sig_score = pToZ((Double) arg_map.get("p_value"));
		System.out.println("Using p-value of " + (Double) arg_map.get("p_value") 
				+ " and significance score of " + sig_score);
		
		File wrk_dir = (File) arg_map.get("wrk_dir");
		wrk_dir = new File(wrk_dir.getAbsolutePath() + File.separator + "final_out");
		if(!wrk_dir.exists())
			wrk_dir.mkdir();
		out_file = new File(wrk_dir.getAbsoluteFile() + File.separator + "significant_loci.tsv");
		int num = 1;
		while(out_file.exists()) {
			out_file = new File(wrk_dir.getAbsoluteFile() + File.separator 
					+ "significant_loci" + num + ".tsv");
			num++;
		}
	}
	
	public void findSignificantSNPs(boolean run_normalization, boolean ignore_mop, boolean ignore_pop) 
			throws IllegalInputException {
		
		if(run_normalization) 
			all_ws = normalizeAllWindows(all_ws);
		
		try {
			
			System.out.println("Extracting significant loci");
			log.addLine("Extracting significant loci");
			
			out_file.createNewFile();
			PrintWriter pw = new PrintWriter(out_file);
			pw.print("snp_id\tposition\tiHS\tXPEHH\tiHH\tdDAF\tDAF\tFst"//\tTajD\tNew
					+ "\tunstd_PoP\tunstd_MoP\twin_PoP\twin_MoP\n");
			
			Collections.sort(all_ws);
			for(int i = 0; i < all_ws.size(); i++) {
				
				WindowStats ws = all_ws.get(i);
				List<SNP> ws_snps = ws.getAllSNPs();
				
				for(int j = 0; j < ws_snps.size(); j++) {
					
					SNP s = ws_snps.get(j);
					
					if(ignore_mop && ignore_pop) {
						if(ws.getStdPopScore(s) >= sig_score
								|| ws.getStdMopScore(s) >= sig_score) 
							pw.print(ws.printSNP(s));
					} else if(ignore_mop) {
						if(ws.getStdPopScore(s) >= sig_score) 
							pw.print(ws.printSNP(s));
					} else if(ignore_pop) {
						if(ws.getStdMopScore(s) >= sig_score)
							pw.print(ws.printSNP(s));
					} else {
						if(ws.getStdPopScore(s) >= sig_score
								&& ws.getStdMopScore(s) >= sig_score)
							pw.print(ws.printSNP(s));
					}
				}
			}
			
			pw.close();
			
		} catch (FileNotFoundException e) {
			String msg = "Error: Could not create output file for final analysis";
			throw new IllegalInputException(log, msg);
		} catch (IOException e) {
			String msg = "Error: Could not create output file for final analysis";
			throw new IllegalInputException(log, msg);
		}
	}
	
	private List<WindowStats> normalizeAllWindows(List<WindowStats> all_stats) {
		
		System.out.println("Running Normalization");
		log.addLine("Running Normalization");
		
		WindowStats comb_ws = new WindowStats(all_stats.get(0).getStPos(), 
				all_stats.get(all_stats.size()-1).getEndPos());
		
		for(int i = 0; i < all_stats.size(); i++) {
			
			WindowStats ws = all_stats.get(i);
			
			comb_ws.addUnstdPoP(ws.getUnstdPoP());
			comb_ws.addUnstdMoP(ws.getUnstdMoP());
		}
		
		comb_ws.normalizeUnstdCompositeScores();
		
		TreeMap<SNP, Double> std_pop = comb_ws.getStdPoP();
		Set<SNP> all_pop_snps = std_pop.keySet();
		for(SNP s : all_pop_snps) {
			for(int i = 0; i < all_stats.size(); i++) {
				if(all_stats.get(i).containsSNP(s))
					all_stats.get(i).addStdPopScore(s, std_pop.get(s));
			}
		}
		
		TreeMap<SNP, Double> std_mop = comb_ws.getStdMoP();
		Set<SNP> all_mop_snps = std_pop.keySet();
		for(SNP s : all_mop_snps) {
			for(int i = 0; i < all_stats.size(); i++) {
				if(all_stats.get(i).containsSNP(s))
					all_stats.get(i).addStdMopScore(s, std_mop.get(s));
			}
		}
		
		return all_stats;
	}
	
	private double pToZ(double p) {
	    double z = Math.sqrt(2) * Erf.erfcInv(2*p);
	    return(z);
	}
}
