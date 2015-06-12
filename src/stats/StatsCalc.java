package stats;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import calc.*;
import errors.StatsCalcException;
import tools.*;

public class StatsCalc {

	/** Hyper-Parallelized Composite of Multiple Signals (CMS) Java implementation GWAS and Local study version 1.0
	 * This program calculated the stats that are used in CMS analysis
	 * @author Hayden Smith
	 * 
	 * @param args	Required	SelecT workspace
	 * @param args	Required	Simulations directory
	 * @param args	Required	Chromosome number
	 * @param args	Required	Window number
	 * @param args	Optional	iHS absolute value probability (--ihs_abs flag)
	 * @param args	Optional	DAF cutoff (--daf_cutoff flag; default is 0.2)
	 * @param args	Optional	Prior Probability (--prior_prob flag; default without flag is number of SNPs in window)
	 */
	public static void main(String[] args) {
		
		HashMap<String, Object> arg_map = setupArgs(args);
		
		File wrk_dir = (File) arg_map.get("wrk_dir");
		Log log = new Log(Log.type.stat, wrk_dir.getName());
		
		try {
			System.out.println("Running Window:\t\t" + arg_map.get("win_num"));
			System.out.println("Chromosome:\t\t" + arg_map.get("chr"));
			System.out.println("SelecT Workspace:\t" + wrk_dir);
			
			StatsCalc sc = new StatsCalc(arg_map, log);
			sc.runStats();
			
			System.out.println("\nWindow Complete!");
			log.addLine(arg_map.get("win_num") + "\tSuccessfulRun\twindow completed without any errors");
			log.close();
			
		 } catch(OutOfMemoryError e) {
			 
			 log.addLine(arg_map.get("win_num") + "\tNoMemoryError\tinsufficient memory for running this window");
			 log.close();
			 
			 System.exit(0);
		 } 
	}
	
	private static HashMap<String, Object> setupArgs(String[] args) {
		
		ArgumentParser parser = ArgumentParsers.newArgumentParser("StatsCalc")
				.defaultHelp(true)
                .description("Run evolution statistics and calculate composite scores");
		
		//Creating required arguments
		parser.addArgument("wrk_dir").type(Arguments.fileType().verifyIsDirectory()
                .verifyCanRead()).help("SelecT workspace directory (created in phase 1)");
		
		parser.addArgument("sim_dir").type(Arguments.fileType().verifyIsDirectory()
                .verifyCanRead()).help("Directory where simulations are saved");
		
		parser.addArgument("chr").type(Integer.class).choices(Arguments.range(1, 22))
				.help("Chromosome number");
		
		parser.addArgument("win_num").type(Integer.class).help("Window number for current analysis");
		
		//Creating optional arguments
		parser.addArgument("-inon", "--nonabs_ihs").action(Arguments.storeFalse())	
				.help("Runs iHS score probabilities where large negative scores ONLY are "
						+ "associated with selection (same as CMS_local). If not included, "
						+ "large positive AND negative iHS scores equate to greater selection");
		
		//na: deflt_prior = true | # of snps in window as prior; prior_prob doesn't matter
		//--prior_prob: deflt_prior = false, prior_prob = 0.0001 | not dynamic, use prior of 1/10000, these are Broad's parameters
		//--prior_prob X: deflt_prior = false, prior_prob = X | completely custom prior probability
		parser.addArgument("-pp", "--prior_prob").nargs("?").setConst(0.0001).setDefault(-1.0)
				.type(Double.class).choices(Arguments.range(0.0, 1.0))
				.help("Sets the prior probability for bayesian score probability analysis. "
						+ "This value must be a float data-type; add a decimal point if not working. "
						+ "If not included, defaults to number of SNPs in window");
		
		parser.addArgument("-dc", "--daf_cutoff").type(Double.class).setDefault(0.2)
				.choices(Arguments.range(0.0, 0.99999))
        		.help("Sets the Derived Allele Frequency cutoff for compose score calculation. "
        				+ "This value must be a float data-type; add a decimal point if not working. "
        				+ "If not included, defaults to a frequency of 0.2. "
        				+ "Must be equal to or greater than zero and less than one");
		
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
		
		File sim_dir = (File) parsedArgs.get("sim_dir");
		String[] sim_files = sim_dir.list();
		boolean contains_neut = false;
		boolean contains_sel = false;
		for(int i = 0; i < sim_files.length; i++) {
			if(sim_files[i].equals("neutral_simulation.tsv"))
				contains_neut= true;
			if(sim_files[i].equals("selection_simulation.tsv"))
				contains_sel = true;
		}
		
		if(!contains_neut || !contains_sel) {
			
			System.out.println("Fatal error in argument parsing: see log");
			System.out.println("Could not find required simulations");
			
			Log err_log = new Log(Log.type.stat);
			err_log.addLine("Error: Could not find specific simulations for analysis");
			err_log.addLine("\t*Go to api for more information or run with -h as first parameter for help");
			err_log.addLine("\t*You will need to redo this entire step--all new data is invalid");
			
			System.exit(0);
		}
		
		return parsedArgs;
	}
	
	
	
	//Threading variables
	private static int WAIT_TIME = 50;
	
	//General information
	private int chr;
	private int win_num;
	private Window tp_win;
	private WindowStats ws;
	
	//IO Files
	private File wrk_dir;
	private File envi_dir;
	private File win_dir;
	private File sim_dir;
	
	//Analysis options
	private boolean deflt_prior;
	private boolean ihs_abs;
	private double daf_cutoff;
	private double prior_prob;
	
	//Evolution calculators
	private iHS i;
	private iHH h;
	private XPEHH x;
	private dDAF d;
	private Fst f;
	//private TajD t;
	//private NewStat new;
	
	private Log log;
	
	public StatsCalc(HashMap<String, Object> argMap, Log log) {
		
		tp_win = null;
		ws = null;
		
		this.log = log;
		
		setArgs(argMap);
	}
	
	public void runStats() {
		
		setupFiles();
		
		createCalculators();
		
		doCalculations();
		
		setWindowStats();
		
		writeOutput();
	}
	
	private void writeOutput() {
		
		System.out.println("Writing Output");
		
		try {
			
			File win_stats_file = new File(wrk_dir.getAbsolutePath() + File.separator 
					+ "win" + win_num + "_" + "chr" + chr + "_s" 
					+ tp_win.getStPos() + "-e" + tp_win.getEndPos() + ".tsv");
			win_stats_file.createNewFile();
			
			PrintWriter pw = new PrintWriter(win_stats_file);
			pw.print("snp_id\tposition\tiHS\tXPEHH\tiHH\tdDAF\tDAF\tFst"//\tTajD\tNew
					+ "\tunstd_PoP\tunstd_MoP\twin_PoP\twin_MoP\n");
			
			pw.print(ws);
			pw.close();
			
		} catch (IOException e) {
			e.printStackTrace();
			log.addLine(win_num + "\tWriteFileError\tCould not create/write the output file for this window\t" + e.getMessage());
			System.exit(0);
		}
	}
	
	private void setWindowStats() {
		
		System.out.println("Saving Evolutionary Stats and Calculating Composite Scores");
		
		ws = new WindowStats(tp_win.getStPos(), tp_win.getEndPos());
		
		ws.setIHS(i.getStats(), i.getSNPs());
		ws.setIHH(h.getStats(), h.getSNPs());
		ws.setXPEHH(x.getStats(), x.getSNPs());
		ws.setDDAF(d.getStats(), d.getSNPs());
		ws.setDAF(d.getDafStats(), d.getSNPs());
		ws.setFST(f.getStats(), f.getSNPs());
		//ws.setTAJD(t.getStats(), t.getSNPs());
		//ws.setNEW(new.getStats(), new.getSNPs());
		
		calcCompositeScores();
		
		ws.normalizeUnstdCompositeScores();
	}
	
	private void calcCompositeScores() {
		
		/*
		 * Structure of score_probs:
		 *  	[0] = iHS data
		 * 		[1] = iHH data
		 * 		[2] = Fst data
		 * 		[3] = dDAF data
		 * 		[4] = XPEHH data
		 * NOTE: new scores (like TajD) do not have simulations thus are not included in these scores
		 */
		final int NUM_TESTS = 5;
		Double[] score_probs = new Double[NUM_TESTS];
		
		List<SNP> all_snps = ws.getAllSNPs();
		for(int j = 0; j < all_snps.size(); j++) {
			
			SNP s = all_snps.get(j);
			
			score_probs[0] = i.getProbAtSNP(s);
			score_probs[1] = h.getProbAtSNP(s);
			score_probs[2] = f.getProbAtSNP(s);
			score_probs[3] = d.getProbAtSNP(s);
			score_probs[4] = x.getProbAtSNP(s);
			
			ws.addUnstdPopScore(s, productOfScores(score_probs, s));
			ws.addUnstdMopScore(s, meanOfScores(score_probs, s));
		}
	}
	
	private Double productOfScores(Double[] score_probs, SNP s) {
		
		Double prod_score = 1.0;
		
		for(int i = 0; i < score_probs.length; i++) {
			if(score_probs[i] != null && ws.getDafScore(s) >= daf_cutoff)
				prod_score = prod_score*score_probs[i];
			else
				return Double.NaN;
		}
		
		return prod_score;
	}
	
	private Double meanOfScores(Double[] score_probs, SNP s) {
		
		int tot_tests = score_probs.length;
		Double score = 0.0;
		for(int i = 0; i < score_probs.length; i++) {
			
			if(score_probs[i] != null 
					&& (daf_cutoff == 0.0 || ws.getDafScore(s) >= daf_cutoff))
				score += score_probs[i];
			else
				tot_tests--;
		}
		
		return score / tot_tests;
	}
	
	private void doCalculations() {
			
		System.out.println("Running Threads");
		
		Object lock = new Object();
		
		try {
		
			StatsThread i_thrd = new StatsThread(log, i, win_num, lock);
			Thread.sleep(WAIT_TIME);
			i_thrd.start();
			
			StatsThread h_thrd = new StatsThread(log, h, win_num, lock);
			Thread.sleep(WAIT_TIME);
			h_thrd.start();
			
			StatsThread x_thrd = new StatsThread(log, x, win_num, lock);
			Thread.sleep(WAIT_TIME);
			x_thrd.start();
			
			StatsThread d_thrd = new StatsThread(log, d, win_num, lock);
			Thread.sleep(WAIT_TIME);
			d_thrd.start();
			
			StatsThread f_thrd = new StatsThread(log, f, win_num, lock); 
			Thread.sleep(WAIT_TIME);
			f_thrd.start();
			
			//StatsThread t_thrd = new StatsThread(t, lock);
			//Thread.sleep(WAIT_TIME);
			//t_thrd.start();
			
			//StatsThread new_thrd = new StatsThread(new, lock);
			//Thread.sleep(WAIT_TIME);
			//new_thrd.start();
			
			synchronize(i_thrd, h_thrd, x_thrd, d_thrd, f_thrd);//add t_thrd
			
			i = (iHS) i_thrd.getTest();
			h = (iHH) h_thrd.getTest();
			x = (XPEHH) x_thrd.getTest();
			d = (dDAF) d_thrd.getTest();
			f = (Fst) f_thrd.getTest();
			//t = (TajD) t_thrd.getTest();
			//new = (NewStat) new_thrd.getTest();
		
		} catch (InterruptedException e) { 
			e.printStackTrace();
			log.addLine(win_num + "\tThreadingError\tThreading for this process did not complete\t" + e.getMessage());
			System.exit(0);
		}
	}
	
	private void synchronize(StatsThread i_thrd, 
								StatsThread h_thrd, 
								StatsThread x_thrd, 
								StatsThread d_thrd, 
								StatsThread f_thrd) {

		for(;;) {
			if(i_thrd.isFinished()
					&& h_thrd.isFinished()
					&& x_thrd.isFinished()
					&& d_thrd.isFinished()
					&& f_thrd.isFinished())
				break;
			else
				continue;
		}
	}
	
	@SuppressWarnings("unchecked")
	private void createCalculators() {
		
		System.out.println("Creating Calculators");
		
		try {
			//Read in saved environment variables from hard disk
			String path = "";
			
			path = getEnviFileName("target_pop_wins.bin");
			List<Window> tp_wins = (List<Window>) getObject(path);
			
			path = getEnviFileName("target_pop_indv.bin");
			Individual[] tp_indv = (Individual[]) getObject(path);
			
			path = getEnviFileName("targetXcross_wins.bin");
			List<Window> txin_wins = (List<Window>) getObject(path);
			
			path = getEnviFileName("targetXcross_indv.bin");
			Individual[] tp_inx_indv = (Individual[]) getObject(path);
			
			path = getEnviFileName("crossXtarget_indv.bin");
			Individual[] xp_int_indv = (Individual[]) getObject(path);
			
			path = getEnviFileName("crossXout_wins.bin");
			List<Window> xoin_wins = (List<Window>) getObject(path);
			
			path = getEnviFileName("crossXout_indv.bin");
			Individual[] xp_ino_indv = (Individual[]) getObject(path);
			
			path = getEnviFileName("outXcross_indv.bin");
			Individual[] op_inx_indv = (Individual[]) getObject(path);
			
			path = getEnviFileName("anc_types.bin");
			List<Window> anc_types = (List<Window>) getObject(path);
			
			path = getEnviFileName("genetic_map.bin");
			GeneticMap gm = (GeneticMap) getObject(path);
			
			path = getTargetWindowFileName();
			tp_win = (Window) getObject(path);
			
			path = getTargetXCrossWindowFileName();
			Window txin_win = (Window) getObject(path);
			
			//Create simulations
			SimulationParser sp = new SimulationParser(sim_dir);
			SimDist[] neut_sim_arr = sp.getNeutralSimulations();
			SimDist[] sel_sim_arr = sp.getSelectedSimulations();
			
			//Instantiate Selection Tests
			i = new iHS(tp_win, tp_indv, anc_types, tp_wins, gm, 
					neut_sim_arr[SimDist.IHS_TYPE], sel_sim_arr[SimDist.IHS_TYPE],
					ihs_abs, deflt_prior, prior_prob);
			
			h = new iHH(tp_win, tp_indv, anc_types, tp_wins, gm, 
					neut_sim_arr[SimDist.IHH_TYPE], sel_sim_arr[SimDist.IHH_TYPE],
					deflt_prior, prior_prob);
			
			x = new XPEHH(txin_win, txin_wins, tp_inx_indv, xp_int_indv, gm, 
					neut_sim_arr[SimDist.XPEHH_TYPE], sel_sim_arr[SimDist.XPEHH_TYPE],
					deflt_prior, prior_prob);
			
			d = new dDAF(tp_win, tp_indv, xoin_wins, xp_ino_indv, op_inx_indv, anc_types, 
					neut_sim_arr[SimDist.DDAF_TYPE], sel_sim_arr[SimDist.DDAF_TYPE],
					deflt_prior, prior_prob);
			
			f = new Fst(txin_win, tp_inx_indv, xp_int_indv, op_inx_indv, 
					neut_sim_arr[SimDist.FST_TYPE], sel_sim_arr[SimDist.FST_TYPE],
					deflt_prior, prior_prob);
			
			//t = new TajD(args0, args1, ..., argsx);
			//new = new NewStat(args0, args1, ..., argsx);
			
		} catch (IOException e) {
			e.printStackTrace();
			log.addLine(win_num + "\tReadFileError\tCould not find the correct file for proper loading of envi; check chr num and api\t" + e.getMessage());
			System.exit(0);
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
			log.addLine(win_num + "\tClassNotFoundError\tCould create the correct object instance while loading of envi\t" + e.getMessage());
			System.exit(0);
		} catch (ClassCastException e) {
			e.printStackTrace();
			log.addLine(win_num + "\tCastingError\tObject casting invalid while loading envi\t" + e.getMessage());
			System.exit(0);
		} catch (StatsCalcException e) {
			e.printStackTrace();
			log.addLine(win_num + "\t" + e.getErrorType() + "\t" + e.getMessage());
			System.exit(0);
		}
	}
	
	@SuppressWarnings("resource")
	private Object getObject(String path) throws IOException, ClassNotFoundException {
		
		File file = new File(path);
		if(!file.exists()) {
			throw new IOException();
		}
		
		ObjectInputStream ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(file)));
		return ois.readObject();
	}
	
	private String getTargetXCrossWindowFileName() throws IOException {
		
		String[] all_files = win_dir.list();
		for(int i = 0; i < all_files.length; i++) {
			
			String file_name = all_files[i];
			if(file_name.contains("chr" + chr)
					&& file_name.contains("win" + win_num + "_")
					&& file_name.contains("x"))
				return win_dir.getAbsolutePath() + File.separator + file_name;
		}
		
		throw new IOException();
	}
	
	private String getTargetWindowFileName() throws IOException {
		
		String[] all_files = win_dir.list();
		for(int i = 0; i < all_files.length; i++) {
			
			String file_name = all_files[i];
			if(file_name.contains("chr" + chr)
					&& file_name.contains("win" + win_num)
					&& !file_name.contains("x")) 
				return win_dir.getAbsolutePath() + File.separator + file_name;
		}
		
		throw new IOException();
	}
	
	private String getEnviFileName(String name) {
		return envi_dir.getAbsolutePath() + File.separator + name;
	}
	
	private void setupFiles() {
		
		envi_dir = new File(wrk_dir.getAbsoluteFile() + File.separator + "envi_files" + File.separator + "envi_var");
		if(!envi_dir.exists() && !envi_dir.isDirectory()) {
			log.addLine(win_num + "\tEnviDirError\tCould not find the evironment directory; check api for path names");
			System.exit(0);
		}
		
		win_dir = new File(wrk_dir.getAbsoluteFile() + File.separator + "envi_files" + File.separator + "all_wins");
		if(!win_dir.exists() && !win_dir.isDirectory()) {
			log.addLine(win_num + "\tWindowDirError\tCould not find the windows directory; check api for path names");
			System.exit(0);
		}
		
		wrk_dir = new File(wrk_dir.getAbsolutePath() + File.separator + "stats_files");
		if(!wrk_dir.exists())
			wrk_dir.mkdirs();
	}
	
	private void setArgs(HashMap<String, Object> args) {
		
		chr = (Integer) args.get("chr");
		win_num = (Integer) args.get("win_num");
		wrk_dir = (File) args.get("wrk_dir");
		sim_dir = (File) args.get("sim_dir");
		
		daf_cutoff = (Double) args.get("daf_cutoff");
		ihs_abs = (Boolean) args.get("nonabs_ihs");
		
		prior_prob = (Double) args.get("prior_prob");
		if(prior_prob == -1.0) 
			deflt_prior = true;
		else
			deflt_prior = false;
	}
}

class StatsThread extends Thread {

	private final Object lock;
	
	private int win_num;
	
	private Log log;
	private Thread thrd;
	private HaplotypeTests tst;
	
	volatile private boolean finished; //volatile means that some other thread could change this value
	
	StatsThread(Log log, HaplotypeTests tst, int win_num, Object lock) {
		
		this.log = log;
		this.tst = tst;
		this.win_num = win_num;
		this.lock = lock;
		
		finished = false;
		
		thrd = new Thread(this);
//		thrd.start();
	}
	
	public void start() {
		thrd.start();
	}
	
	@Override
	public void run() {
		
		try {
			tst.runStat();
		} catch (StatsCalcException e) {
			e.printStackTrace();
			log.addLine(win_num + "\t" + e.getErrorType() + "\t" + e.getMessage());
			System.exit(0);
		}	
		
		synchronized(lock) {
			finished = true;
		}
		
		thrd.interrupt();
	}
	
	public HaplotypeTests getTest() {
		return tst;
	}
	
	public boolean isFinished() {
		
		synchronized(lock) {
			return finished;
		}
	}
	
}
