package envi;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import errors.IllegalInputException;
import errors.UnknownFileException;
import tools.GeneticMap;
import tools.Individual;
import tools.Log;
import tools.Window;

public class SetupDriver {

	private static int MEGABASE_CONVERSION = 1000000;
	private static String LEGEND_TYPE = ".legend";
	private static String EMF_TYPE = "emf";
	private static String ANCESTRAL_TYPE = "ancestral";
	private static String HAP_TYPE = ".hap";
	private static String VCF_TYPE = ".vcf";
	
	private boolean run_intersect = true;
	private int win_size;
	private int chr_st;
	private int chr_end;
	
	//population declarations
	private String t_pop;
	private String x_pop;
	private String o_pop;
	
	//directories and files for accessing and writing data
	private File data_dir;
	private File map_dir;
	private File wrk_dir;
	
	//target population variables (tp)
	private List<Window> tp_wins;
	private Individual[] tp_indv;
	
	//cross population variables (xp)
	private List<Window> xp_wins;
	private Individual[] xp_indv;
	
	//out-group population variable (op)
	private List<Window> op_wins;
	private Individual[] op_indv;
	
	//intersection of target and cross populations (txin)
	private List<Window> txin_wins;//for listing the windows in the intersection
	private Individual[] tp_inx_indv;//for listing the individual alleles in the intersection
	private Individual[] xp_int_indv;
	
	//intersection of cross population and out-group population (xoin)
	private List<Window> xoin_wins;
	private Individual[] xp_ino_indv;
	private Individual[] op_inx_indv;
	
	//universal variables
	private GeneticMap gm;
	private List<Window> anc_types;
	
	//progress log
	private Log log;

	public SetupDriver(HashMap<String, Object> arg_map, Log log) throws Exception {
		
		this.log = log;
		
		setArgs(arg_map);
	}
	
	public void runSetup() throws Exception {
		
		for(int i = chr_st; i <= chr_end; i++) {
			
			parseFiles(i);	
			intersectPopulations();
			createEnviFiles(i);
		}
	}
	
	private void createEnviFiles(int chr) throws Exception {
		
		System.out.println("Creating Envi Files");
		
		writeEnviVariables();
		writeWindows(chr);
	}
	
	private void writeWindows(int chr) throws IllegalInputException {
		
		//create unique window directory
		try {
		
			String all_wins_path = wrk_dir.getAbsolutePath() + File.separator + "all_wins" + File.separator;
			File all_wins = new File(all_wins_path);
			all_wins.mkdirs();
			
			//add all windows to window directory
			for(int i = 0; i < tp_wins.size() && i < txin_wins.size(); i++) {
				
				//unique to each file
				Window tp_w = tp_wins.get(i);
				Window txin_w = txin_wins.get(i);
				
				//file identifier string
				String tp_file_id = String.format("%s_win%d_chr%d_ind%d-%d.bin", 
						t_pop, i, chr, tp_w.getStIndex(), tp_w.getEndIndex());
				String txin_file_id = String.format("%sx%s_win%d_chr%d_ind%d-%d.bin", 
						t_pop, x_pop, i, chr, tp_w.getStIndex(), tp_w.getEndIndex());
				
				//add job info
				writeObj(all_wins_path + File.separator + tp_file_id, tp_w);
				writeObj(all_wins_path + File.separator + txin_file_id, txin_w);
			}
			
		} catch (FileNotFoundException e) {
			String msg = "Problem with file structure in output directory";
			throw new IllegalInputException(log, msg);
		} catch (IOException e) {
			String msg = "Could not read/write environment variables in output directory";
			throw new IllegalInputException(log, msg);
		}
		
	}
	
	private void writeEnviVariables() throws IllegalInputException {
		
		try {
			File var_dir = new File(wrk_dir.getAbsolutePath() + File.separator + "envi_var");
			var_dir.mkdirs();
			
			if(var_dir.isDirectory()) {
				String path = "";
				
				path = var_dir.getAbsoluteFile() + File.separator + "target_pop_wins.bin";
				writeObj(path, tp_wins);
				
				path = var_dir.getAbsoluteFile() + File.separator + "target_pop_indv.bin";
				writeObj(path, tp_indv);
				
				path = var_dir.getAbsoluteFile() + File.separator + "targetXcross_wins.bin";
				writeObj(path, txin_wins);
				
				path = var_dir.getAbsoluteFile() + File.separator + "targetXcross_indv.bin";
				writeObj(path, tp_inx_indv);
				
				path = var_dir.getAbsoluteFile() + File.separator + "crossXtarget_indv.bin";
				writeObj(path, xp_int_indv);
				
				path = var_dir.getAbsoluteFile() + File.separator + "crossXout_wins.bin";
				writeObj(path, xoin_wins);
				
				path = var_dir.getAbsoluteFile() + File.separator + "crossXout_indv.bin";
				writeObj(path, xp_ino_indv);
				
				path = var_dir.getAbsoluteFile() + File.separator + "outXcross_indv.bin";
				writeObj(path, op_inx_indv);
				
				path = var_dir.getAbsoluteFile() + File.separator + "anc_types.bin";
				writeObj(path, anc_types);
				
				path = var_dir.getAbsoluteFile() + File.separator + "genetic_map.bin";
				writeObj(path, gm);
			}
			else {
				String msg = "Could not create directory for environment variables; check output directory path";
				throw new IllegalInputException(log, msg);
			}
		} catch (FileNotFoundException e) {
			String msg = "Problem with file structure in output directory";
			throw new IllegalInputException(log, msg);
		} catch (IOException e) {
			e.printStackTrace();
			
			String msg = "Could not read/write environment variables in output directory";
			throw new IllegalInputException(log, msg);
		}
		
	}
	
	private void writeObj(String file_path, Object obj) throws FileNotFoundException, IOException {
		
		File file = new File(file_path);
		file.createNewFile();
		
		ObjectOutputStream oos = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(file)));
		oos.writeObject(obj);
		oos.close();
	}
	
	private void intersectPopulations() {
		
		if(run_intersect) {
			System.out.println("Running Intersections");
			log.add("\nRunning Intersections....");
			
			PopIntersector pi = new PopIntersector();
			pi.intersectCrossWithTargetPopulations(tp_wins, xp_wins, tp_indv, xp_indv);
			log.add("....");
			pi.intersectCrossWithOutgroupPopulations(op_wins, xp_wins, op_indv, op_indv);
			
			txin_wins = pi.getTargetXCrossWins();
			tp_inx_indv = pi.getTargetXCrossIndv();
			xp_int_indv = pi.getCrossXTargetIndv();
			
			xoin_wins = pi.getCrossXOutWins();
			xp_ino_indv = pi.getCrossXOutIndv();
			op_inx_indv = pi.getOutXCrossIndv();
			
			log.addLine(" complete!");
		}
		else {
			System.out.println("Skipping Intersection");
			log.addLine("\nSkipping Intersections");
			
			txin_wins = tp_wins;
			tp_inx_indv = tp_indv;
			xp_int_indv = xp_indv;
			
			xoin_wins = xp_wins;
			xp_ino_indv = xp_indv;
			op_inx_indv = op_indv;
		}
	}
	
	private void parseFiles(int chr) throws Exception {
		log.addLine("\nLoading referenced data into memory for chromosome " + chr);
		System.out.println("Loading Data");
		
		if(containsVCF(data_dir))
			runVcfParcers(chr);
		else
			runHapsLegendParcers(chr);
	}
	
	private void runHapsLegendParcers(int chr) throws Exception {
		
		//========Find Path Variables========
		String lg_tp_path = getPhasedPath(data_dir, LEGEND_TYPE, chr, t_pop);
		String ph_tp_path = getPhasedPath(data_dir, HAP_TYPE, chr, t_pop);//for target population
		
		String lg_xp_path = getPhasedPath(data_dir, LEGEND_TYPE, chr, x_pop);
		String ph_xp_path = getPhasedPath(data_dir, HAP_TYPE, chr, x_pop);//for cross population
		
		String lg_op_path = getPhasedPath(data_dir, LEGEND_TYPE, chr, o_pop);
		String ph_op_path = getPhasedPath(data_dir, HAP_TYPE, chr, o_pop);
		
		String anc_path = getAncestralPath(data_dir, chr);
		String map_path = getMapPath(map_dir, chr);
		
		if(lg_tp_path.equals(lg_xp_path) && lg_tp_path.equals(lg_op_path))
			run_intersect = false;
		
		//=======Instantiate Parsers==========
		PhasedParser tp_pp = new PhasedParser(lg_tp_path, ph_tp_path, log);
		PhasedParser xp_pp = new PhasedParser(lg_xp_path, ph_xp_path, log);
		PhasedParser op_pp = new PhasedParser(lg_op_path, ph_op_path, log);
		MapParser mp = new MapParser(map_path, log);
		
		AncestralParser ap = null;
		if(!anc_path.equals(".:MISSING:."))
			ap = new AncestralParser(anc_path, chr, log);
		
		//========Import Phased Data===========
		tp_wins = tp_pp.parseLegend(win_size);
		tp_indv = tp_pp.parsePhased(chr);
		
		xp_wins = xp_pp.parseLegend(win_size);
		xp_indv = xp_pp.parsePhased(chr);
		
		op_wins = op_pp.parseLegend(win_size);
		op_indv = op_pp.parsePhased(chr);
		
		//=======Import Environment Data========
		gm = mp.parseGeneMap();
		if(!anc_path.equals(".:MISSING:."))
			anc_types = ap.parseAncestralTypes();
		else
			anc_types = new ArrayList<Window>();
	}
	
	private void runVcfParcers(int chr) throws Exception {
		
		//========Find Path Variables==========
		String tp_vcf_path = getPhasedPath(data_dir, VCF_TYPE, chr, t_pop);
		String xp_vcf_path = getPhasedPath(data_dir, VCF_TYPE, chr, x_pop);
		String op_vcf_path = getPhasedPath(data_dir, VCF_TYPE, chr, o_pop);
		
		String map_path = getMapPath(map_dir, chr);
		
		//=====Run Parsers and Save Data=======
		VcfParser tp_vp = new VcfParser(tp_vcf_path, chr, log);
		tp_vp.parseVCF(win_size, true);
		tp_wins = tp_vp.getWindows();
		tp_indv = tp_vp.getIndividuals();
		
		VcfParser xp_vp = new VcfParser(xp_vcf_path, chr, log);
		xp_vp.parseVCF(win_size, false);
		xp_wins = xp_vp.getWindows();
		xp_indv = xp_vp.getIndividuals();
		
		VcfParser op_vp = new VcfParser(op_vcf_path, chr, log);
		op_vp.parseVCF(win_size, false);
		op_wins = op_vp.getWindows();
		op_indv = op_vp.getIndividuals();
		
		MapParser mp = new MapParser(map_path, log);
		
		//=======Import Environment Data========
		gm = mp.parseGeneMap();
		anc_types = tp_vp.getAncestralTypes();
	}
	
	private String getPhasedPath(File dir, String type, int chr, String pop) 
			throws UnknownFileException {
		
		String chr_check = "chr" + chr;
		
		String[] all_files = dir.list();
		for(int i = 0; i < all_files.length; i++) {
			
			String file_name = all_files[i];
			
			if(file_name.contains(type)
					&& file_name.contains(chr_check)
					&& file_name.contains(pop)
					&& file_name.charAt(0) != '.') {
				return dir.getAbsolutePath() + File.separator + file_name;
			}
		}
		
		throw new UnknownFileException(log, dir);
	}
	
	private String getAncestralPath(File dir, int chr) 
			throws UnknownFileException {
		
		String chr_check = "chr" + chr + "_";
		
		String[] all_files = dir.list();
		for(int i = 0; i < all_files.length; i++) {
			
			String file_name = all_files[i];
			if(file_name.contains(LEGEND_TYPE) 
					&& file_name.contains(chr_check)
					&& file_name.charAt(0) != '.'
					&& file_name.contains(ANCESTRAL_TYPE))
				return dir.getAbsolutePath() + File.separator + file_name;
			if(file_name.contains(EMF_TYPE) 
					&& file_name.contains(chr_check)
					&& file_name.charAt(0) != '.')
				return dir.getAbsolutePath() + File.separator + file_name;
		}
		
		return ".:MISSING:.";
//		String msg = "the issue is with your ancestral data";
//		throw new UnknownFileException(log, dir, msg);
	}
	
	private String getMapPath(File dir, int chr) 
			throws UnknownFileException {
		
		String chr_check = "chr" + chr + "_";
		
		String[] all_files = dir.list();
		
		for(int i = 0; i < all_files.length; i++) {
			
			String file_name = all_files[i];
			if(file_name.contains(chr_check)
					&& file_name.charAt(0) != '.')
				return dir.getAbsolutePath() + File.separator + file_name;
		}
		
		String msg = "the issue is with your map data";
		throw new UnknownFileException(log, dir, msg);
	}
	
	private boolean containsVCF(File dir) {
		String[] all_files = dir.list();
		
		for(int i = 0; i < all_files.length; i++) {
			if(all_files[i].contains(VCF_TYPE))
				return true;
		}
		
		return false;
	}
	
    private void setArgs(HashMap<String, Object> args) throws IllegalInputException, IOException {
        
        log.add("\nParameter Check");

        data_dir = (File) args.get("data_dir");
        log.add(".");

        map_dir = (File) args.get("map_dir");
        log.add(".");

        File temp_wrk_dir = new File(args.get("working_dir") + File.separator + "SelecT_workspace" + File.separator);
        int n = 1;
        while(temp_wrk_dir.exists()) {
                temp_wrk_dir = new File(args.get("working_dir") + File.separator + "SelecT_workspace" + n + File.separator);
                n++;
        }   
        temp_wrk_dir.mkdir();
        wrk_dir = new File(temp_wrk_dir.getAbsolutePath() + File.separator + "envi_files" + File.separator);
        wrk_dir.mkdirs();
        log.add(".");

        t_pop = (String)args.get("target_pop");
        log.add(".");

        x_pop = (String)args.get("cross_pop");
        log.add(".");

        o_pop = (String)args.get("out_pop");
        log.add(".");

        chr_st = (Integer) args.get("start_chr");
        log.add(".");

        chr_end = (Integer) args.get("end_chr");
        log.add(".");

        win_size = getWindowSize( (Double) args.get("win_size") );

        log.addLine(" complete!");
        System.out.println("Paramater Check Complete!");
        System.out.println("I'm praying your analysis works too...");
    }   

    private int getWindowSize(Double in) throws IllegalInputException {

        int in_size = -1; 

        try {
                double win_size_in = in * MEGABASE_CONVERSION;

                in_size = (int) win_size_in;

        } catch (NumberFormatException e) {
                String msg = "Error: Window size invalid format";
                throw new IllegalInputException(log, msg);
        }   

        if (in_size <= 0 || in_size > (100 * MEGABASE_CONVERSION)) {
                String msg = "Error: Window size declaration invalid";
                throw new IllegalInputException(log, msg);
        }   

        return in_size;
    } 
}
