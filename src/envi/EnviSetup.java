package envi;

import java.util.HashMap;

import envi.SetupDriver;
import errors.IllegalInputException;
import tools.Log;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;

public class EnviSetup {

	/** Hyper-Parallelized Composite of Multiple Signals (CMS) Java implementation GWAS and Local study version 1.0
	 * This program sets up the environment needed for stats calculation
	 * @author Hayden Smith
	 * 
	 * @param Required	Phased and Ancestral Data Files Directory Path (.hap/.legend -OR- .vcf)
	 * @param Required	Genetic Map directory
	 * @param Required	Start chromosome number
	 * @param Required	End chromosome number 
	 * @param Required	Target Population
	 * @param Required	Cross Population
	 * @param Optional	Outgroup Population (-op flag; default is same as cross population)
	 * @param Optional	Place for SelecT workspace (-wd flag; default is current directory)
	 * @param Optional	Window Size in Mb (-ws flag; default is 0.5Mb)
	 */
	public static void main(String[] args) {
		
		Log log = new Log(Log.type.envi);
		
		try {
			HashMap<String, Object> arg_map = setupArgs(args, log);
			
			String select_logo = " ____       _          _____ \n"
					+ "/ ___|  ___| | ___  __|_   _|\n"
					+ "\\___ \\ / _ \\ |/ _ \\/ __|| |  \n"
					+ " ___) |  __/ |  __/ (__ | |  \n"
					+ "|____/ \\___|_|\\___|\\___||_|  \n";
			System.out.println(select_logo);
		
			SetupDriver dv = new SetupDriver(arg_map, log);
			dv.runSetup();
			
			System.out.println("Setup Finished");
			log.addLine("\n\nYour environment has been set up successfully!");
			log.addLine("Phase one of SelecT is complete.");
			
		} catch (Exception e) {
			
			System.out.println("CMS Died Prematurely." 
					+ " Check log output for troubleshooting.");
			
			log.addLine("\n\nCMS Died Prematurely. Error in computation.");
			
			e.printStackTrace();
		}
	}
	
    private static HashMap<String, Object> setupArgs(String[] args, Log log)
            throws IllegalInputException {

    	ArgumentParser parser = ArgumentParsers.newArgumentParser("EnviSetup")
    				.defaultHelp(true)
                    .description("Set up the environment for SelecT");
    	
    	//Creating required arguments
    	parser.addArgument("data_dir").type(Arguments.fileType().verifyIsDirectory()
                    .verifyCanRead()).help("Directory with data files");

    	parser.addArgument("map_dir").type(Arguments.fileType().verifyIsDirectory()
                    .verifyCanRead()).help("Directory with map file");

    	parser.addArgument("start_chr").type(Integer.class).choices(Arguments.range(1, 22))
    				.help("Starting chromosome number");

    	parser.addArgument("end_chr").type(Integer.class).choices(Arguments.range(1, 22))
    				.help("Ending chromosome number. "
    						+ "Must be equal to or greater than starting chromosome number.");

    	parser.addArgument("target_pop").choices("ACB", "ASW", "BEB", "CDX", "CEU","CHB", "CHS", "CLM",
                    "ESN", "ESN", "FIN", "GBR", "GIH","GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL",
                    "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "TST", "YRI").help("Target population");

    	parser.addArgument("cross_pop").choices("ACB", "ASW", "BEB", "CDX", "CEU","CHB", "CHS", "CLM",
                    "ESN", "ESN", "FIN", "GBR", "GIH","GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL",
                    "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "TST", "YRI").help("Cross population.");

    	//Creating optional arguments
    	parser.addArgument("-op", "--out_pop").choices("ACB", "ASW", "BEB", "CDX", "CEU","CHB", "CHS",
                    "CLM", "ESN", "ESN", "FIN", "GBR", "GIH","GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL",
                    "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "TST", "YRI").help("Outgroup population. If not "
                    + "included, defaults to equal cross population");

    	parser.addArgument("-wd", "--working_dir").type(Arguments.fileType().verifyIsDirectory()
                    .verifyCanRead()).setDefault(System.getProperty("user.dir"))
                    .help("Defines the directory where SelecT will create a new working directory. "
                    		+ "If not set, defaults to the current working directory");


    	parser.addArgument("-ws", "--win_size").type(Double.class).setDefault(0.5)
                    .help("Window size in megabases. If not included, defaults to 0.5 megabases");


    	//Parsing user-inputed arguments
    	HashMap<String, Object> parsedArgs = new HashMap<String, Object>();
    	
    	//Checking to make sure input is correct
    	try {parser.parseArgs(args, parsedArgs);}
    	catch (ArgumentParserException e) {
            e.printStackTrace();
            String msg = "Error: Failed to parse arguments.\n" 
            		+ "\t*" + e.getMessage();
            throw new IllegalInputException (log, msg);
        }

	    //default end_chr to start_chr if not set
	    if (parsedArgs.get("end_chr") == null) {
	        Object start = parsedArgs.get("start_chr");
	        parsedArgs.put("end_chr", start);
	    } else if ((Integer) parsedArgs.get("end_chr") < (Integer) parsedArgs.get("start_chr")){
	    	String msg = "Error: End Chromosome comes before start chromosome";
	        throw new IllegalInputException(log, msg);
	    }
	
	    //default out_pop to cross population if not set
	    if (parsedArgs.get("out_pop") == null) {
	        Object XP = parsedArgs.get("cross_pop");
	        parsedArgs.put("out_pop", XP);
	    }
	
	    //make sure the out-group pop and the cross pop are not the same as the target pop
	    if (parsedArgs.get("cross_pop").equals(parsedArgs.get("target_pop") ) ) {
	    	String msg = "Error: Cross population is same as target population";
	        throw new IllegalInputException(log, msg);
	    } else if (parsedArgs.get("out_pop").equals(parsedArgs.get("target_pop"))) {
	    	String msg = "Error: Out-group population (-op) is same as target population";
	        throw new IllegalInputException(log, msg);
	    }
    
	    log.addLine("Working Parameters");
	    log.addLine("Data Dir:\t\t" + parsedArgs.get("data_dir"));
	    log.addLine("Map Dir:\t\t" + parsedArgs.get("map_dir"));
	    log.addLine("Envi Output Dir:\t" + parsedArgs.get("working_dir"));
	    log.addLine("Target Pop:\t\t" + parsedArgs.get("target_pop"));
	    log.addLine("Cross Pop:\t\t" + parsedArgs.get("cross_pop"));
	    log.addLine("Outgroup Pop:\t\t" + parsedArgs.get("out_pop"));
	    log.addLine("Chr Range:\t\t" + parsedArgs.get("start_chr") + "-" + parsedArgs.get("end_chr"));
	    log.addLine("Window Size:\t\t" + parsedArgs.get("win_size") + " Mb");
	
	    return parsedArgs;
    }
}
