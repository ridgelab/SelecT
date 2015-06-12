package tools;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class Log {
	
	public static enum type {envi, stat, analysis, combine};
	
	private static String ENVI_ST = "//SelecT environment setup version 1.0\n";
	private static String STAT_ST = "//SelecT multi-threaded statistical analysis version 1.0\n";
	private static String ANAL_ST = "//SelecT synthesis and calculation\n";
	private static String COMB_ST = "//SelecT raw score printout\n";
	private static String DEF_ST = "//SelecT version 1.0\n";
	private static String BEGIN1 = "//Log created to produce compilation messages and to aid in "
					+ "debugging and configuration\n";
	private static String STAT_BEGIN1 = "//Log created to produce error messages and to "
					+ "explain the nature of those errors\n";
	private static String STAT_BEGIN2 = "window\terror_type\tdescription\n";
	private static String BEGIN2 = "//Software created by BYU Bioinformatics\n\n";
	
	private File log_file;
	private PrintWriter wr;
	
	public Log() {
		
		try {
			
			log_file = createLogFile("Log", true);
			
			wr = new PrintWriter(log_file);
			addLine(DEF_ST);
			add(BEGIN1);
			add(BEGIN2);
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public Log(Log.type type) {
		
		try {
			
			if(type == Log.type.envi) {
				
				log_file = createLogFile("envi_log", true);
				wr = new PrintWriter(log_file);
				add(ENVI_ST);
				
			} else if(type == Log.type.stat) {
				
				log_file = createLogFile("stats_ERROR", false);
				wr = new PrintWriter(log_file);
				add(STAT_ST);
				
			} else if(type == Log.type.analysis) {
				
				log_file = createLogFile("analysis_log", true);
				wr = new PrintWriter(log_file);
				add(ANAL_ST);
				
			} else if(type == Log.type.combine) {
				
				log_file = createLogFile("combine_log", true);
				wr = new PrintWriter(log_file);
				add(COMB_ST);
				
			} else {
				
				log_file = createLogFile("log", true);
				wr = new PrintWriter(log_file);
				add(DEF_ST);
			}

			add(BEGIN1);
			add(BEGIN2);
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public Log(Log.type type, String name) {
		
		try {
			
			if(type != Log.type.stat) {
				//TODO: throw new error
			}
			
			log_file = new File("stats_log_" + name + ".txt");
			if(!log_file.exists()) {
				
				log_file.createNewFile();
				
				wr = new PrintWriter(new FileWriter(log_file, true));
				add(STAT_ST);
				add(STAT_BEGIN1);
				add(BEGIN2);
				add(STAT_BEGIN2);
			}
			else
				wr = new PrintWriter(new FileWriter(log_file, true));
			
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void addLine(String args) {
		wr.write(args + "\n");
		wr.flush();
	}
	
	public void add(String args) {
		wr.write(args);
		wr.flush();
	}
	
	public void close() {
		wr.close();
	}
	
	private File createLogFile(String name, boolean rename) throws IOException {
		
		int num = 1;
		File file = new File(name + ".txt");
		
		while(file.exists() && rename) {
			file = new File(name + num + ".txt");
			num++;
		}
		
		file.createNewFile();
		return file;
	}
}

