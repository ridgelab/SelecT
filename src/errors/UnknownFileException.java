package errors;

import java.io.File;

import tools.Log;

/**
 * Exception throw when a directory error is thrown while parsing inputs,
 * usually thrown when a file cannot be found
 */
@SuppressWarnings("serial")
public class UnknownFileException extends Exception {

	/**
	 * Simple default error when files cannot be read from the directory
	 * 
	 * @param log	universal log file
	 * @param dir	directory searched for files
	 */
	public UnknownFileException(Log log, File dir) {
		
		log.addLine("There is an error with reading files from " 
				+ dir.getAbsolutePath());
		log.addLine("\t*check that you have the correct flags in your file names");
		//TODO: Include correct wiki page reference
		log.addLine("\t*and go to wiki (https://github.com/jbelyeu/mySelecT/wiki) for parameter descriptions");
	}
	
	/**
	 * Error with message when files cannot be read from the directory
	 * 
	 * @param log	universal log file
	 * @param dir	directory searched for files
	 * @param msg	error message
	 */
	public UnknownFileException(Log log, File dir, String msg) {
		log.addLine("There is an error with reading files from " 
				+ dir.getAbsolutePath());
		log.addLine("\t*" + msg);
		log.addLine("\t*check that you have the correct flags in your file names");
		//TODO: Include correct wiki page reference 
		log.addLine("\t*and go to wiki (https://github.com/jbelyeu/mySelecT/wiki) for parameter descriptions");
	}
}
