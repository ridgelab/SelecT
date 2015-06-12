package errors;

import java.io.File;

import tools.Log;

@SuppressWarnings("serial")
public class UnknownFileException extends Exception {

	public UnknownFileException(Log log, File dir) {
		
		log.addLine("There is an error with reading files from " 
				+ dir.getAbsolutePath());
		log.addLine("\t*check that you have the correct flags in your file names");
		log.addLine("\t*and go to api for parameter descriptions");
	}
	
	public UnknownFileException(Log log, File dir, String msg) {
		log.addLine("There is an error with reading files from " 
				+ dir.getAbsolutePath());
		log.addLine("\t*" + msg);
		log.addLine("\t*check that you have the correct flags in your file names");
		log.addLine("\t*and go to api for parameter descriptions");
	}
}
