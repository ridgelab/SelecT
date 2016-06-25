package errors;

import tools.Log;

/**
 * Exception thrown when there are errors within a file,
 * rather than with the filename or path. 
 */
@SuppressWarnings("serial")
public class FileParsingException extends Exception {
	
	/**
	 * Default constructor
	 */
	public FileParsingException() {}
	
	/**
	 * Constructor allowing log and error message
	 * 
	 * @param log		universal log file
	 * @param message	error message
	 */
	public FileParsingException(Log log, String message) {
		log.addLine("\n");
		log.addLine(message);
		log.addLine("\t*Problem with file structure or file formatting" );
		//TODO: Include correct wiki page reference
		log.addLine("\t*Go to wiki (https://github.com/jbelyeu/mySelecT/wiki) for more information");
	}
}
