package errors;

import tools.Log;

/**
 * Exception thrown when arguments cannot be parsed
 */
@SuppressWarnings("serial")
public class IllegalInputException extends Exception {
	
	/**
	 * Default constructor
	 */
	public IllegalInputException() {}
	
	
	/**
	 * Constructor allowing log and error message
	 * 
	 * @param log		universal log file
	 * @param message	error message
	 */
	public IllegalInputException(Log log, String message) {
		log.addLine("\n");
		log.addLine(message);
		log.addLine("\t*Make sure all arguments are correct" );
		//TODO: Include correct wiki page reference
		log.addLine("\t*Go to wiki (https://github.com/jbelyeu/mySelecT/wiki) for parameter descriptions");
	}

}
