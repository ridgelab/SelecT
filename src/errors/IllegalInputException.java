package errors;

import tools.Log;

@SuppressWarnings("serial")
public class IllegalInputException extends Exception {
	
	public IllegalInputException() {}
	
	public IllegalInputException(Log log, String message) {
		log.addLine("\n");
		log.addLine(message);
		log.addLine("\t*Make sure all arguments are correct" );
		log.addLine("\t*Go to api for parameter descriptions");
	}

}
