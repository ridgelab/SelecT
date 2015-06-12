package errors;

import tools.Log;

@SuppressWarnings("serial")
public class FileParsingException extends Exception {
	
	public FileParsingException() {}
	
	public FileParsingException(Log log, String message) {
		log.addLine("\n");
		log.addLine(message);
		log.addLine("\t*Problem with file structure or file formatting" );
		log.addLine("\t*Go to api for more information");
	}
}
