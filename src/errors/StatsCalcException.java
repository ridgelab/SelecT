package errors;

/**
 * Simple exception for stats calculation errors
 */
@SuppressWarnings("serial")
public class StatsCalcException extends Exception {
	
	private String error_type = "";
	
	/**
	 * Constructor which allows an error type to be set
	 * 
	 * @param err_type
	 */
	public StatsCalcException(String err_type) {
		this.error_type = err_type;
	}
	
	/**
	 * Gets the type of error thrown
	 */
	public String getErrorType() {
		return error_type;
	}

}
