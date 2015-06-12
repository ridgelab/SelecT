package errors;

@SuppressWarnings("serial")
public class StatsCalcException extends Exception {
	
	private String error_type = "";
	
	public StatsCalcException(String err_type) {
		this.error_type = err_type;
	}
	
	public String getErrorType() {
		return error_type;
	}

}
