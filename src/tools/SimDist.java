package tools;

import java.util.ArrayList;
import java.util.List;

/**
 * Class to hold simulation distance data
 *
 */
public class SimDist {
	
	public static final int IHS_TYPE = 0;
	public static final int IHH_TYPE = 1;
	public static final int FST_TYPE = 2;
	public static final int DDAF_TYPE = 3;
	public static final int XPEHH_TYPE = 4;
	
	private final int BIN_NUM = 60;
	
	private int up_bndry;
	private int low_bndry;
	private double total_prob;
	private int mean_indx;
	
	private List<Double> sim_vals;
	
	/**
	 * Constructor with lower and upper bounds
	 * 
	 * @param low_bndry		lower bound for simulation data
	 * @param up_bndry		upper bound for simulation data
	 */
	public SimDist(int low_bndry, int up_bndry) {
		
		this.up_bndry = up_bndry;
		this.low_bndry = low_bndry;
		this.total_prob = 0.0;
		mean_indx = -1;
		
		sim_vals = new ArrayList<Double>();
	}

	/**
	 * Adds a new simulation point
	 * 
	 * @param val	new simulation point to add
	 */
	public void addSimValue(double val) {
		
		sim_vals.add(val);
		total_prob += val;
		
		if (total_prob >= 0.5 && mean_indx == -1) {
			mean_indx = sim_vals.size();
		}
	}
	
	public List<Double> getSimVals() {
		return sim_vals;
	}
	
	public double getTotalProb() {
		return total_prob;
	}
	
	public int getMeanIndex() {
		return mean_indx;
	}
	
	/**
	 * Returns index of where the score given appears within simulation data
	 * 
	 * @param score		score value of interest
	 * @return			index where the score is found
	 */
	public int getScoreIndex(Double score) {
		
		double rng = Math.abs((double) up_bndry - (double) low_bndry);
		double bin_size = rng / BIN_NUM;
		
		int indx = (int)((score - low_bndry) / bin_size);
		
		if (indx <= 0) {
			return 0;
		}
		if (indx >= BIN_NUM-1) {
			return BIN_NUM-1;
		}
		
		return indx;
	}
	
	/**
	 * Returns probability at the given index within simulation data
	 * 
	 * @param index		index of interest
	 * @return			probability of that index
	 */
	public Double getProbAtIndex(int index) {
		
		if (index >= 0 && index < BIN_NUM) {
			return sim_vals.get(index);
		}
		else {
			return Double.NaN;
		}
	}	
}
