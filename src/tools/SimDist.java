package tools;

import java.util.ArrayList;
import java.util.List;

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
	
	public SimDist(int low_bndry, int up_bndry) {
		
		this.up_bndry = up_bndry;
		this.low_bndry = low_bndry;
		this.total_prob = 0.0;
		mean_indx = -1;
		
		sim_vals = new ArrayList<Double>();
	}

	public void addSimValue(double val) {
		
		sim_vals.add(val);
		total_prob += val;
		
		if(total_prob >= 0.5 && mean_indx == -1)
			mean_indx = sim_vals.size();
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
	
	public int getScoreIndex(Double score) {
		
		double rng = Math.abs((double) up_bndry - (double) low_bndry);
		double bin_size = rng / (double) BIN_NUM;
		
		int indx = (int)((score - low_bndry) / bin_size);
		
		if(indx <= 0)
			return 0;
		if(indx >= BIN_NUM-1) 
			return BIN_NUM-1;
		
		return indx;
	}
	
	public Double getProbAtIndex(int indx) {
		
		if(indx >= 0 && indx < BIN_NUM)
			return sim_vals.get(indx);
		else
			return Double.NaN;
	}
	
	
	
//===================Artifact of a prior way of analyzing simulation distributions============================
//	public Double get1SidedProb(Double score, boolean foreword) throws StatsCalcException {
//	
//		if(sim_vals.size() != BIN_NUM) {
//			String err_type = "BinError\tBin numbers don't coincide, error in reading simulated data; redo window";
//			throw new StatsCalcException(err_type);
//		}
//	
//		int s_indx = getScoreIndex(score, up_bndry, low_bndry);
//	
//		if(foreword)
//			return calc1SidedProbAtBinForeword(s_indx);
//		else
//			return calc1SidedProbAtBinReverse(s_indx);
//		}
//
//	public Double get2SidedProb(Double score) throws StatsCalcException {
//	
//		if(sim_vals.size() != BIN_NUM) {
//			String err_type = "BinError\tBin numbers don't coincide, error in reading simulated data; redo window";
//			throw new StatsCalcException(err_type);
//		}
//	
//		int real_indx = getScoreIndex(score, up_bndry, low_bndry);
//		int other_indx = -1;
//		if(real_indx <= mean_indx) {
//			other_indx = mean_indx + (mean_indx - real_indx);
//			return calc2SidedProbAtBin(other_indx, real_indx);
//		}
//		else {
//			other_indx = mean_indx - (real_indx - mean_indx);
//			return calc2SidedProbAtBin(real_indx, other_indx);
//		}
//	}
//
//	public Double getIHS2SidedSelProb(Double score) throws StatsCalcException {
//	
//		if(sim_vals.size() != BIN_NUM) {
//			String err_type = "BinError\tBin numbers don't coincide, error in reading simulated data; redo window";
//			throw new StatsCalcException(err_type);
//		}
//		
//		if(score > 0)
//			score = -1 * Math.abs(score);
//		
//		int real_indx = getScoreIndex(score, up_bndry, low_bndry);
//		int other_indx = -1;
//		if(real_indx <= mean_indx) {
//			other_indx = mean_indx + (mean_indx - real_indx);
//			return calc2SidedProbAtBin(other_indx, real_indx);
//		}
//		else {
//			other_indx = mean_indx - (real_indx - mean_indx);
//			return calc2SidedProbAtBin(real_indx, other_indx);
//		}
//	
//	}
//
//	public Double getProb(Double score, boolean two_sided, boolean forward) throws StatsCalcException {
//	
//		if(score == Double.NaN)
//			return Double.NaN;
//		
//		if(two_sided)
//			return get2SidedProb(score);
//		else
//			return get1SidedProb(score, forward);
//	}
//
//	public Double calcTotalProb() {
//	
//		Double prob = 0.0;
//		
//		for(int i = 0; i < sim_vals.size(); i++)
//			prob += sim_vals.get(i);
//		
//		return prob;
//	}
//
//	private Double calc1SidedProbAtBinForeword(int indx) {
//	
//		Double prob = 0.0;
//		
//		for(int i = 0; i <= indx; i++)
//			prob += sim_vals.get(i);
//		
//		return prob;
//	}
//
//	private Double calc1SidedProbAtBinReverse(int indx) {
//	
//		Double prob = 0.0;
//		
//		for(int i = indx; i < sim_vals.size(); i++)
//			prob += sim_vals.get(i);
//		
//		return prob;
//	}
//
//	private Double calc2SidedProbAtBin(int up_indx, int low_indx)  {
//		
//		Double prob = 0.0;
//		
//		//Lower Bound
//		for(int i = 0; i <= low_indx; i++) 
//			prob += sim_vals.get(i);
//		
//		//Upper Bound
//		for(int i = up_indx; i <= BIN_NUM-1; i++)
//			prob += sim_vals.get(i);
//		
//		return prob;
//	}
//
//	private int getScoreIndex(Double score, int up, int dwn) {
//		
//		double rng = (double) up - (double) dwn;
//		double bin_size = rng / (double) BIN_NUM;
//		
//		for(int i = 0; i < BIN_NUM; i++) {
//			if((dwn + bin_size*i) >= score)
//				return i;
//		}
//		
//		return BIN_NUM-1;//because we are looking at indexes
//	}
	
	
}
