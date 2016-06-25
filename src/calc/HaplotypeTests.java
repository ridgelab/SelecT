package calc;

import java.util.ArrayList;
import java.util.List;

import errors.StatsCalcException;
import tools.ExtendedHaplotype;
import tools.GeneticMap;
import tools.Individual;
import tools.SNP;
import tools.SimDist;
import tools.Window;

/**
 * Parent class for each of the haplotype-based statistical tests,
 * meaning Fst, dDAF, iHS, iHH, and XP-EHH.
 */
public abstract class HaplotypeTests {
	
	/**
	 * Standardizes the data by calculating the z-score at every position
	 * follows this equation: z = (x-mean[x])/std[x];
	 * z = standardized data point
	 * x = unstandardized data point
	 * mean[x] = mean of all data in set
	 * std[x] = standard deviation of all data in set 
	 * 
	 * @param all_unstd_data 	the unstandardized values
	 * @return 					returns the standardized values
	 */
	public static List<Double> normalizeData(List<Double> all_unstd_data) {
		
		List<Double> all_data = new ArrayList<Double>();
		
		//==========Find Mean==========
		double sum = 0.0;
		for (Double unstd_data : all_unstd_data) {
			sum += unstd_data;
		}
		
		double win_expect = sum / all_unstd_data.size();
		//=============================
		
		//===Find Standard Deviation===
		double[] all_deviations = new double[all_unstd_data.size()];
		
		for (int i = 0; i < all_unstd_data.size(); i++) {
			
			double dev = all_unstd_data.get(i) - win_expect; 
			dev = Math.pow(dev, 2);//or dev^2
			all_deviations[i] = dev;
		}
		
		double dev_sum = 0.0;
		for(int i = 0; i < all_deviations.length; i++) {
			dev_sum += all_deviations[i];
		}
		
		double dev_mean = dev_sum / all_deviations.length;
		
		double win_stddv = Math.sqrt(dev_mean);
		//==============================
		
		//==========Calc Z-Score========
		for (int i = 0; i < all_unstd_data.size(); i++) {
			all_data.add((all_unstd_data.get(i) - win_expect) / win_stddv);
		}
		//==============================
		
		return all_data;
	}
	
	public abstract void runStat() throws StatsCalcException;
	public abstract Double getScoreAtSNP(SNP s);
	public abstract Double getProbAtSNP(SNP s);
	public abstract List<SNP> getSNPs();
	public abstract List<Double> getStats();
	public abstract List<Double> getProbs();
	public abstract void printStats();

	/**
	 * Calculates the Bayesian posterior probability of each score
	 *  
	 * @param std_scores		list of scores (of type dDAF, Fst, iHH, iHS, or XPEHH)
	 * @param neut_sim			neutral simulation distances
	 * @param sel_sim			simulation distances with selection
	 * @param default_prior		True if the default probability score (1 / number of scores) should be used instead of the prior_prob.
	 * @param prior				prior probability score
	 * @return					Returns a list of bayesian score probabilities
	 */
	protected List<Double> calcScoreProbabilities(List<Double> std_scores, 
													SimDist neut_sim, 
													SimDist sel_sim,
													boolean default_prior,
													double prior) {
		
		
		List<Double> bayes_post = new ArrayList<Double>();
		double prior_prob = prior;
		
		if (default_prior) {
			prior_prob = 1 / (double) std_scores.size();
		}
				
		for (int i = 0; i < std_scores.size(); i++) {
			
			Double score = std_scores.get(i);
			
			int score_indx = neut_sim.getScoreIndex(score);
			Double post_prob = 0.0;
			
			Double neut_prob = neut_sim.getProbAtIndex(score_indx);//j
			Double sel_prob = sel_sim.getProbAtIndex(score_indx);//j
						
			if (!neut_prob.equals(Double.NaN) && !sel_prob.equals(Double.NaN)) {
				
				double cms_num = sel_prob * prior_prob;
				double cms_denom = ((sel_prob*prior_prob) + (neut_prob*(1-prior_prob)));
				post_prob = cms_num / cms_denom;
			}

			if (post_prob != 0.0) {
				bayes_post.add(post_prob);
			}
			else {
				bayes_post.add(Double.NaN);
			}
		}
		
		return bayes_post;
	}
	
	protected void setHaplotypeGroups(ExtendedHaplotype anc_eh, 
				ExtendedHaplotype der_eh, 
				Individual[] individuals,
				int snp_index, 
				SNP anc_snp, 
				SNP win_snp) {
		
		//When win_snps's a0 = ancestral type
		if (win_snp.getAllele0().equals(anc_snp.getAllele0())) {
			for (int i = 0; i < individuals.length; i++) {
				
				//adding the index of the individual with the corresponding strand (1)
				int st1_allele = individuals[i].getAlleleFromStrand(snp_index, true);
				if (st1_allele == 0) {
					anc_eh.add(i, 1);
				}
				else {
					der_eh.add(i, 1);
				}
				
				//adding the index of the individual with the corresponding strand (2)
				int st2_allele = individuals[i].getAlleleFromStrand(snp_index, false);
				if (st2_allele == 0) { 
					anc_eh.add(i, 2);
				}
				else {
					der_eh.add(i, 2);
				}
			}
		}
		//When win_snps's a1 = ancestral type
		else if (win_snp.getAllele1().equals(anc_snp.getAllele0())) {
			for (int i = 0; i < individuals.length; i++) {
				
				int st1_allele = individuals[i].getAlleleFromStrand(snp_index, true);
				if (st1_allele == 0){
					der_eh.add(i, 1);
				}
				else {
					anc_eh.add(i, 1);
				}
				
				int st2_allele = individuals[i].getAlleleFromStrand(snp_index, false);
				if (st2_allele == 0) {
					der_eh.add(i, 2);
				}
				else {
					anc_eh.add(i, 2);
				}
			}
		}
	}
	
	protected ExtendedHaplotype setHaplotypeGroup(Individual[] all_indv) {
		
		ExtendedHaplotype eh = new ExtendedHaplotype();
		
		for (int i = 0; i < all_indv.length; i++) {
			eh.add(i, 1);
			eh.add(i, 2);
		}
		
		return eh;
	}
 	
	protected double integrateEhhValues(double[] ehh_vals, int[] ehh_pos, SNP core_snp, GeneticMap geneticMap) {
		
		int core_snp_pos = core_snp.getPosition();
		int core_index = -1;
		
		for (int i = 0; i < ehh_pos.length; i++) {
			if (ehh_pos[i] == core_snp_pos) {
				core_index = i;
			}
		}
		
		
		
		double ihh_lft = 0.0;
		for (int i = core_index - 1; i >= 0; i--) {
			
			double rate = geneticMap.getRecombRate(ehh_pos[i+1], ehh_pos[i]);
			
			int dist = ehh_pos[i+1] - ehh_pos[i];
			
			double wt_dist = dist * rate;
			
			ihh_lft += wt_dist * ehh_vals[i];
		}
		
		double ihh_rt = 0.0;
		for (int i = core_index + 1; i < ehh_pos.length; i++) {
			
			double rate = geneticMap.getRecombRate(ehh_pos[i], ehh_pos[i-1]);
			
			int dist = ehh_pos[i] - ehh_pos[i-1];
			
			double wt_dist = dist * rate;
			
			ihh_rt += wt_dist * ehh_vals[i];
			
		}
		
		return ihh_lft + ihh_rt;
	}
	
	protected boolean checkValidSnpComparison(SNP core_snp, SNP anc_snp) {
		
		if (anc_snp == null || core_snp == null) {
			return false;
		}
		else if (anc_snp.getAllele0().equals(core_snp.getAllele0())) {
			return true;
		}
		else if (anc_snp.getAllele0().equals(core_snp.getAllele1())) {
			return true;
		}
		else {
			return false;
		}
	}
	
	protected SNP getAncestralSNP(SNP s, List<Window> anc_types) {
		
		for (Window w : anc_types) {			
			if (w.getStPos() < s.getPosition() && w.getEndPos() >= s.getPosition()) {
				List<SNP> win_snps = w.getSNPs();
				for (SNP anc_snp : win_snps) {
					if (s.getPosition() == anc_snp.getPosition()) {
						return anc_snp;
					}
				}
			}
		}
		
		return null;
	}
	
	protected Individual[] combineIndvArrays(Individual[] a, Individual[] b) {
		Individual[] tot = new Individual[a.length + b.length];
		
		for (int i = 0; i < a.length; i++) {
			tot[i] = a[i];
		}
		
		for (int j = 0; j < b.length; j++) {
			tot[j + a.length] = b[j];
		}
		
		return tot;
	}
	
	protected List<Double> standardizeData(List<Double> all_unstd_data) {
		//Standardizes the data by calculating the z-score at every position
		//	follows this equation: z = (x-mean[x])/std[x]; 
		//		z = standardized data point
		//		x = unstandardized data point
		//		mean[x] = mean of all data in set
		//		std[x] = standard deviation of all data in set
		
		List<Double> all_data = new ArrayList<Double>();
		
		double win_expect = findMean(all_unstd_data);
		double win_stddv = findStandardDeviation(all_unstd_data, win_expect);
		
		for(Double unstd_data : all_unstd_data) {
			all_data.add(calcZScore(unstd_data, win_expect, win_stddv));
		}
		
		return all_data;
	}
	
	private double calcZScore(double unstd_data, double expect, double stddv) {
		
		return (unstd_data - expect) / stddv;
	}
	
	protected double findStandardDeviation(List<Double> all_unstd_data, double mean) {
		
		double[] all_deviations = new double[all_unstd_data.size()];
		
		for (int i = 0; i < all_unstd_data.size(); i++) {			
			double dev = all_unstd_data.get(i) - mean; 
			dev = Math.pow(dev, 2);//or dev^2
			all_deviations[i] = dev;
		}
		
		double dev_sum = 0.0;
		for (int i = 0; i < all_deviations.length; i++) {
			dev_sum += all_deviations[i];
		}
		
		double dev_mean = dev_sum / all_deviations.length;
		
		return Math.sqrt(dev_mean);
	}
	
	protected double findMean(List<Double> all_unstd_data) {
		
		double sum = 0.0;
		
		for (Double unstd_data : all_unstd_data) {
			sum += unstd_data;
		}
		
		return sum / all_unstd_data.size();
	}
}



