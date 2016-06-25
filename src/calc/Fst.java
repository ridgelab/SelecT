package calc;

import java.util.ArrayList;
import java.util.List;

import errors.StatsCalcException;
import tools.Individual;
import tools.SNP;
import tools.SimDist;
import tools.Window;

/**
 * Calculates the Fst (f-statistic) score as presented 
 * by Wiley and Cockerham (1984)
 */
public class Fst extends HaplotypeTests {
	
	static int NUM_POPS = 3;
	static int NUM_ALLELES = 4;
	
	//Same for all populations
	private Window win;
	
	//All three populations intersected with each other
	private Individual[] tp_indv;
	private Individual[] xp_indv;
	private Individual[] op_indv;
	
	//Simulations
	private SimDist neut_sim;
	private SimDist sel_sim;
	
	//Analysis options
	private boolean default_prior;
	private double prior_prob;
	
	//Fst statistic information
	private List<SNP> all_Fst_snps;
	private List<Double> all_Fst;
	private List<Double> bayes_probs;
	
	
	/**
	 * For setting up the environment to run the Fst statistic
	 * See supplemental material for more detail.
	 * 
	 * @param txin_win		current Window for the intersection of target-cross populations
	 * @param tp_inx_indv	all Individuals of the target population after the intersection of target-cross populations
	 * @param xp_int_indv	all Individuals of the cross population after the intersection of target-cross populations
	 * @param op_inx_indv	all Individuals of the cross population after the intersection of target-cross populations
	 * @param neut_sim		neutral simulation distances
	 * @param sel_sim		simulation distances with selection
	 * @param default_prior	True if the default probability score (1 / number of scores) should be used instead of the prior_prob.
	 * @param prior_prob	prior probability score
	 * 
	 */
	public Fst(Window txin_win, 
				Individual[] tp_inx_indv,
				Individual[] xp_int_indv,
				Individual[] op_inx_indv,
				SimDist neut_sim,
				SimDist sel_sim,
				boolean default_prior,
				double prior_prob) {
		
		this.win = txin_win;
		
		this.tp_indv = tp_inx_indv;
		this.xp_indv = xp_int_indv;
		this.op_indv = op_inx_indv;
		
		this.neut_sim = neut_sim;
		this.sel_sim = sel_sim;
		
		this.default_prior = default_prior;
		this.prior_prob = prior_prob;
		
		all_Fst_snps = new ArrayList<SNP>();
		all_Fst = new ArrayList<Double>();
		bayes_probs = new ArrayList<Double>();
	}
	
	/**
	 * Runs the Fst statistic using the environment setup by the constructor. 
	 * See supplemental for details
	 * @throws StatsCalcException 
	 */
	@Override
	public void runStat() throws StatsCalcException {
		
		//Starting Fst Analysis
		List<SNP> win_snps = win.getSNPs();
		for (int i = 0; i < win_snps.size(); i++) {
			
			SNP core_snp = win_snps.get(i);
			int index = win.getSnpIndex(core_snp);
			
			//Calculate frequencies
			double tp_size = (double) tp_indv.length*2;
			int tp_instance = getInstanceOfAllele(tp_indv, index);
			double tp_freq = tp_instance / tp_size;
			
			double xp_size = (double) xp_indv.length*2;
			int xp_instance = getInstanceOfAllele(xp_indv, index);
			double xp_freq = xp_instance / xp_size;
			
			double op_size = (double) op_indv.length*2;
			int op_instance = getInstanceOfAllele(op_indv, index);
			double op_freq = op_instance /op_size;
			
			double avg_size = (tp_size + xp_size + op_size) / NUM_POPS;
			
			//All values are averaged between the populations
			double sqrd_coef_var = getSqrdCoefficientOfVariation(tp_size,
																	xp_size,
																	op_size,
																	avg_size);
			
			double avg_freq = getAverageFrequencyOfAllele(tp_freq, 
															xp_freq, 
															op_freq,
															tp_size,
															xp_size,
															op_size,
															avg_size);
			
			double sample_var = calcSampleVariance(tp_freq, 
													xp_freq, 
													op_freq,
													avg_freq,
													tp_size,
													xp_size,
													op_size,
													avg_size);
			
			double avg_hetero_freq = getAvgHeterozygoteFreq(index,
															tp_size,
															xp_size,
															op_size,
															avg_size);
			
			//Weir and Cockerham's Fst 
			double fst = calcFst(avg_size,
									sqrd_coef_var,
									avg_freq,
									sample_var,
									avg_hetero_freq);
			
			
			all_Fst_snps.add(core_snp);
			all_Fst.add(fst);
		}
		
		//calculates the bayesian posterior probability of each given score
		bayes_probs = calcScoreProbabilities(all_Fst, neut_sim, sel_sim, default_prior, prior_prob);

	}
	
	/**
	 * Gets the Fst score at a given SNP
	 * 
	 * @param s	the SNP whose score is desired
	 * @return the Fst score 
	 */
	@Override
	public Double getScoreAtSNP(SNP s) {
		for (int i = 0; i < all_Fst_snps.size(); i++) {
	  		if (s.sameAs(all_Fst_snps.get(i))) {
	  			return all_Fst.get(i);
	  		}
	  	}
	  
	  	return Double.NaN;
	}
	
	/**
	 * Gets the bayesian probability score at a given SNP
	 * 
	 * @param s	the SNP whose score is desired
	 * @return the probability score 
	 */
	@Override
	public Double getProbAtSNP(SNP s) {
	  	for (int i = 0; i < all_Fst_snps.size(); i++) {
	  		if (s.sameAs(all_Fst_snps.get(i))) {
	  			return bayes_probs.get(i);
	  		}
	  	}
	  
	  	return null;
	}
	
	@Override
	public List<SNP> getSNPs() {
		return all_Fst_snps;
	}
	
	@Override
	public List<Double> getStats() {
		return all_Fst;
	}
	
	@Override
	public List<Double> getProbs() {
		return bayes_probs;
	}
	
	@Override
	public void printStats() {
		
		System.out.println("\nShowing Fst Data");
		for (int i = 0; i < all_Fst.size(); i++) {
			System.out.print("Fst =\t");
			System.out.print(all_Fst_snps.get(i) + "\t");
			System.out.println(all_Fst.get(i));	
		}
	}

	
	public void printRStats() {
		
		double mean  = findMean(all_Fst);
		double st_dev = findStandardDeviation(all_Fst, mean);
		
		StringBuilder fst_sb = new StringBuilder();
		StringBuilder pos_sb = new StringBuilder();
		
		System.out.println("\nShowing R output: Fst");
		System.out.println("\tMean:\t" + mean);
		System.out.println("\tSt Dev:\t" + st_dev);
		
		for (int i = 0; i < all_Fst.size(); i++) {
			
			fst_sb.append(all_Fst.get(i) + ",");
			pos_sb.append(all_Fst_snps.get(i).getPosition() + ",");
		}
		System.out.println("Fst =\t" + fst_sb.toString());
		System.out.println("Pos =\t" + pos_sb.toString());
	}
	
	private double calcFst(double avg_size,
							double sqrd_coef_var,
							double avg_freq,
							double sample_var,
							double avg_hetero_freq) {
		
		//calculating the numerator to Weir and Cockerham's Fst statistic
		double n1 = sample_var;
		double n2 = 1 / (avg_size - 1);
		double n3 = avg_freq * (1 - avg_freq);
		double n4 = sample_var * ((NUM_POPS - 1) / NUM_POPS);
		double n5 = avg_hetero_freq / NUM_ALLELES;
		
		double numerator = n1 - n2*(n3 - n4 - n5);
		
		//calculating the denominator to Weir and Cockerham's Fst statistic
		double d1 = 1 - (avg_size*sqrd_coef_var) / (NUM_POPS*(avg_size - 1));
		double d2 = avg_freq * (1 - avg_freq);
		double d3 = 1 + (((NUM_POPS - 1)*avg_size*sqrd_coef_var) / (NUM_POPS*(avg_size - 1)));
		double d4 = sample_var / NUM_POPS;
		double d5 = sqrd_coef_var / (NUM_POPS*(avg_size - 1));
		double d6 = avg_hetero_freq / NUM_ALLELES;
		
		double denominator = d1*d2 + d3*d4 + d5*d6;
		
		return numerator / denominator;
	}
	
	private double getAvgHeterozygoteFreq(int index,
											double tp_s,
											double xp_s,
											double op_s,
											double s_avg) {
		
		int tp_het_inst = getInstanceOfHeterozygosity(tp_indv, index);
		int xp_het_inst = getInstanceOfHeterozygosity(xp_indv, index);
		int op_het_inst = getInstanceOfHeterozygosity(op_indv, index);
		
		double tp_het_freq = tp_het_inst / tp_s;
		double xp_het_freq = xp_het_inst / xp_s;
		double op_het_freq = op_het_inst / op_s;
		
		double val1 = (tp_het_freq * tp_s) / (NUM_POPS * s_avg);
		double val2 = (xp_het_freq * xp_s) / (NUM_POPS * s_avg);
		double val3 = (op_het_freq * op_s) / (NUM_POPS * s_avg);
		
		double avg_het_freq = val1 + val2 + val3;
		
		return avg_het_freq;
	}
	
	private double calcSampleVariance(double f1, 
										double f2, 
										double f3,
										double f_avg,
										double s1,
										double s2,
										double s3,
										double s_avg) {
		
		double val1 = (s1 * (Math.pow((f1 - f_avg), 2))) / ((NUM_POPS - 1) * s_avg);
		double val2 = (s2 * (Math.pow((f2 - f_avg), 2))) / ((NUM_POPS - 1) * s_avg);
		double val3 = (s3 * (Math.pow((f3 - f_avg), 2))) / ((NUM_POPS - 1) * s_avg);
		
		double sample_var = val1 + val2 + val3;
		
		return sample_var;
	}
	
	private double getAverageFrequencyOfAllele(double f1, 
												double f2, 
												double f3,
												double s1,
												double s2,
												double s3,
												double s_avg) {
		
		double val1 = (f1 * s1) / (NUM_POPS * s_avg);
		double val2 = (f2 * s2) / (NUM_POPS * s_avg);
		double val3 = (f3 * s3) / (NUM_POPS * s_avg);
		
		double avg_freq = val1 + val2 + val3;
		
		return avg_freq;
	}
	
	private double getSqrdCoefficientOfVariation(double s1,
													double s2,
													double s3,
													double s_avg) {
		
		double dev1 = Math.pow((s1 - s_avg), 2);
		double dev2 = Math.pow((s2 - s_avg), 2);
		double dev3 = Math.pow((s3 - s_avg), 2);
		
		double s_sd = Math.sqrt((dev1 + dev2 + dev3) / NUM_POPS);
		
		double coef_var = s_sd / s_avg;
		
		return Math.pow((coef_var), 2);
	}
	
	private int getInstanceOfHeterozygosity(Individual[] indv, int index) {
		
		int instance = 0;
		for (int i = 0; i <indv.length; i++) {
			if (indv[i].getAlleleFromStrand(index, true) != indv[i].getAlleleFromStrand(index, false)) {
				instance++;
			}
		}
		
		return instance;
	}
	
	//returns the instance of a0 (it doesn't matter which allele because Fst is a measure of freq, not instance)
	private int getInstanceOfAllele(Individual[] indv, int index) {
		
		int instance = 0;
		for (int i = 0; i < indv.length; i++) {
			if (indv[i].getAlleleFromStrand(index, true) == 0) {
				instance++;
			}
			if (indv[i].getAlleleFromStrand(index, false) == 0) {
				instance++;
			}
		}
		
		return instance;
	}
}

