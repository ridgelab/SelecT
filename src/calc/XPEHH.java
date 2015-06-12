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

public class XPEHH extends HaplotypeTests {
	
	private static double SIGNIFICANT_EHH_VALUE = 0.045;
	
	//General population information
	private Window win;
	private Individual[] tp_individuals;//target population (tp) intersected with the cross population
	private Individual[] xp_individuals;//cross population (xp) intersected with the target population
	private GeneticMap gm;
	private List<Window> all_win;
	
	//Simulations
	private SimDist neut_sim;
	private SimDist sel_sim;
	
	//Analysis options
	private boolean deflt_prior;
	private double prior_prob;
	
	//XPEHH statistic information
	private List<SNP> unused_snps;
	private List<Double> all_unstd_XPEHH;
	private List<SNP> all_XPEHH_snps;
	private List<Double> all_XPEHH;
	private List<Double> bayes_probs;
	
	/**
	 * For setting up the environment to run the XPEHH statistic
	 * 
	 * @param log				universal log for progress and error output
	 * @param win				current Window for the intersection of target-cross populations
	 * @param all_win			all Windows for the intersection of target-cross populations
	 * @param tp_individuals	all Individuals of the target population after the intersection of target-cross populations
	 * @param xp_individuals	all Individuals of the cross population after the intersection of target-cross populations
	 * @param gm				Genetic Map for the tested region, usually the chr
	 */
	public XPEHH(Window win,
					List<Window> all_win,
					Individual[] tp_individuals,
					Individual[] xp_individuals,
					GeneticMap gm,
					SimDist neut_sim,
					SimDist sel_sim,
					boolean deflt_prior,
					double prior_prob) {
		
		this.win = win;
		this.tp_individuals = tp_individuals;
		this.xp_individuals = xp_individuals;
		this.gm = gm;
		this.neut_sim = neut_sim;
		this.sel_sim = sel_sim;
		
		this.all_win = all_win;
		
		this.deflt_prior = deflt_prior;
		this.prior_prob = prior_prob;
		
		unused_snps = new ArrayList<SNP>();
		all_XPEHH_snps = new ArrayList<SNP>();
		all_unstd_XPEHH = new ArrayList<Double>();
		all_XPEHH = new ArrayList<Double>();
		bayes_probs = new ArrayList<Double>();
	}

	/**
	 * Runs the XPEHH statistic using the environment setup by the constructor. The
	 * below series of evens spans multiple private methods but outlines the 
	 * logic for calculating XPEHH. This is done for every SNP in the Window (core_snp)
	 * 		-Step 1: Combine all the Individuals from both target and cross populations
	 * 		-Step 2: Calculate EHH values for the combined group of Individuals until you reach an EHH value of ~0.045
	 * 		-Step 3: Figure out the position boundaries were for the combined EHH analysis 
	 * 		-Step 4: For both populations calculate EHH values from core_snp to position boundaries
	 * 		-Step 5: Integrate EHH values from core to ends for both target and cross populations; weight value based upon Genetic Map
	 * 		-Step 6: Calculate unstandardized XPEHH from previously calculated iHH values
	 * 		-Step 7: Repeat and save these unstandard XPEHH values for all SNPs in Window
	 * 		-Step 8: Standardize the all the XPEHH values within the Window
	 * 
	 * Note that many of these functions are extended from HaplotypeTests and 
	 * can't be found in this class.
	 * @throws StatsCalcException 
	 */
	@Override
	public void runStat() throws StatsCalcException {
		
		//Starting XPEHH Analysis
		Individual[] all_indv = combineIndvArrays(tp_individuals, xp_individuals);
		
		List<SNP> win_snps = win.getSNPs();
		for(int i = 0; i < win_snps.size(); i++) {
			
			SNP core_snp = win_snps.get(i);
			
			//calculate EHH scores for the combined populations (tp with xp)
			EHH comb_ehh = getCombinedEHH(all_indv, core_snp);
			double last_ehh = comb_ehh.getLastEhhValue();
			
			if(last_ehh < SIGNIFICANT_EHH_VALUE && last_ehh > 0.0) {
				
				//for defining the EHH range to end EHH calculations
				SNP last_snp = comb_ehh.getLastSNP();
				
				//find the area under the curve created by the EHH data
				Double tp_integral = calcUnstandardEhhIntegral(core_snp, last_snp, tp_individuals);
				Double xp_integral = calcUnstandardEhhIntegral(core_snp, last_snp, xp_individuals);
				
				if(tp_integral != null && xp_integral != null) {
					
					//main XPEHH function; unstandardized
					double unstd_XPEHH = Math.log(tp_integral / xp_integral);
					
					//saving both the successful XPEHH SNP and unstandardized XPEHH value
					all_XPEHH_snps.add(core_snp);
					all_unstd_XPEHH.add(unstd_XPEHH);
				}
				else {
					//Error with calculation of the integral
					unused_snps.add(core_snp);
				}
			}
			else {
				//insignificant EHH value at the end of the boundaries
				unused_snps.add(core_snp);
			}
		}
		
		//calculating and saving all standardized XPEHH values
		all_XPEHH = standardizeData(all_unstd_XPEHH);
		
		//calculates the bayesian posterior probability of each given score
//		bayes_probs = calcScoreProbabilities(all_XPEHH, neut_sim, sel_sim, false);
		bayes_probs = calcScoreProbabilities(all_XPEHH, neut_sim, sel_sim, deflt_prior, prior_prob);
		
//		printStats();
//		logRunStats();
	} 
	
	@Override
	public Double getScoreAtSNP(SNP s) {
		for(int i = 0; i < all_XPEHH_snps.size(); i++) {
	  		if(s.sameAs(all_XPEHH_snps.get(i)))
	  			return all_XPEHH.get(i);
	  	}
	  
	  	return Double.NaN;
	}
	
	@Override
	public Double getProbAtSNP(SNP s) {
	  	for(int i = 0; i < all_XPEHH_snps.size(); i++) {
	  		if(s.sameAs(all_XPEHH_snps.get(i)))
	  			return bayes_probs.get(i);
	  	}
	  
	  	return null;
	}
	
	@Override
	public List<SNP> getSNPs() {
		return all_XPEHH_snps;
	}
	
	@Override
	public List<Double> getStats() {
		return all_XPEHH;
	}
	
	@Override
	public List<Double> getProbs() {
		return bayes_probs;
	}
	
	@Override
	public void printStats() {
		
		System.out.println("\nShowing XPEHH Data");
		for(int i = 0; i < all_XPEHH.size(); i++) {
			System.out.print("XPEHH =\t");
			System.out.print(all_XPEHH_snps.get(i) + "\t");
			System.out.print(all_unstd_XPEHH.get(i) + "\t");
			System.out.println(all_XPEHH.get(i));
		}
	}

//	@Override
//	public void logRunStats() {
//		
//		log.addLine("Out of " + win.getSNPs().size() + " SNPs, " 
//				+ all_XPEHH.size() + " were successful and " + unused_snps.size() 
//				+ " SNPs were unsuccessful");
//	}
	
	public void printRStats() {
		
		double mean  = findMean(all_XPEHH);
		double st_dev = findStandardDeviation(all_XPEHH, mean);
		
		StringBuilder xpehh_sb = new StringBuilder();
		StringBuilder pos_sb = new StringBuilder();
		
		System.out.println("\nShowing R output: XPEHH");
		System.out.println("\tMean:\t" + mean);
		System.out.println("\tSt Dev:\t" + st_dev);
		
		for(int i = 0; i < all_XPEHH.size(); i++) {
			
			xpehh_sb.append(all_XPEHH.get(i) + ",");
			pos_sb.append(all_XPEHH_snps.get(i).getPosition() + ",");
		}
		System.out.println("XPEHH =\t" + xpehh_sb.toString());
		System.out.println("Pos =\t" + pos_sb.toString());
	}
	
	private Double calcUnstandardEhhIntegral(SNP core_snp, SNP last_snp, Individual[] indv) {
		
		ExtendedHaplotype pop_eh = setHaplotypeGroup(indv);
		EHH pop_ehh = new EHH(win, indv, core_snp, pop_eh, all_win);
		boolean significant = pop_ehh.calcEhhToPosition(last_snp.getPosition());
		if(!significant)
			return null;
		
		double[] ehh_vals = pop_ehh.getEhhValues();
		int[] ehh_pos = pop_ehh.getEhhPositions();
		
		return integrateEhhValues(ehh_vals, ehh_pos, core_snp, gm);
	}
	
	private EHH getCombinedEHH(Individual[] all_indv, SNP core_snp) {
		
		ExtendedHaplotype all_eh = setHaplotypeGroup(all_indv);
		EHH comb_ehh = new EHH(win, all_indv, core_snp, all_eh, all_win);
		
		comb_ehh.calcSignificantEhhValues(SIGNIFICANT_EHH_VALUE);
		
		return comb_ehh;
	}
}
