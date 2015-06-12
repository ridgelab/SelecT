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

public class iHS extends HaplotypeTests {
	
	//General population information
	private Window win;
	private Individual[] individuals;
	private ExtendedHaplotype anc_eh;
	private ExtendedHaplotype der_eh;
	private List<Window> all_win;
	private List<Window> anc_types;
	private GeneticMap gm;
	
	//Simulations
	private SimDist neut_sim;
	private SimDist sel_sim;
	
	//Analysis options
	private boolean ihs_abs;
	private boolean deflt_prior;
	private double prior_prob;
	
	//iHS statistic information
	private List<SNP> unused_snps;
	private List<SNP> all_iHS_snp;
	private List<Double> all_unstd_iHS;
	private List<Double> all_std_iHS;
	private List<Double> bayes_post;//bayesian posterior probability
	
	/**
	 * For setting up the environment to run the iHS statistic
	 * 
	 * @param log			universal log for progress and error output
	 * @param win			current Window within the target population (tp)
	 * @param individuals	all Individuals of the target population
	 * @param anc_types		all Ancestral types in the form of SNPs; ancestral type is a0
	 * @param all_win		all Windows in the tested region, usually the chr
	 * @param gm			Genetic Map for the tested region, usually the chr
	 */
	public iHS(Window win, 
				Individual[] individuals, 
				List<Window> anc_types,
				List<Window> all_win, 
				GeneticMap gm,
				SimDist neut_sim,
				SimDist sel_sim,
				boolean ihs_abs,
				boolean deflt_prior,
				double prior_prob){
		
		this.win = win;
		this.individuals = individuals;
		this.gm = gm;
		this.neut_sim = neut_sim;
		this.sel_sim = sel_sim;
		
		this.anc_types = anc_types;
		this.all_win = all_win;
		
		this.ihs_abs = ihs_abs;
		this.deflt_prior = deflt_prior;
		this.prior_prob = prior_prob;
		
		anc_eh = new ExtendedHaplotype();
		der_eh = new ExtendedHaplotype();
		
		unused_snps = new ArrayList<SNP>();
		all_iHS_snp = new ArrayList<SNP>();
		all_unstd_iHS = new ArrayList<Double>();
		all_std_iHS = new ArrayList<Double>();
		bayes_post = new ArrayList<Double>();
	}
	
	/**
	 * Runs the iHS statistic using the environment setup by the constructor. The
	 * below series of events spans multiple private methods but outlines the 
	 * logic for calculating iHS. This is done for every SNP in the Window (core_snp)
	 * 		-Step 1: Create extended haplotypes for all Individuals (2 per Individual)
	 * 		-Step 2: Separate the pool of haplotypes based upon whether or not they have the Ancestral allele at the core position
	 * 		-Step 3: Calculated EHH values for both groups until reaching a significantly insignificant EHH value (EHH = 0.05)
	 * 		-Step 4: Integrate EHH values from core to ends for both Ancestral and Derived groups; weight value based upon Genetic Map
	 * 		-Step 5: Calculate iHS from previously calculated iHH values
	 * 		-Step 6: Repeat and save these unstandard iHS values for all SNPs in Window
	 * 		-Step 7: Standardize all the iHS values within the Window
	 * 
	 * Note that many of these functions are extended from HaplotypeTests and 
	 * can't be found in this class.
	 * @throws StatsCalcException 
	 */
	@Override
	public void runStat() throws StatsCalcException {
		
		//Starting iHS Analysis
		int st_index = win.getStIndex();
		for(int i = 0; i < win.getSNPs().size(); i++) {
			
			Double unstd_iHS = getUnstandardizedIHS(win.getSNPs().get(i), (st_index + i));
			
			//saving the successful unstandardized iHS 
			if(unstd_iHS != null)
				all_unstd_iHS.add(unstd_iHS);

		}
		
		//calculating and saving all standardized iHS values	
		all_std_iHS = standardizeData(all_unstd_iHS);
		
		//calculates the bayesian posterior probability of each given score
//		bayes_probs = calcScoreProbabilities(all_std_iHS, neut_sim, sel_sim, true);
		calcScoreProbabilities();
		
//		printStats();
//		logRunStats();	
	}
	
	@Override
	public Double getScoreAtSNP(SNP s) {
		for(int i = 0; i < all_iHS_snp.size(); i++) {
	  		if(s.sameAs(all_iHS_snp.get(i)))
	  			return all_std_iHS.get(i);
	  	}
	  
	  	return Double.NaN;
	}
	
	@Override
	public Double getProbAtSNP(SNP s) {
	  	for(int i = 0; i < all_iHS_snp.size(); i++) {
	  		if(s.sameAs(all_iHS_snp.get(i)))
	  			return bayes_post.get(i);
	  	}
	  
	  	return null;
	}
	
	@Override
	public List<SNP> getSNPs() {
		return all_iHS_snp;
	}
	
	@Override
	public List<Double> getStats() {
		return all_std_iHS;
	}
	
	@Override
	public List<Double> getProbs() {
		return bayes_post;
	}
	 
	@Override
	public void printStats() {
		
		System.out.println("\nShowing iHS Data");
		for(int i = 0; i < all_std_iHS.size(); i++) {
			System.out.print("iHS =\t");
			System.out.print(all_iHS_snp.get(i) + "\t");
			System.out.print(all_unstd_iHS.get(i) + "\t");
			System.out.println(all_std_iHS.get(i));	
		}
	}

//	@Override
//	public void logRunStats() {
//		
//		log.addLine("Out of " + win.getSNPs().size() + " SNPs, " 
//				+ all_std_iHS.size() + " were successful and " + unused_snps.size() 
//				+ " SNPs were unsuccessful");
//	}
	
	public void printRStats() {
		
		double mean  = findMean(all_std_iHS);
		double st_dev = findStandardDeviation(all_std_iHS, mean);
		
		StringBuilder ihs_sb = new StringBuilder();
		StringBuilder pos_sb = new StringBuilder();
		
		System.out.println("\nShowing R output: iHS");
		System.out.println("\tMean:\t" + mean);
		System.out.println("\tSt Dev:\t" + st_dev);
		
		for(int i = 0; i < all_std_iHS.size(); i++) {
			
			ihs_sb.append(all_std_iHS.get(i) + ",");
			pos_sb.append(all_iHS_snp.get(i).getPosition() + ",");
		}
		System.out.println("iHS =\t" + ihs_sb.toString());
		System.out.println("Pos =\t" + pos_sb.toString());
	}
	
	private void calcScoreProbabilities() throws StatsCalcException {
		
		if(ihs_abs) {
			List<Double> all_std_iHS_abs = new ArrayList<Double>();
			for(int i = 0; i < all_std_iHS.size(); i++) {
				if(all_std_iHS.get(i) > 0)
					all_std_iHS_abs.add((-1)*all_std_iHS.get(i));
				else
					all_std_iHS_abs.add(all_std_iHS.get(i));
			}
			bayes_post = calcScoreProbabilities(all_std_iHS_abs, neut_sim, sel_sim, deflt_prior, prior_prob);
		}
		else
			bayes_post = calcScoreProbabilities(all_std_iHS, neut_sim, sel_sim, deflt_prior, prior_prob);
		
	}
	
	private Double getUnstandardizedIHS(SNP core_snp, int snp_index) {
		
		double unstd_iHS = 0.0;
		
		
		SNP anc_snp = getAncestralSNP(core_snp, anc_types);
		
		
		if(checkValidSnpComparison(core_snp, anc_snp)) {
				
			//Initial Grouping (according to ancestral or derived type)
			setHaplotypeGroups(anc_eh, der_eh, individuals, snp_index, anc_snp, core_snp);
			
			//No variance and thus no EHH pattern can be found
			if(anc_eh.size() <= 1 || der_eh.size() <= 1) {
				unused_snps.add(core_snp);
				return null;
			}
			
			//Starting EHH Analysis
			EHH anc_ehh = new EHH(win, individuals, core_snp, anc_eh, all_win);
			EHH der_ehh = new EHH(win, individuals, core_snp, der_eh, all_win);
			
			boolean significant = false; 
			
			
			//Running Ancestral Analysis
			significant = anc_ehh.calcSignificantEhhValues();
			if(!significant)
				return null;
			
			double[] ehh_values_anc = anc_ehh.getEhhValues();
			int[] ehh_pos_anc = anc_ehh.getEhhPositions();
			
			
			//Running Derived Analysis
			significant = der_ehh.calcSignificantEhhValues();
			if(!significant)
				return null;
			
			double[] ehh_values_der = der_ehh.getEhhValues();
			int[] ehh_pos_der = der_ehh.getEhhPositions();
			
			//find the area under the curve created by the EHH data
			double anc_ihh = integrateEhhValues(ehh_values_anc, ehh_pos_anc, core_snp, gm);
			double der_ihh = integrateEhhValues(ehh_values_der, ehh_pos_der, core_snp, gm);
			
			//main iHS function; unstandardized
			unstd_iHS = Math.log(anc_ihh / der_ihh);
		}
		else {
			//No ancestral allele for proper comparison
			unused_snps.add(core_snp);
			return null;
		}
		
		if(Double.isNaN(unstd_iHS) || Double.isInfinite(unstd_iHS)) {
			//Irregular iHS values
			unused_snps.add(core_snp);
			return null;
		}
		
		//saving the successful iHS SNP
		all_iHS_snp.add(core_snp);
		return unstd_iHS;
	}
	
	public void printUnusedSnps() {
		
		System.out.println("\nUnused SNPs:");
		for(SNP s : unused_snps) {
			System.out.println(s);
		}
	}
}

