package calc;

import java.util.ArrayList;
import java.util.List;

import errors.StatsCalcException;
import tools.Individual;
import tools.SNP;
import tools.SimDist;
import tools.Window;

public class dDAF extends HaplotypeTests {
	
	//General population information
	private Window tp_win;
	private Individual[] tp_indv;
	private List<Window> xoin_wins;
	private Individual[] xp_ino_indv;//previously intersected with op
	private Individual[] op_inx_indv;//previously intersected with xp
	private List<Window> anc_types;
	
	//Simulations
	private SimDist neut_sim;
	private SimDist sel_sim;
	
	//Analysis options
	private boolean deflt_prior;
	private double prior_prob;
	
	//DAF statistic information
	private List<SNP> unused_snps;
	private List<SNP> all_delta_DAF_snps;
	private List<Double> all_DAF;
	private List<Double> all_delta_DAF;
	private List<Double> bayes_probs;
	
	public dDAF(Window tp_win, 
				Individual[] tp_indv,
				List<Window> xoin_wins,
				Individual[] xp_ino_indv,
				Individual[] op_inx_indv,
				List<Window> anc_types,
				SimDist neut_sim,
				SimDist sel_sim,
				boolean deflt_prior,
				double prior_prob){
		
		this.tp_win = tp_win;
		this.tp_indv = tp_indv;
		
		this.xoin_wins = xoin_wins;
		this.xp_ino_indv = xp_ino_indv;
		this.op_inx_indv = op_inx_indv;
		this.neut_sim = neut_sim;
		this.sel_sim = sel_sim;
		
		this.anc_types = anc_types;
		
		this.deflt_prior = deflt_prior;
		this.prior_prob = prior_prob;
		
		unused_snps = new ArrayList<SNP>();
		all_delta_DAF_snps = new ArrayList<SNP>();
		all_DAF = new ArrayList<Double>();
		all_delta_DAF = new ArrayList<Double>();
		bayes_probs = new ArrayList<Double>();
	}
	
	@Override
	public void runStat() throws StatsCalcException {
		
		//Starting dDAF Analysis
		Individual[] all_xo_indv = combineIndvArrays(xp_ino_indv, op_inx_indv);
		
		List<SNP> win_snps = tp_win.getSNPs();
		for(int i = 0; i < win_snps.size(); i++) {
			
			SNP core_snp = win_snps.get(i);
			SNP anc_snp = getAncestralSNP(core_snp, anc_types);
			
			if(checkValidSnpComparison(core_snp, anc_snp)) {
				
				//Get proper intersected window for cross and outgroup populations
				Window xo_win = getEquivalentWindow(xoin_wins, tp_win);
				if(xo_win != null && xo_win.containsSNP(core_snp)) {
				
					SNP xo_snp = xo_win.getSNP(core_snp.getPosition(), 
												core_snp.getAllele0(), 
												core_snp.getAllele1());
					
					int tp_indx = tp_win.getSnpIndex(core_snp);
					int xo_indx = xo_win.getSnpIndex(xo_snp);
				
					int tp_instance_der = getInstanceOfDerivedAllele(tp_indv,
														core_snp,
														anc_snp,
														tp_indx);
					int xo_instance_der = getInstanceOfDerivedAllele(all_xo_indv,
														xo_snp,
														anc_snp,
														xo_indx);
					
					double daf_tp = (double) tp_instance_der / ((double) tp_indv.length*2);
					double daf_xo = (double) xo_instance_der / ((double) all_xo_indv.length*2);
					
					double delta_daf = daf_xo - daf_tp;
					
					all_delta_DAF_snps.add(core_snp);
					all_DAF.add(daf_tp);
					all_delta_DAF.add(delta_daf);
				}
			}
			else {
				//No ancestral allele for proper comparison
				unused_snps.add(core_snp);
			}
		}
		
		//calculates the bayesian posterior probability of each given score
//		bayes_probs = calcScoreProbabilities(all_delta_DAF, neut_sim, sel_sim, true);
		bayes_probs = calcScoreProbabilities(all_delta_DAF, neut_sim, sel_sim, deflt_prior, prior_prob);
		
//		printStats();
//		logRunStats();
	}
	
	@Override
	public Double getScoreAtSNP(SNP s) {
		for(int i = 0; i < all_delta_DAF_snps.size(); i++) {
	  		if(s.sameAs(all_delta_DAF_snps.get(i)))
	  			return all_delta_DAF.get(i);
	  	}
	  
	  	return Double.NaN;
	}
	
	@Override
	public Double getProbAtSNP(SNP s) {
	  	for(int i = 0; i < all_delta_DAF_snps.size(); i++) {
	  		if(s.sameAs(all_delta_DAF_snps.get(i)))
	  			return bayes_probs.get(i);
	  	}
	  
	  	return null;
	}
	
	@Override
	public List<SNP> getSNPs() {
		return all_delta_DAF_snps;
	}
	
	@Override
	public List<Double> getStats() {
		return all_delta_DAF;
	}
	
	@Override
	public List<Double> getProbs() {
		return bayes_probs;
	}
	
	@Override
	public void printStats() {
		
		System.out.println("\nShowing DAF Data");
		for(int i = 0; i < all_delta_DAF.size(); i++) {
			System.out.print("DAF =\t");
			System.out.print(all_delta_DAF_snps.get(i) + "\t");
			System.out.print(all_DAF.get(i) + "\t");
			System.out.println(all_delta_DAF.get(i));
		}
	}
	
//	@Override
//	public void logRunStats() {
//		
//		log.addLine("Out of " + tp_win.getSNPs().size() + " SNPs, " 
//				+ all_delta_DAF.size() + " were successful and " + unused_snps.size() 
//				+ " SNPs were unsuccessful");
//	}
	
	public List<Double> getDafStats() {
		return all_DAF;
	}
	
	public void printRStats() {
		
		double mean  = findMean(all_delta_DAF);
		double st_dev = findStandardDeviation(all_delta_DAF, mean);
		
		StringBuilder daf_sb = new StringBuilder();
		StringBuilder pos_sb = new StringBuilder();
		
		System.out.println("\nShowing R output: DAF");
		System.out.println("\tMean:\t" + mean);
		System.out.println("\tSt Dev:\t" + st_dev);
		
		for(int i = 0; i < all_delta_DAF.size(); i++) {
			
			daf_sb.append(all_delta_DAF.get(i) + ",");
			pos_sb.append(all_delta_DAF_snps.get(i).getPosition() + ",");
		}
		System.out.println("DAF =\t" + daf_sb.toString());
		System.out.println("Pos =\t" + pos_sb.toString());
	}

	
	private int getInstanceOfDerivedAllele(Individual[] indv, SNP core_snp, SNP anc_snp, int snp_index) {
		
		int count = 0;
		
		//When core_snps's a1 = derived type
		if(core_snp.getAllele0().equals(anc_snp.getAllele0())) {
			for(int i = 0; i < indv.length; i++) {
				
				//adding the index of the individual with the corresponding strand (1)
				int st1_allele = indv[i].getStrand1Allele(snp_index);
				if(st1_allele == 1)
					count++;
				
				//adding the index of the individual with the corresponding strand (2)
				int st2_allele = indv[i].getStrand2Allele(snp_index);
				if(st2_allele == 1)
					count++;
			}
		}
		//When core_snps's a0 = derived type
		else if(core_snp.getAllele1().equals(anc_snp.getAllele0())) {
			for(int i = 0; i < indv.length; i++) {
				
				int st1_allele = indv[i].getStrand1Allele(snp_index);
				if(st1_allele == 0)
					count++;
				
				int st2_allele = indv[i].getStrand2Allele(snp_index);
				if(st2_allele == 0)
					count++;
			}
		}
	
		return count;
	}
	
	private Window getEquivalentWindow(List<Window> cross_wins, Window target_win) {
		
		for(Window w : cross_wins) {
			if(w.getStPos() == target_win.getStPos()
					&& w.getEndPos() == target_win.getEndPos())
				return w;
		}
		
		return null;
	}
}
