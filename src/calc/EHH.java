package calc;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import tools.Individual;
import tools.SNP;
import tools.Window;
import tools.ExtendedHaplotype;

/*
 * Uses the equation:
 * 
 * 			Summation(i=1)=>s (eti _combine_ 2)
 * EHH	=	---------------
 * 				(Ct _combine_ 2)
 * 
 * 
 */

public class EHH {
	
	private static double SIGNIFICANT_EHH = 0.05;
	private static int MAX_DISTANCE = 3000000;
	
	private double cur_ehh_value;
	private double last_ehh_value;
	
	private LinkedList<ExtendedHaplotype> group;
	private LinkedList<Double> all_ehh_values;
	private LinkedList<Integer> all_ehh_positions;
	private List<Window> all_win;
	
	private SNP core_snp;//the core haplotype
	private SNP dwnstrm_snp;
	private SNP upstrm_snp;
	private SNP last_snp;
	private Window core_win;
	private Window dwnstrm_win;
	private Window upstrm_win;
	private Individual[] individuals;
	private ExtendedHaplotype all_haplo;
	
	private List<SNP> dups_test = new ArrayList<SNP>();
	
	public EHH(Window core_win, Individual[] individuals, SNP core_snp, 
			ExtendedHaplotype all_haplo, List<Window> all_win) {
		
		this.core_win = core_win;
		this.individuals = individuals;
		this.core_snp = core_snp;
		this.all_win = all_win;
		this.all_haplo = all_haplo;
		
		cur_ehh_value = 1.0; 
		last_ehh_value = 1.0;
		
		dwnstrm_snp = this.core_snp;
		upstrm_snp = this.core_snp;
		last_snp = this.core_snp;
		dwnstrm_win = this.core_win;
		upstrm_win = this.core_win;
		
		all_ehh_values = new LinkedList<Double>();
		all_ehh_positions = new LinkedList<Integer>();
		group = new LinkedList<ExtendedHaplotype>();
		group.add(all_haplo);
	}
	
	public boolean calcEhhToPosition(int end_pos) {
		
		dups_test = new ArrayList<SNP>();
		
		//This is the EHH eqn denominator
		int ct_comb_2 = combineSetBy2(all_haplo);
		
		all_ehh_values.add(1.0);
		all_ehh_positions.add(core_snp.getPosition());
		
		//Boundary position check for EHH calculation
		SNP nxt_snp = new SNP();
		while(isValidPosition(end_pos, last_snp.getPosition())) {
			
			nxt_snp = getClosestSNP(nxt_snp);
			
			//=========Duplicate Data Exception==========
			if(dups_test.contains(nxt_snp)) {
				System.out.println("\nWarning: CORE_" + core_snp + " has unhandled duplicate data");
				System.out.println("\t-Unexpected duplicate with " + nxt_snp);
				System.out.println("\t-Consider removing duplicate and rerun window");
				return false;
			}
			else
				dups_test.add(nxt_snp);
			//===========================================
			
			//both boundaries are hit; not enough variance for significance
			if(nxt_snp == null)
				return false;
			
			//a 3Mb max distance enforced (bigger bounds are unlikely to have significant LD)
			if(Math.abs(nxt_snp.getPosition() - core_snp.getPosition()) > MAX_DISTANCE)
				return false;
			
			//incorporates the new SNP into all extended haplotypes (increase length by 1)
			group = createNewExtHaploGroup(nxt_snp);
			
			cur_ehh_value = calcEHH(ct_comb_2);
			
			saveEHH(cur_ehh_value, nxt_snp);
		}
		
		return true;
	}
	
	/**
	 * Does not consider the core SNP for EHH analysis and grouping.
	 * Organizes all Extended Haplotypes until there are insignificant 
	 * homozygosity levels (EHH <= 0.05)
	 * 
	 * @return return true if the analysis generated significant results
	 * @throws EhhComputationException 
	 */
	public boolean calcSignificantEhhValues(double ehh_cutoff) {
		
		dups_test = new ArrayList<SNP>();
		
		//This is the EHH eqn denominator
		int ct_comb_2 = combineSetBy2(all_haplo); 
		
		all_ehh_values.add(1.0);
		all_ehh_positions.add(core_snp.getPosition());
		
		//Significance check of EHH value
		SNP nxt_snp = new SNP();
		while(cur_ehh_value > ehh_cutoff) { 
			
			nxt_snp = getClosestSNP(nxt_snp);
			
			//=========Duplicate Data Exception==========
			if(dups_test.contains(nxt_snp)) {
				System.out.println("*Warning: CORE_" + core_snp + " has unhandled duplicate data");
				System.out.println("\t-Unexpected duplicate with " + nxt_snp);
				System.out.println("\t-Consider removing duplicate and rerun window");
				return false;
			}
			else
				dups_test.add(nxt_snp);
			//===========================================
			
			//both boundaries are hit; not enough variance for significance
			if(nxt_snp == null)
				return false;
			
			//a 3Mb max distance enforced (bigger bounds are unlikely to have significant LD)
			if(Math.abs(nxt_snp.getPosition() - core_snp.getPosition()) > MAX_DISTANCE) 
				return false;
			
			//incorporates the new SNP into all extended haplotypes (increase length by 1)
			group = createNewExtHaploGroup(nxt_snp);
			
			cur_ehh_value = calcEHH(ct_comb_2);
			
			saveEHH(cur_ehh_value, nxt_snp);
		}
		
		return true;
	}
	
	public boolean calcSignificantEhhValues() {
		return calcSignificantEhhValues(SIGNIFICANT_EHH);
	}
	
	public double[] getEhhValues() {
		return convertAllEhhValuesToArray();
	}
	
	public int[] getEhhPositions() {
		return convertAllEhhPositionsToArray();
	}
	
	public double getLastEhhValue() {
		return last_ehh_value;
	}
	
	public SNP getLastSNP() {
		return last_snp;
	}
	
	private boolean isValidPosition(int boundary_pos, int prev_pos) {
		
		int boundary = Math.abs(core_snp.getPosition() - boundary_pos);
		int prev_range = Math.abs(core_snp.getPosition() - prev_pos);
		
		if(boundary >= prev_range)
			return true;
		
		return false;
	}
	
	private double[] convertAllEhhValuesToArray() {
		
		double[] all_ehh = new double[all_ehh_values.size()];
		for(int i = 0; i < all_ehh_values.size(); i++)
			all_ehh[i] = all_ehh_values.get(i);
		
		return all_ehh;
	}
	
	private int[] convertAllEhhPositionsToArray() {
		
		int[] all_pos = new int[all_ehh_positions.size()];
		for(int i = 0; i < all_ehh_positions.size(); i++) 
			all_pos[i] = all_ehh_positions.get(i);
		
		return all_pos;
	}
	
	private void saveEHH(double ehh, SNP nxt_snp) {
		
		if(nxt_snp.getPosition() >= core_snp.getPosition()) {
			all_ehh_values.addLast(ehh);
			all_ehh_positions.addLast(nxt_snp.getPosition());
		}
		else {
			all_ehh_values.addFirst(ehh);
			all_ehh_positions.addFirst(nxt_snp.getPosition());
		}
		
		last_ehh_value = ehh;
		last_snp = nxt_snp;
	}
	
	private double calcEHH(int ct_comb_2) {
		
		int sum = 0;
		for(ExtendedHaplotype eh : group) {
			sum += combineSetBy2(eh);
		}
		
		double ehh = (double) sum / (double) ct_comb_2;
		
		return ehh;
	}
	
	private int combineSetBy2(ExtendedHaplotype set) {
		
		int size = set.size();
		
		int score = 0;
		for(int i = 0; i < size; i++)
			score += size - (i + 1);
		
		return score;
	}
	
	private LinkedList<ExtendedHaplotype> createNewExtHaploGroup(SNP nxt_snp) {
		
		if(nxt_snp == null)
			return null;
		
		LinkedList<ExtendedHaplotype> updated_group = new LinkedList<ExtendedHaplotype>();
		
		int nxt_snp_index = -1;
		if(nxt_snp.getPosition() >= core_snp.getPosition())
			nxt_snp_index = upstrm_win.getSnpIndex(nxt_snp);
		else
			nxt_snp_index = dwnstrm_win.getSnpIndex(nxt_snp);
		
		for(ExtendedHaplotype eh : group) {
			
			//if the ExtHaplo is size = 1 it is by definition already completely unique 
			//and doen't need to be checked for more uniqueness
			if(eh.size() == 1){
				updated_group.add(eh);
			}
			else {
			
				ExtendedHaplotype eh_0 = new ExtendedHaplotype();
				ExtendedHaplotype eh_1 = new ExtendedHaplotype();
				
				while(!eh.isEmpty()) {
					
					int id = eh.getNextID();
					int strand = eh.getNextStrand();
					
					Individual indv = individuals[id];
					
					int allele = -1;
					if(strand == 1)
						allele = indv.getStrand1Allele(nxt_snp_index);
					if(strand == 2)
						allele = indv.getStrand2Allele(nxt_snp_index);
					
					if(allele == 0) 
						eh_0.add(id, strand);
					if(allele == 1) 
						eh_1.add(id, strand);
					
					eh.removeFirst();
				}
				
				if(eh_0.size() > 0)
					updated_group.add(eh_0);
				if(eh_1.size() > 0)
					updated_group.add(eh_1);
			}
		}
		
		return updated_group;
	}
	
	private SNP getClosestSNP(SNP prev_snp) {
		
		SNP temp_dwnstrm_snp = getNextDownstreamSNP(prev_snp);
		SNP temp_upstrm_snp = getNextUpstreamSNP(prev_snp);
		
		int dwnstrm_snp_length = -1;
		if(temp_dwnstrm_snp != null)
			dwnstrm_snp_length = Math.abs(core_snp.getPosition() - temp_dwnstrm_snp.getPosition());
		
		int upstrm_snp_length = -1;
		if(temp_upstrm_snp != null)
			upstrm_snp_length = Math.abs(core_snp.getPosition() - temp_upstrm_snp.getPosition());
			
		if(upstrm_snp_length == -1 && dwnstrm_snp_length == -1) {
			return null;
		} else if(upstrm_snp_length == -1) {
			return incrementDownstream(temp_dwnstrm_snp);
		} else if(dwnstrm_snp_length == -1) {
			return incrementUpstream(temp_upstrm_snp);
		} else if(dwnstrm_snp_length < upstrm_snp_length) {
			return incrementDownstream(temp_dwnstrm_snp);
		} else if(upstrm_snp_length <= dwnstrm_snp_length){
			return incrementUpstream(temp_upstrm_snp);
		}
		
		return null;
	}
	
	private SNP incrementUpstream(SNP temp_upstrm_snp) {
		
		int old_index = upstrm_win.getSnpIndex(upstrm_snp);// a check for incrementing the window
		
		upstrm_snp = temp_upstrm_snp;
		
		//this is done only if the new upstrm_snp is outside of the current window
		if(!upstrm_win.containsIndex(old_index + 1))
			upstrm_win = findWindow(old_index + 1);
		
		return upstrm_snp;
	}
	
	private SNP incrementDownstream(SNP temp_dwnstrm_snp) {
		
		int old_index = dwnstrm_win.getSnpIndex(dwnstrm_snp);// a check for incrementing the window
		
		dwnstrm_snp = temp_dwnstrm_snp;
		
		//this is done only if the new dwnstrm_snp is outside of the current window
		if(!dwnstrm_win.containsIndex(old_index - 1))
			dwnstrm_win = findWindow(old_index - 1);
		
		return dwnstrm_snp;
	}
	
	private SNP getNextDownstreamSNP(SNP prev_snp) {
		
		//to skip this function if you have reached the boundary of the chr
		if(dwnstrm_snp == null || dwnstrm_win == null)
			return null;
		
		SNP nxt_dwn_snp = new SNP();
		int nxt_index = dwnstrm_win.getSnpIndex(dwnstrm_snp) - 1;
		
		if(dwnstrm_win.containsIndex(nxt_index))
			nxt_dwn_snp = dwnstrm_win.getSNP(nxt_index);
		else {
			
			Window temp_dwnstrm_win = findWindow(nxt_index);
			
			if(temp_dwnstrm_win == null) {
				dwnstrm_win = null;
				dwnstrm_snp = null;
				return null;
			}
			
			nxt_dwn_snp = temp_dwnstrm_win.getSNP(nxt_index);
		}
		
		return nxt_dwn_snp;
	}
	
	private SNP getNextUpstreamSNP(SNP prev_snp) {
		
		//to skip this function if you have reached the boundary of the chr
		if(upstrm_snp == null || upstrm_win == null)
			return null;
		
		SNP nxt_up_snp = new SNP();
		int nxt_index = upstrm_win.getSnpIndex(upstrm_snp) + 1;
		
		if(upstrm_win.containsIndex(nxt_index))
			nxt_up_snp = upstrm_win.getSNP(nxt_index);
		else {
			
			Window temp_upstrm_win = findWindow(nxt_index);
			
			if(temp_upstrm_win == null) {
				upstrm_win = null;
				upstrm_snp = null;
				return null;
			}
			
			nxt_up_snp = temp_upstrm_win.getSNP(nxt_index);
		}
		
		return nxt_up_snp;
	}
	
	private Window findWindow(int index) {
		
		for(Window w : all_win) {
			if(w.containsIndex(index)) 
				return w;
		}
		
		return null;
	}
}
