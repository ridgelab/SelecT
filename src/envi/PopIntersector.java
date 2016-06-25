package envi;

import java.util.ArrayList;
import java.util.List;

import tools.Individual;
import tools.SNP;
import tools.Window;

/**
 * Performs intersections of the populations to find overlapping positions between all three populations.
 * Populations are compared by position and not variantID except in cases of irregular data types.
 * The overlapping is used for deltaDAF, XPEHH, and Fst, which compare target populations to cross
 * populations and therefore require the overlap locations.
 */
public class PopIntersector {
	
	//intersection of target and cross populations (txin)
	private List<Window> txin_wins;//for listing the windows in the intersection
	private Individual[] tp_inx_indv;//for listing the individual alleles in the intersection
	private Individual[] xp_int_indv;
	
	//intersection of cross population and out-group population (xoin)
	private List<Window> xoin_wins;
	private Individual[] xp_ino_indv;
	private Individual[] op_inx_indv;
	
	/**
	 * Simple constructor for PopIntersector class
	 */
	public PopIntersector() {
		
		txin_wins = null;
		xoin_wins = null;
		
		tp_inx_indv = null;
		xp_int_indv = null;
		xp_ino_indv = null;
		op_inx_indv = null;
	}

	/**
	 * Performs intersection between the target population and the cross population.
	 * Inputs are lists of Windows for the target and cross populations, 
	 * as well as arrays of individuals. The function creates new Individual 
	 * arrays from the inputs for the intersects and calls functions to compare the 
	 * populations to each other
	 * 
	 * @param tp_wins	windows for target population
	 * @param xp_wins	windows for cross population
	 * @param tp_indv	individuals for target population
	 * @param xp_indv	individuals for cross population
	 */
	public void intersectCrossWithTargetPopulations(List<Window> tp_wins, 
												List<Window> xp_wins,
												Individual[] tp_indv, 
												Individual[] xp_indv) {
		
		Individual[] xp_indv_intersect = new Individual[xp_indv.length];
		Individual[] tp_indv_intersect = new Individual[tp_indv.length];
		
		//TODO: isn't it inefficient to copy the individuals like this?
		for (int i = 0; i < xp_indv_intersect.length; i++) {
			xp_indv_intersect[i] = new Individual(xp_indv[i].getID(), xp_indv[i].getChr());
		}
		for (int i = 0; i < tp_indv_intersect.length; i++){
			tp_indv_intersect[i] = new Individual(tp_indv[i].getID(), tp_indv[i].getChr());
		}
		
		List<Window> wins_intersect = new ArrayList<Window>();
		
		compareWindows(wins_intersect, xp_wins, xp_indv, xp_indv_intersect, tp_wins, tp_indv, tp_indv_intersect);

		//set the global variables
		txin_wins = wins_intersect;
		xp_int_indv = xp_indv_intersect;
		tp_inx_indv = tp_indv_intersect;
	}
	
	/**
	 * Performs intersection between the cross population and the outgroup population.
	 * Inputs are lists of Windows for the cross and outgroup populations, 
	 * as well as arrays of individuals. The function creates new Individual 
	 * arrays from the inputs for the intersects and calls functions to compare the 
	 * populations to each other
	 * 
	 * @param op_wins	windows for outgroup population
	 * @param xp_wins	windows for cross population
	 * @param op_indv	individuals for outgroup population
	 * @param xp_indv	individuals for cross population
	 */
	public void intersectCrossWithOutgroupPopulations(List<Window> op_wins,
													List<Window> xp_wins, 
													Individual[] op_indv,	
													Individual[] xp_indv) {
		
		Individual[] xp_indv_intersect = new Individual[xp_indv.length];
		Individual[] op_indv_intersect = new Individual[op_indv.length];
		
		for (int i = 0; i < xp_indv_intersect.length; i++) {
			xp_indv_intersect[i] = new Individual(xp_indv[i].getID(), xp_indv[i].getChr());
		}
		for (int i = 0; i < op_indv_intersect.length; i++) {
			op_indv_intersect[i] = new Individual(op_indv[i].getID(), op_indv[i].getChr());
		}
		
		List<Window> wins_intersect = new ArrayList<Window>();
		
		compareWindows(wins_intersect, xp_wins, xp_indv, xp_indv_intersect, op_wins, op_indv, op_indv_intersect);

		//set the global variables
		xoin_wins = wins_intersect;
		xp_ino_indv = xp_indv_intersect;
		op_inx_indv = op_indv_intersect;
	}
	
	/**
	 * Gets the target population-by-cross population intersection windows
	 */
	public List<Window> getTargetXCrossWins() {
		return txin_wins;
	}

	/**
	 * Gets the target population-by-cross population intersection individuals
	 */
	public Individual[] getTargetXCrossIndv() {
		return tp_inx_indv;
	}

	/**
	 * Gets the cross population-by-target population intersection individuals
	 */
	public Individual[] getCrossXTargetIndv() {
		return xp_int_indv;
	}

	/**
	 * Gets the cross population-by-outgroup population intersection windows
	 */
	public List<Window> getCrossXOutWins() {
		return xoin_wins;
	}

	/**
	 * Gets the cross population-by-outgroup population intersection individuals
	 */
	public Individual[] getCrossXOutIndv() {
		return xp_ino_indv;
	}

	/**
	 * Gets the outgroup population-by-cross population intersection individuals
	 */
	public Individual[] getOutXCrossIndv() {
		return op_inx_indv;
	}

	private void compareWindows(List<Window> wins_intersect, 
								List<Window> p1_wins,	
								Individual[] p1_indv,
								Individual[] p1_indv_intersect, 
								List<Window> p2_wins, 
								Individual[] p2_indv, 
								Individual[] p2_indv_intersect) {

		wins_intersect.add(new Window(0, 0, 0));//to initialize the comparator list

		for (int i = 0; i < p1_wins.size(); i++) {
			for (int j = 0; j < p2_wins.size(); j++) {
				Window p1_win = p1_wins.get(i);
				Window p2_win = p2_wins.get(j);

				int p1_win_st = p1_win.getStPos();
				int p1_win_end = p1_win.getEndPos();
				int p2_win_st = p2_win.getStPos();
				int p2_win_end = p2_win.getEndPos();

				if (p1_win_st == p2_win_st
						&& p1_win_end == p2_win_end) {			

					List<SNP> p1_win_snps = p1_wins.get(i).getSNPs();
					List<SNP> p2_win_snps = p2_wins.get(j).getSNPs();

					compareSNPs(p1_win_snps, 
								p2_win_snps, 
								wins_intersect, 
								p1_win_st, 
								p2_win_end, 
								p1_win, 
								p2_win, 
								p1_indv,
								p2_indv,
								p1_indv_intersect, 
								p2_indv_intersect);	
				}
			}	
		}

		wins_intersect.remove(0);//to get rid of the initial window
	}
	
	private void compareSNPs(List<SNP> p1_win_snps, 
								List<SNP> p2_win_snps,
								List<Window> wins_intersect,
								int p1_win_st,
								int p1_win_end,
								Window p1_win,
								Window p2_win,
								Individual[] p1_indv,
								Individual[] p2_indv,
								Individual[] p1_indv_intersect,
								Individual[] p2_indv_intersect) {

		for (int k = 0; k < p1_win_snps.size(); k++) {
			for (int l = 0; l < p2_win_snps.size(); l++) {
				if (p1_win_snps.get(k).sameAs(p2_win_snps.get(l))) { 
					if (!containsWindow(wins_intersect, p1_win_st, p1_win_end)) {
	
						//make and put a new window window in wins_intersect
						Window last_win = wins_intersect.get(wins_intersect.size() - 1);
						int last_win_indx = last_win.getStIndex() + last_win.getSnpListSize() - 1;
						last_win.setEndIndex(last_win_indx);
	
						Window new_win = new Window(p1_win_st, p1_win_end, (last_win_indx + 1));
						wins_intersect.add(new_win);
					}

					Window cur_win = getCurWindow(wins_intersect, p1_win_st, p1_win_end);
					int cur_win_indx = wins_intersect.indexOf(cur_win);

					SNP p1_snp = p1_win_snps.get(k);
					SNP p2_snp = p2_win_snps.get(l);

					int p1_indx = p1_win.getSnpIndex(p1_snp);
					int p2_indx = p2_win.getSnpIndex(p2_snp);

					addAllelesToIndividuals(p1_indx, 
											p2_indx, 
											p1_snp, 
											p2_snp, 
											p1_indv,
											p2_indv,
											p1_indv_intersect, 
											p2_indv_intersect);

					cur_win.addSNP(p1_snp);
					cur_win.setEndIndex(cur_win.getStIndex() + cur_win.getSnpListSize() - 1);

					wins_intersect.set(cur_win_indx, cur_win);

				}
			}	
		}
	}
	
	private void addAllelesToIndividuals(int p1_indx, 
											int p2_indx, 
											SNP p1_snp, 
											SNP p2_snp,
											Individual[] p1_indv,
											Individual[] p2_indv,
											Individual[] p1_indv_intersect,
											Individual[] p2_indv_intersect) {

		//Adding alleles to p1 population's individuals
		for (int m = 0; m < p1_indv_intersect.length; m++) {
			//TODO: Causing death
			Integer str_1 = p1_indv[m].getAlleleFromStrand(p1_indx, true);
			Integer str_2 = p1_indv[m].getAlleleFromStrand(p1_indx, false);

			p1_indv_intersect[m].addAlleleToStrand(str_1.toString().charAt(0), true);
			p1_indv_intersect[m].addAlleleToStrand(str_2.toString().charAt(0), false);

		}

		//Adding alleles to p2 population's individuals
		for (int i = 0; i < p2_indv_intersect.length; i++) {

			Integer str_1 = p2_indv[i].getAlleleFromStrand(p2_indx, true);
			Integer str_2 = p2_indv[i].getAlleleFromStrand(p2_indx, false);
			
			//switch allele types because they are reported on opposite a0 or a1 column
			if (p1_snp.getAllele0().equals(p2_snp.getAllele1())) {

				if (str_1 == 0) {
					str_1 = 1;
				}
				else {
					str_1 = 0;
				}

				if (str_2 == 0){
					str_2 = 1;
				}
				else {
					str_2 = 0;
				}
			}

			p2_indv_intersect[i].addAlleleToStrand(str_1.toString().charAt(0), true);
			p2_indv_intersect[i].addAlleleToStrand(str_2.toString().charAt(0), false);
		}
	}
	
	private Window getCurWindow(List<Window> wins, int st, int end) {
		
		for (Window w : wins) {
			if(w.getStPos() == st && w.getEndPos() == end) {
				return w;
			}
		}
		
		return null;
	}
	
	private boolean containsWindow(List<Window> wins, int st, int end) {
		
		
		for (Window w : wins) {
			if (w.getStPos() == st && w.getEndPos() == end) {
				return true;
			}
		}
		
		return false;
	}
	
}

