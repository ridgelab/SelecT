package envi;

import java.util.ArrayList;
import java.util.List;

import tools.Individual;
import tools.SNP;
import tools.Window;

public class PopIntersector {
	
	//intersection of target and cross populations (txin)
	private List<Window> txin_wins;//for listing the windows in the intersection
	private Individual[] tp_inx_indv;//for listing the individual alleles in the intersection
	private Individual[] xp_int_indv;
	
	//intersection of cross population and out-group population (xoin)
	private List<Window> xoin_wins;
	private Individual[] xp_ino_indv;
	private Individual[] op_inx_indv;
	
	public PopIntersector() {
		
		txin_wins = null;
		xoin_wins = null;
		
		tp_inx_indv = null;
		xp_int_indv = null;
		xp_ino_indv = null;
		op_inx_indv = null;
	}

	public void intersectCrossWithTargetPopulations(List<Window> tp_wins, 
													List<Window> xp_wins, 
													Individual[] tp_indv,
													Individual[] xp_indv) {
		
		Individual[] xp_indv_insect = new Individual[xp_indv.length];
		Individual[] tp_indv_insect = new Individual[tp_indv.length];
		
		for(int i = 0; i < xp_indv_insect.length; i++)
			xp_indv_insect[i] = new Individual(xp_indv[i].getID(), xp_indv[i].getChr());
		for(int i = 0; i < tp_indv_insect.length; i++)
			tp_indv_insect[i] = new Individual(tp_indv[i].getID(), tp_indv[i].getChr());
		
		List<Window> wins_insect = new ArrayList<Window>();
		
		compareWindows(wins_insect, xp_wins, xp_indv, xp_indv_insect, tp_wins, tp_indv, tp_indv_insect);

		//set the global variables
		txin_wins = wins_insect;
		xp_int_indv = xp_indv_insect;
		tp_inx_indv = tp_indv_insect;
	}
	
	public void intersectCrossWithOutgroupPopulations(List<Window> op_wins, 
														List<Window> xp_wins, 
														Individual[] op_indv,
														Individual[] xp_indv) {
		
		Individual[] xp_indv_insect = new Individual[xp_indv.length];
		Individual[] op_indv_insect = new Individual[op_indv.length];
		
		for(int i = 0; i < xp_indv_insect.length; i++)
			xp_indv_insect[i] = new Individual(xp_indv[i].getID(), xp_indv[i].getChr());
		for(int i = 0; i < op_indv_insect.length; i++)
			op_indv_insect[i] = new Individual(op_indv[i].getID(), op_indv[i].getChr());
		
		List<Window> wins_insect = new ArrayList<Window>();
		
		compareWindows(wins_insect, xp_wins, xp_indv, xp_indv_insect, op_wins, op_indv, op_indv_insect);

		//set the global variables
		xoin_wins = wins_insect;
		xp_ino_indv = xp_indv_insect;
		op_inx_indv = op_indv_insect;
	}
	
	public List<Window> getTargetXCrossWins() {
		return txin_wins;
	}

	public Individual[] getTargetXCrossIndv() {
		return tp_inx_indv;
	}

	public Individual[] getCrossXTargetIndv() {
		return xp_int_indv;
	}

	public List<Window> getCrossXOutWins() {
		return xoin_wins;
	}

	public Individual[] getCrossXOutIndv() {
		return xp_ino_indv;
	}

	public Individual[] getOutXCrossIndv() {
		return op_inx_indv;
	}

	private void compareWindows(List<Window> wins_insect,
								List<Window> p1_wins,
								Individual[] p1_indv,
								Individual[] p1_indv_insect,
								List<Window> p2_wins,
								Individual[] p2_indv,
								Individual[] p2_indv_insect) {

		wins_insect.add(new Window(0, 0, 0));//to initialize the comparator list

		for(int i = 0; i < p1_wins.size(); i++) {
			for(int j = 0; j < p2_wins.size(); j++) {
				Window p1_win = p1_wins.get(i);
				Window p2_win = p2_wins.get(j);

				int p1_win_st = p1_win.getStPos();
				int p1_win_end = p1_win.getEndPos();
				int p2_win_st = p2_win.getStPos();
				int p2_win_end = p2_win.getEndPos();

				if(p1_win_st == p2_win_st
						&& p1_win_end == p2_win_end) {			

					List<SNP> p1_win_snps = p1_wins.get(i).getSNPs();
					List<SNP> p2_win_snps = p2_wins.get(j).getSNPs();

					compareSNPs(p1_win_snps, 
								p2_win_snps, 
								wins_insect, 
								p1_win_st, 
								p2_win_end, 
								p1_win, 
								p2_win, 
								p1_indv,
								p2_indv,
								p1_indv_insect, 
								p2_indv_insect);	
				}
			}	
		}

		wins_insect.remove(0);//to get rid of the initial window
	}
	
	private void compareSNPs(List<SNP> p1_win_snps, 
								List<SNP> p2_win_snps,
								List<Window> wins_insect,
								int p1_win_st,
								int p1_win_end,
								Window p1_win,
								Window p2_win,
								Individual[] p1_indv,
								Individual[] p2_indv,
								Individual[] p1_indv_insect,
								Individual[] p2_indv_insect) {

		for(int k = 0; k < p1_win_snps.size(); k++) {
			for(int l = 0; l < p2_win_snps.size(); l++) {
				if(p1_win_snps.get(k).sameAs(p2_win_snps.get(l))) { 
					if(!containsWindow(wins_insect, p1_win_st, p1_win_end)) {
	
						//make and put a new window window in wins_insect
						Window last_win = wins_insect.get(wins_insect.size() - 1);
						int last_win_indx = last_win.getStIndex() + last_win.getSnpListSize() - 1;
						last_win.setEndIndex(last_win_indx);
	
						Window new_win = new Window(p1_win_st, p1_win_end, (last_win_indx + 1));
						wins_insect.add(new_win);
					}

					Window cur_win = getCurWindow(wins_insect, p1_win_st, p1_win_end);
					int cur_win_indx = wins_insect.indexOf(cur_win);

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
											p1_indv_insect, 
											p2_indv_insect);

					cur_win.addSNP(p1_snp);
					cur_win.setEndIndex(cur_win.getStIndex() + cur_win.getSnpListSize() - 1);

					wins_insect.set(cur_win_indx, cur_win);

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
											Individual[] p1_indv_insect,
											Individual[] p2_indv_insect) {

		//Adding alleles to p1 population's individuals
		for(int m = 0; m < p1_indv_insect.length; m++) {
			Integer str_1 = p1_indv[m].getStrand1Allele(p1_indx);
			Integer str_2 = p1_indv[m].getStrand2Allele(p1_indx);

			p1_indv_insect[m].addAlleleToStrand1(str_1.toString());
			p1_indv_insect[m].addAlleleToStrand2(str_2.toString());
		}

		//Adding alleles to p2 population's individuals
		for(int i = 0; i < p2_indv_insect.length; i++) {

			Integer str_1 = p2_indv[i].getStrand1Allele(p2_indx);
			Integer str_2 = p2_indv[i].getStrand2Allele(p2_indx);

			//switch allele types because they are reported on opposite a0 or a1 column
			if(p1_snp.getAllele0().equals(p2_snp.getAllele1())) {

				if(str_1 == 0)
					str_1 = 1;
				else
					str_1 = 0;

				if(str_2 == 0)
					str_2 = 1;
				else
					str_2 = 0;
			}

			p2_indv_insect[i].addAlleleToStrand1(str_1.toString());
			p2_indv_insect[i].addAlleleToStrand2(str_2.toString());
		}
	}
	
	private Window getCurWindow(List<Window> wins, int st, int end) {
		
		for(Window w : wins) {
			if(w.getStPos() == st && w.getEndPos() == end) 
				return w;
		}
		
		return null;
	}
	
	private boolean containsWindow(List<Window> wins, int st, int end) {
		
		
		for(Window w : wins) {
			if(w.getStPos() == st && w.getEndPos() == end)
				return true;
		}
		
		return false;
	}
	
}
