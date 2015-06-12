package tools;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeMap;

import calc.HaplotypeTests;

public class WindowStats implements Comparable<WindowStats>{
	
	private int st_pos;
	private int end_pos;
	
	private List<SNP> ihs_snps;
	private List<SNP> xpehh_snps;
	private List<SNP> ihh_snps;
	private List<SNP> ddaf_snps;
	private List<SNP> daf_snps;
	private List<SNP> fst_snps;
	//private List<SNP> tajd_snps;
	//private List<SNP>	new_snps;
	
	private List<Double> ihs_stats;
	private List<Double> xpehh_stats;
	private List<Double> ihh_stats;
	private List<Double> ddaf_stats;
	private List<Double> daf_stats;
	private List<Double> fst_stats;
	//private List<Double> tajd_stats;
	//private List<Double> new_stats;
	
	private TreeMap<SNP, Double> win_scores_prod_unstd;
	private TreeMap<SNP, Double> win_scores_mean_unstd;
	private TreeMap<SNP, Double> win_scores_prod_std;
	private TreeMap<SNP, Double> win_scores_mean_std;
	
	public WindowStats() {
		
		this(-1, -1);
	}
	
	public WindowStats(int st_pos, int end_pos) {
		
		this.st_pos = st_pos;
		this.end_pos = end_pos;
		
		ihs_snps = new ArrayList<SNP>();
		xpehh_snps = new ArrayList<SNP>();
		ihh_snps = new ArrayList<SNP>();
		ddaf_snps = new ArrayList<SNP>();
		daf_snps = new ArrayList<SNP>();
		fst_snps = new ArrayList<SNP>();
		//tajd_snps = new ArrayList<SNP>();
		//new_snps = new ArrayList<SNP>();
		
		ihs_stats = new ArrayList<Double>();
		xpehh_stats = new ArrayList<Double>();
		ihh_stats = new ArrayList<Double>();
		ddaf_stats = new ArrayList<Double>();
		daf_stats = new ArrayList<Double>();
		fst_stats = new ArrayList<Double>();
		//tajd_stats = new ArrayList<Double>();
		//new_stats = new ArrayList<Double>();
		
		win_scores_prod_unstd = new TreeMap<SNP, Double>();
		win_scores_mean_unstd = new TreeMap<SNP, Double>();
		win_scores_prod_std = new TreeMap<SNP, Double>();
		win_scores_mean_std = new TreeMap<SNP, Double>();
	}
	
	public List<SNP> getAllSNPs() {
		
		List<SNP> all_snps = new LinkedList<SNP>();
		
		all_snps = buildAllSNPs(all_snps, ihs_snps);
		all_snps = buildAllSNPs(all_snps, xpehh_snps);
		all_snps = buildAllSNPs(all_snps, ihh_snps);
		all_snps = buildAllSNPs(all_snps, ddaf_snps);
		all_snps = buildAllSNPs(all_snps, daf_snps);
		all_snps = buildAllSNPs(all_snps, fst_snps);
		//all_snps = buildAllSNPs(all_snps, tajd_snps);
		//all_snps = buildAllSNPs(all_snps, new_snps);
		
		Collections.sort(all_snps);
			
		return all_snps;
	}
	
	public int getTotNumSNPs() {
		
		List<SNP> all_snps = getAllSNPs();
		
		return all_snps.size();
	}
	
	public int getNextPosition(int prev_pos) {
		
		int nxt_pos = end_pos;
		
		nxt_pos = comparePositions(nxt_pos, prev_pos, ihs_snps);
		nxt_pos = comparePositions(nxt_pos, prev_pos, xpehh_snps);
		nxt_pos = comparePositions(nxt_pos, prev_pos, ihh_snps);
		nxt_pos = comparePositions(nxt_pos, prev_pos, ddaf_snps);
		nxt_pos = comparePositions(nxt_pos, prev_pos, daf_snps);
		nxt_pos = comparePositions(nxt_pos, prev_pos, fst_snps);
		//nxt_pos = comparePositions(nxt_pos, prev_pos, tajd_snps);
		//nxt_pos = comparePositions(nxt_pos, prev_pos, new_snps);
		
		if(nxt_pos == prev_pos)
			return -1;
		
		return nxt_pos;
	}
	
	public int getStPos() {
		return st_pos;
	}
	
	public int getEndPos() {
		return end_pos;
	}

	public List<SNP> getIHSsnps() {
		return ihs_snps;
	}

	public List<Double> getIHSstats() {
		return ihs_stats;
	}

	public void setIHS(List<Double> ihs_stats, List<SNP> ihs_snps) {
		this.ihs_stats = ihs_stats;
		this.ihs_snps = ihs_snps;
	}
	
	public void addIHS(List<Double> ihs_stats, List<SNP> ihs_snps) {
		this.ihs_stats.addAll(ihs_stats);
		this.ihs_snps.addAll(ihs_snps);
	}

	public List<SNP> getXPEHHsnps() {
		return xpehh_snps;
	}

	public List<Double> getXPEHHstats() {
		return xpehh_stats;
	}

	public void setXPEHH(List<Double> xpehh_stats, List<SNP> xpehh_snps) {
		this.xpehh_stats = xpehh_stats;
		this.xpehh_snps = xpehh_snps;
	}
	
	public void addXPEHH(List<Double> xpehh_stats, List<SNP> xpehh_snps) {
		this.xpehh_stats.addAll(xpehh_stats);
		this.xpehh_snps.addAll(xpehh_snps);
	}

	public List<SNP> getIHHsnps() {
		return ihh_snps;
	}
	
	public List<Double> getIHHstats() {
		return ihh_stats;
	}

	public void setIHH(List<Double> ihh_stats, List<SNP> ihh_snps) {
		this.ihh_stats = ihh_stats;
		this.ihh_snps = ihh_snps;
	}
	
	public void addIHH(List<Double> ihh_stats, List<SNP> ihh_snps) {
		this.ihh_stats.addAll(ihh_stats);
		this.ihh_snps.addAll(ihh_snps);
	}

	public List<SNP> getDDAFsnps() {
		return ddaf_snps;
	}

	public List<Double> getDDAFstats() {
		return ddaf_stats;
	}

	public void setDDAF(List<Double> daf_stats, List<SNP> daf_snps) {
		this.ddaf_stats = daf_stats;
		this.ddaf_snps = daf_snps;
	}
	
	public void addDDAF(List<Double> ddaf_stats, List<SNP> ddaf_snps) {
		this.ddaf_stats.addAll(ddaf_stats);
		this.ddaf_snps.addAll(ddaf_snps);
	}
	
	public List<SNP> getDAFsnps() {
		return daf_snps;
	}

	public List<Double> getDAFstats() {
		return daf_stats;
	}

	public void setDAF(List<Double> daf_stats, List<SNP> daf_snps) {
		this.daf_stats = daf_stats;
		this.daf_snps = daf_snps;
	}
	
	public void addDAF(List<Double> daf_stats, List<SNP> daf_snps) {
		this.daf_stats.addAll(daf_stats);
		this.daf_snps.addAll(daf_snps);
	}
	
//	public List<SNP> getTAJDsnps() {
//		return daf_snps;
//	}
//
//	public List<Double> getTAJDstats() {
//		return daf_stats;
//	}
//
//	public void setTAJD(List<Double> tajd_stats, List<SNP> tajd_snps) {
//		this.tajd_stats = tajd_stats;
//		this.tajd_snps = tajd_snps;
//	}
//	
//	public void addTAFD(List<Double> tajd_stats, List<SNP> tajd_snps) {
//		this.tajd_stats.addAll(tajd_stats);
//		this.tajd_snps.addAll(tajd_snps);
//	}
	
//	public List<SNP> getNEWsnps() {
//		return new_snps;
//	}
//
//	public List<Double> getNEWstats() {
//		return new_stats;
//	}
//
//	public void setNEW(List<Double> new_stats, List<SNP> new_snps) {
//		this.new_stats = new_stats;
//		this.new_snps = new_snps;
//	}
//	
//	public void addNEW(List<Double> new_stats, List<SNP> new_snps) {
//		this.new_stats.addAll(new_stats);
//		this.new_snps.addAll(new_snps);
//	}

	public List<SNP> getFSTsnps() {
		return fst_snps;
	}
	
	public List<Double> getFSTstats() {
		return fst_stats;
	}

	public void setFST(List<Double> fst_stats, List<SNP> fst_snps) {
		this.fst_stats = fst_stats;
		this.fst_snps = fst_snps;
	}
	
	public void addFST(List<Double> fst_stats, List<SNP> fst_snps) {
		this.fst_stats.addAll(fst_stats);
		this.fst_snps.addAll(fst_snps);
	}
	
	public Double getScore(List<SNP> snps, List<Double> stats, SNP snp) {
		
		for(int i = 0; i < snps.size(); i++) {
			if(snps.get(i).sameAs(snp))
				return stats.get(i);
		}
		
		return Double.NaN;
	}
	
	public Double getIhsScore(SNP snp) {
		
		for(int i = 0; i < ihs_snps.size(); i++) {
			if(ihs_snps.get(i).sameAs(snp))
				return ihs_stats.get(i);
		}
		
		return Double.NaN;
	}
	
	public Double getIhhScore(SNP snp) {
		
		for(int i = 0; i < ihh_snps.size(); i++) {
			if(ihh_snps.get(i).sameAs(snp))
				return ihh_stats.get(i);
		}
		
		return Double.NaN;
	}
	
	public Double getXpehhScore(SNP snp) {
		
		for(int i = 0; i < xpehh_snps.size(); i++) {
			if(xpehh_snps.get(i).sameAs(snp))
				return xpehh_stats.get(i);
		}
		
		return Double.NaN;
	}
	
	public Double getDDafScore(SNP snp) {
		
		for(int i = 0; i < ddaf_snps.size(); i++) {
			if(ddaf_snps.get(i).sameAs(snp))
				return ddaf_stats.get(i);
		}
		
		return Double.NaN;
	}
	
	public Double getDafScore(SNP snp) {
		
		for(int i = 0; i < daf_snps.size(); i++) {
			if(daf_snps.get(i).sameAs(snp))
				return daf_stats.get(i);
		}
		
		return Double.NaN;
	}

//	public Double getTajDScore(SNP snp) {
//	
//		for(int i = 0; i < tajd_snps.size(); i++) {
//			if(tajd_snps.get(i).sameAs(snp))
//				return tajd_stats.get(i);
//		}
//	
//		return Double.NaN;
//	}
	
//	public Double getNewScore(SNP snp) {
//	
//		for(int i = 0; i < new_snps.size(); i++) {
//			if(new_snps.get(i).sameAs(snp))
//				return new_stats.get(i);
//		}
//	
//		return Double.NaN;
//	}
	
	public Double getFstScore(SNP snp) {
		
		for(int i = 0; i < fst_snps.size(); i++) {
			if(fst_snps.get(i).sameAs(snp))
				return fst_stats.get(i);
		}
		
		return Double.NaN;
	}
	
	public boolean containsSNP(SNP s) {
		if(st_pos <= s.getPosition() && end_pos >= s.getPosition())
			return true;
		else
			return false;
	}
	
	public void addUnstdPopScore(SNP s, Double score) {
		
		if(!score.equals(Double.NaN))
			win_scores_prod_unstd.put(s, score);
	}
	
	public void addUnstdMopScore(SNP s, Double score) {
		
		if(!score.equals(Double.NaN))
			win_scores_mean_unstd.put(s, score);
	}
	
	public void addUnstdPoP(TreeMap<SNP, Double> unstd_pop) {
		win_scores_prod_unstd.putAll(unstd_pop);
	}
	
	public void addUnstdMoP(TreeMap<SNP, Double> unstd_mop) {
		win_scores_mean_unstd.putAll(unstd_mop);
	}
	
	public void addStdPopScore(SNP s, Double score) {
		
		if(!score.equals(Double.NaN))
			win_scores_prod_std.put(s, score);
	}
	
	public void addStdMopScore(SNP s, Double score) {
		
		if(!score.equals(Double.NaN))
			win_scores_mean_std.put(s, score);
	}
	
	public void addStdPoP(TreeMap<SNP, Double> std_pop) {
		win_scores_prod_std.putAll(std_pop);
	}
	
	public void addStdMoP(TreeMap<SNP, Double> std_mop) {
		win_scores_mean_std.putAll(std_mop);
	}
	
	public TreeMap<SNP, Double> getUnstdPoP() {
		return win_scores_prod_unstd;
	}
	
	public TreeMap<SNP, Double> getUnstdMoP() {
		return win_scores_mean_unstd;
	}
	
	public Double getUnstdPopScore(SNP s) {
		
		Double score = win_scores_prod_unstd.get(s);
		if(score == null)
			return Double.NaN;
		else
			return score;
	}
	
	public TreeMap<SNP, Double> getStdPoP() {
		return win_scores_prod_std;
	}
	
	public TreeMap<SNP, Double> getStdMoP() {
		return win_scores_mean_std;
	}
	
	public Double getUnstdMopScore(SNP s) {
		
		Double score = win_scores_mean_unstd.get(s);
		if(score == null)
			return Double.NaN;
		else
			return score;
	}
	
	public Double getStdPopScore(SNP s) {
		
		Double score = win_scores_prod_std.get(s);
		if(score == null)
			return Double.NaN;
		else
			return score;
	}
	
	public Double getStdMopScore(SNP s) {
		
		Double score = win_scores_mean_std.get(s);
		if(score == null)
			return Double.NaN;
		else
			return score;
	}
	
	public void normalizeUnstdCompositeScores() {
		
		win_scores_prod_std = normalizeData(win_scores_prod_unstd);
		win_scores_mean_std = normalizeData(win_scores_mean_unstd);
	}
	
	private TreeMap<SNP, Double> normalizeData(TreeMap<SNP, Double> unstd_cms) {
		
		List<Double> all_values = new LinkedList<Double>();
		TreeMap<SNP, Double> std_cms = new TreeMap<SNP, Double>();
		
		Iterator<SNP> itr = unstd_cms.navigableKeySet().iterator();
		while(itr.hasNext()) 
			all_values.add(unstd_cms.get(itr.next()));
		
		all_values = HaplotypeTests.normalizeData(all_values);
		
		itr = unstd_cms.navigableKeySet().iterator();
		int indx = 0;
		while(itr.hasNext()) {
			std_cms.put(itr.next(), all_values.get(indx));
			indx++;
		}
		
		return std_cms;
	}
	
	private List<SNP> buildAllSNPs(List<SNP> all_snps, List<SNP> snps) {
		
		for(int i = 0; i < snps.size(); i++) {
			if(!containsSNP(all_snps, snps.get(i)))
					all_snps.add(snps.get(i));
		}
		
		return all_snps;
	}
	
	private boolean containsSNP(List<SNP> all_snps, SNP snp) {
		
		for(SNP s : all_snps) {
			if(s.sameAs(snp))
				return true;
		}
		
		return false;
	}
	
	private int comparePositions(int nxt_pos, int prev_pos, List<SNP> snps) {
		
		for(int i = 0; i < snps.size(); i++) {
			SNP s = snps.get(i);
			if(s.getPosition() > prev_pos && s.getPosition() <= nxt_pos) {
				return s.getPosition();
			}
		}
		
		return nxt_pos;
	}
	
	/*
	 * For testing the output of the different WindowStats objects
	 * Not called at all hence the SuppressWarning
	 */
	@SuppressWarnings("unused")
	private String printLists(List<SNP> snps, List<Double> stats) {
		
		StringBuilder sb = new StringBuilder();
		
		sb.append("SNPS:\t" + snps.size() + "\n");
		for(int i = 0; i < snps.size(); i++) {
			sb.append(snps.get(i) + "\n");
		}
		
		sb.append("Scores:\t" + stats.size() + "\n");
		for(int i = 0; i < stats.size(); i++) {
			sb.append(stats.get(i) + "\n");
		}
		
		return sb.toString();
	}
	
	@Override
	public String toString() {
		
		StringBuilder sb = new StringBuilder();
		
		List<SNP> all_snps = getAllSNPs();
		
		for(int i = 0; i < all_snps.size(); i++) {
			
			SNP cur_snp = all_snps.get(i);

			sb.append(printSNP(cur_snp));
		}
		
		return sb.toString();
	}
	
	public String printSNP(SNP s) {
		
		StringBuilder sb = new StringBuilder();
		
		Double iHS_score = getScore(ihs_snps, ihs_stats, s);
		Double XPEHH_score = getScore(xpehh_snps, xpehh_stats, s);
		Double iHH_score = getScore(ihh_snps, ihh_stats, s);
		Double DDAF_score = getScore(ddaf_snps, ddaf_stats, s);
		Double DAF_score = getScore(daf_snps, daf_stats, s);
		Double Fst_score = getScore(fst_snps, fst_stats, s);
		//Double TAJD_score = getScore(tajd_snps, tajd_stats, s);
		//Double NEW_score = getScore(new_snps, new_stats, cur_snp);
		
		Double pop_score_std = win_scores_prod_std.get(s);
		Double mop_score_std = win_scores_mean_std.get(s);
		Double pop_score = win_scores_prod_unstd.get(s);
		Double mop_score = win_scores_mean_unstd.get(s);
		
		if(pop_score_std == null)
			pop_score_std = Double.NaN;
		if(mop_score_std == null)
			mop_score_std = Double.NaN;
		if(pop_score == null)
			pop_score = Double.NaN;
		if(mop_score == null)
			mop_score = Double.NaN;
		
		sb.append(s.getSnpID() + "\t");
		sb.append(s.getPosition() + "\t");
		sb.append(iHS_score + "\t");
		sb.append(XPEHH_score + "\t");
		sb.append(iHH_score + "\t");
		sb.append(DDAF_score + "\t");
		sb.append(DAF_score + "\t");
		sb.append(Fst_score + "\t");
		//sb.append(TAJD_score + "\t");
		//sb.append(NEW_score + "\t");
		
		sb.append(pop_score + "\t");
		sb.append(mop_score + "\t");
		sb.append(pop_score_std + "\t");
		sb.append(mop_score_std + "\n");
		
		return sb.toString();
	}

	@Override
	public int compareTo(WindowStats ws) {
		
		if(this.getStPos() < ws.getStPos())
			return -1;
		if(this.getStPos() > ws.getStPos())
			return 1;
		
		return 0;
	}
	
	
}
