package tools;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import calc.HaplotypeTests;

public class WindowStats implements Comparable<WindowStats>{
	
	private int st_pos;
	private int end_pos;
	
	private TreeMap<SNP, Double> ihs;
	private TreeMap<SNP, Double> xpehh;
	private TreeMap<SNP, Double> ihh;
	private TreeMap<SNP, Double> ddaf;
	private TreeMap<SNP, Double> daf;
	private TreeMap<SNP, Double> fst;
	//private TreeMap<SNP, Double> tajd;
	//private TreeMap<SNP, Double> new;
	
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
		
		ihs = new TreeMap<SNP, Double>();
		xpehh = new TreeMap<SNP, Double>();
		ihh = new TreeMap<SNP, Double>();
		ddaf = new TreeMap<SNP, Double>();
		daf = new TreeMap<SNP, Double>();
		fst = new TreeMap<SNP, Double>();
		//tajd = new TreeMap<SNP, Double>();
		//new = new TreeMap<SNP, Double>();
		
		win_scores_prod_unstd = new TreeMap<SNP, Double>();
		win_scores_mean_unstd = new TreeMap<SNP, Double>();
		win_scores_prod_std = new TreeMap<SNP, Double>();
		win_scores_mean_std = new TreeMap<SNP, Double>();
	}
	
	public List<SNP> getAllSNPs() {
		
		TreeSet<SNP> snps_set = new TreeSet<SNP>();
		
		snps_set.addAll(ihs.keySet());
		snps_set.addAll(xpehh.keySet());
		snps_set.addAll(ihh.keySet());
		snps_set.addAll(ddaf.keySet());
		snps_set.addAll(daf.keySet());
		snps_set.addAll(fst.keySet());
		//snps_set.addAll(tajd.keySet());
		//snps_set.addAll(new.keySet());
		
		List<SNP> all_snps = new LinkedList<SNP>();
		
		for(SNP s : snps_set)
			all_snps.add(s);
		
		return all_snps;
	}
	
	public int getTotNumSNPs() {
		
		List<SNP> all_snps = getAllSNPs();
		
		return all_snps.size();
	}
	
	public int getNextPosition(int prev_pos) {
		
		int nxt_pos = end_pos;
		
		nxt_pos = comparePositions(nxt_pos, prev_pos, ihs.keySet());
		nxt_pos = comparePositions(nxt_pos, prev_pos, xpehh.keySet());
		nxt_pos = comparePositions(nxt_pos, prev_pos, ihh.keySet());
		nxt_pos = comparePositions(nxt_pos, prev_pos, ddaf.keySet());
		nxt_pos = comparePositions(nxt_pos, prev_pos, daf.keySet());
		nxt_pos = comparePositions(nxt_pos, prev_pos, fst.keySet());
		//nxt_pos = comparePositions(nxt_pos, prev_pos, tajd.keySet());
		//nxt_pos = comparePositions(nxt_pos, prev_pos, new.keySet());
		
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
		
		List<SNP> ihs_snps = new ArrayList<SNP>();
		
		for(SNP s : ihs.keySet())
			ihs_snps.add(s);
		
		return ihs_snps;
	}

	public List<Double> getIHSstats() {
		
		List<Double> ihs_stats = new ArrayList<Double>();
		
		for(SNP s : ihs.keySet())
			ihs_stats.add(ihs.get(s));
		
		return ihs_stats;
	}

	/**
	 * Sets the WindowStats ihs tree to the Lists of stats and snps. Assumes List indexes indicate
	 * corresponding values of SNP and ihs score. If lists are of unequal length nothing happens.
	 * 
	 * @param ihs_stats
	 * @param ihs_snps
	 */
	public void setIHS(List<Double> ihs_stats, List<SNP> ihs_snps) {
		
		if(ihs_stats.size() == ihs_snps.size()) {
			
			ihs.clear();
			ihs = new TreeMap<SNP, Double>();
			
			for(int i = 0; i < ihs_snps.size(); i++) 
				ihs.put(ihs_snps.get(i), ihs_stats.get(i));
		}	
	}
	
	/**
	 * Adds to the WindowStats ihs tree to the Lists of stats and snps. Assumes List indexes indicate
	 * corresponding values of SNP and ihs score. If lists are of unequal length nothing happens.
	 * 
	 * @param ihs_stats
	 * @param ihs_snps
	 */
	public void addIHS(List<Double> ihs_stats, List<SNP> ihs_snps) {
		
		if(ihs_stats.size() == ihs_snps.size()) {
			
			for(int i = 0; i < ihs_snps.size(); i++) 
				ihs.put(ihs_snps.get(i), ihs_stats.get(i));
		}
	}

	public List<SNP> getXPEHHsnps() {
		
		List<SNP> xpehh_snps = new ArrayList<SNP>();
		
		for(SNP s : xpehh.keySet())
			xpehh_snps.add(s);
		
		return xpehh_snps;
	}

	public List<Double> getXPEHHstats() {
		
		List<Double> xpehh_stats = new ArrayList<Double>();
		
		for(SNP s : xpehh.keySet())
			xpehh_stats.add(xpehh.get(s));
		
		return xpehh_stats;
	}

	public void setXPEHH(List<Double> xpehh_stats, List<SNP> xpehh_snps) {
		
		if(xpehh_stats.size() == xpehh_snps.size()) {
			
			xpehh.clear();
			xpehh = new TreeMap<SNP, Double>();
			
			for(int i = 0; i < xpehh_snps.size(); i++) 
				xpehh.put(xpehh_snps.get(i), xpehh_stats.get(i));
		}	
	}
	
	public void addXPEHH(List<Double> xpehh_stats, List<SNP> xpehh_snps) {
		
		if(xpehh_stats.size() == xpehh_snps.size()) {
			
			for(int i = 0; i < xpehh_snps.size(); i++) 
				xpehh.put(xpehh_snps.get(i), xpehh_stats.get(i));
		}
	}

	public List<SNP> getIHHsnps() {
		
		List<SNP> ihh_snps = new ArrayList<SNP>();
		
		for(SNP s : ihh.keySet())
			ihh_snps.add(s);
		
		return ihh_snps;
	}
	
	public List<Double> getIHHstats() {
		
		List<Double> ihh_stats = new ArrayList<Double>();
		
		for(SNP s : ihh.keySet())
			ihh_stats.add(ihh.get(s));
		
		return ihh_stats;
	}

	public void setIHH(List<Double> ihh_stats, List<SNP> ihh_snps) {
		
		if(ihh_stats.size() == ihh_snps.size()) {
			
			ihh.clear();
			ihh = new TreeMap<SNP, Double>();
			
			for(int i = 0; i < ihh_snps.size(); i++)
				ihh.put(ihh_snps.get(i), ihh_stats.get(i));
		}
	}
	
	public void addIHH(List<Double> ihh_stats, List<SNP> ihh_snps) {
		
		if(ihh_stats.size() == ihh_snps.size()) {
			
			for(int i = 0; i < ihh_snps.size(); i++)
				ihh.put(ihh_snps.get(i), ihh_stats.get(i));
		}
	}

	public List<SNP> getDDAFsnps() {
		
		List<SNP> ddaf_snps = new ArrayList<SNP>();
		
		for(SNP s : ddaf.keySet())
			ddaf_snps.add(s);
		
		return ddaf_snps;
	}

	public List<Double> getDDAFstats() {
		
		List<Double> ddaf_stats = new ArrayList<Double>();
		
		for(SNP s : ddaf.keySet())
			ddaf_stats.add(ddaf.get(s));
		
		return ddaf_stats;
	}

	public void setDDAF(List<Double> ddaf_stats, List<SNP> ddaf_snps) {
		
		if(ddaf_stats.size() == ddaf_snps.size()) {
			
			ddaf.clear();
			ddaf = new TreeMap<SNP, Double>();
			
			for(int i = 0; i < ddaf_snps.size(); i++)
				ddaf.put(ddaf_snps.get(i), ddaf_stats.get(i));
		}
	}
	
	public void addDDAF(List<Double> ddaf_stats, List<SNP> ddaf_snps) {
		
		if(ddaf_stats.size() == ddaf_snps.size()) {
			
			for(int i = 0; i < ddaf_snps.size(); i++)
				ddaf.put(ddaf_snps.get(i), ddaf_stats.get(i));
		}
	}
	
	public List<SNP> getDAFsnps() {
		
		List<SNP> daf_snps = new ArrayList<SNP>();
		
		for(SNP s : daf.keySet())
			daf_snps.add(s);
		
		return daf_snps;
	}

	public List<Double> getDAFstats() {
		
		List<Double> daf_stats = new ArrayList<Double>();
		
		for(SNP s : daf.keySet())
			daf_stats.add(daf.get(s));
		
		return daf_stats;
	}

	public void setDAF(List<Double> daf_stats, List<SNP> daf_snps) {
		
		if(daf_stats.size() == daf_snps.size()) {
			
			daf.clear();
			daf = new TreeMap<SNP, Double>();
			
			for(int i = 0; i < daf_snps.size(); i++)
				daf.put(daf_snps.get(i), daf_stats.get(i));
		}
	}
	
	public void addDAF(List<Double> daf_stats, List<SNP> daf_snps) {
		
		if(daf_stats.size() == daf_snps.size()) {
			
			for(int i = 0; i < daf_snps.size(); i++)
				daf.put(daf_snps.get(i), daf_stats.get(i));
		}
	}
	
	public List<SNP> getFSTsnps() {
		
		List<SNP> fst_snps = new ArrayList<SNP>();
		
		for(SNP s : fst.keySet())
			fst_snps.add(s);
		
		return fst_snps;
	}
	
	public List<Double> getFSTstats() {
		
		List<Double> fst_stats = new ArrayList<Double>();
		
		for(SNP s : fst.keySet())
			fst_stats.add(fst.get(s));
		
		return fst_stats;
	}

	public void setFST(List<Double> fst_stats, List<SNP> fst_snps) {
		
		if(fst_stats.size() == fst_snps.size()) {
			
			fst.clear();
			fst = new TreeMap<SNP, Double>();
			
			for(int i = 0; i < fst_snps.size(); i++)
				fst.put(fst_snps.get(i), fst_stats.get(i));
		}
	}
	
	public void addFST(List<Double> fst_stats, List<SNP> fst_snps) {
		
		if(fst_stats.size() == fst_snps.size()) {
			
			for(int i = 0; i < fst_snps.size(); i++)
				fst.put(fst_snps.get(i), fst_stats.get(i));
		}
	}
	
//	public List<SNP> getTAJDsnps() {
//		
//		List<SNP> tajd_snps = new ArrayList<SNP>();
//		
//		for(SNP s : tajd.keySet())
//			tajd_snps.add(s);
//		
//		return tajd_snps;
//	}
//
//	public List<Double> getTAJDstats() {
//		
//		List<Double> tajd_stats = new ArrayList<Double>();
//		
//		for(SNP s : tajd.keySet())
//			tajd_stats.add(tajd.get(s));
//		
//		return tajd_stats;
//	}
//
//	public void setTAJD(List<Double> tajd_stats, List<SNP> tajd_snps) {
//		
//		if(tajd_stats.size() == tajd_snps.size()) {
//			
//			tajd.clear();
//			tajd = new TreeMap<SNP, Double>();
//			
//			for(int i = 0; i < tajd_snps.size(); i++)
//				tajd.put(tajd_snps.get(i), tajd_stats.get(i));
//		}
//	}
//	
//	public void addTAFD(List<Double> tajd_stats, List<SNP> tajd_snps) {
//		
//		if(tajd_stats.size() == tajd_snps.size()) {
//			
//			for(int i = 0; i < tajd_snps.size(); i++)
//				tajd.put(tajd_snps.get(i), tajd_stats.get(i));
//		}
//	}
	
//	public List<SNP> getNEWsnps() {
//		
//		List<SNP> new_snps = new ArrayList<SNP>();
//		
//		for(SNP s : new.keySet())
//			new_snps.add(s);
//		
//		return new_snps;
//	}
//
//	public List<Double> getNEWstats() {
//		List<Double> new_stats = new ArrayList<Double>();
//		
//		for(SNP s : new.keySet())
//			new_stats.add(new.get(s));
//		
//		return new_stats;
//	}
//
//	public void setNEW(List<Double> new_stats, List<SNP> new_snps) {
//		
//		if(new_stats.size() == new_snps.size()) {
//			
//			new.clear();
//			new = new TreeMap<SNP, Double>();
//			
//			for(int i = 0; i < new_snps.size(); i++)
//				new.put(new_snps.get(i), new_stats.get(i));
//		}
//	}
//	
//	public void addNEW(List<Double> new_stats, List<SNP> new_snps) {
//		
//		if(new_stats.size() == new_snps.size()) {
//			
//			for(int i = 0; i < new_snps.size(); i++)
//				new.put(new_snps.get(i), new_stats.get(i));
//		}
//	}
	
	public Double getIhsScore(SNP snp) {
		
		Double score = ihs.get(snp);
		
		if(score == null)
			return Double.NaN;
		else
			return score;
	}
	
	public Double getIhhScore(SNP snp) {
		
		Double score = ihh.get(snp);
		
		if(score == null)
			return Double.NaN;
		else
			return score;
	}
	
	public Double getXpehhScore(SNP snp) {
		
		Double score = xpehh.get(snp);
		
		if(score == null)
			return Double.NaN;
		else
			return score;
	}
	
	public Double getDDafScore(SNP snp) {
		
		Double score = ddaf.get(snp);
		
		if(score == null)
			return Double.NaN;
		else
			return score;
	}
	
	public Double getDafScore(SNP snp) {
		
		Double score = daf.get(snp);
		
		if(score == null)
			return Double.NaN;
		else
			return score;
	}
	
	public Double getFstScore(SNP snp) {
		
		Double score = fst.get(snp);
		
		if(score == null)
			return Double.NaN;
		else
			return score;
	}

//	public Double getTajdScore(SNP snp) {
//		
//		Double score = tajd.get(snp);
//		
//		if(score == null)
//			return Double.NaN;
//		else
//			return score;
//	}
	
//	public Double getNewScore(SNP snp) {
//		
//		Double score = new.get(snp);
//		
//		if(score == null)
//			return Double.NaN;
//		else
//			return score;
//	}
	
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
	
	private int comparePositions(int nxt_pos, int prev_pos, Set<SNP> snps) {
		
		for(SNP s : snps) {
			if(s.getPosition() > prev_pos && s.getPosition() <= nxt_pos)
				nxt_pos = s.getPosition();
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
		
		Double iHS_score = getIhsScore(s);
		Double XPEHH_score = getXpehhScore(s);
		Double iHH_score = getIhhScore(s);
		Double DDAF_score = getDDafScore(s);
		Double DAF_score = getDafScore(s);
		Double Fst_score = getFstScore(s);
		//Double TAJD_score = getTajdScore(s);
		//Double NEW_score = getNewScore(s);
		
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
