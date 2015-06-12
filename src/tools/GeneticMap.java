package tools;

import java.io.Serializable;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

@SuppressWarnings("serial")
public class GeneticMap implements Serializable {
	
	int rate_sum;
	
	private List<Range> st_pos;
	private Map<Range, Double> gen_map;
	
	public GeneticMap() {
		
		rate_sum = 0;
		
		gen_map = new HashMap<Range, Double>();
		st_pos = new LinkedList<Range>();
		
	}
	
	/**
	 * For use in adding to the wrapped Map called GeneticMap. This wrapped Map
	 * is necessary because the key must be a range of positions (start and end)
	 * with the key being the recombination rate between those two positions. 
	 * 
	 * It is important to note that the Genetic Map values MUST be inputed in
	 * an order fashion, meaning, you can not randomly insert recombination 
	 * rates but have to start from the smallest physical point to the largest 
	 * physical point.
	 * 
	 * @param st		Start position in base pair physical position (hg19)
	 * @param end		End position in base pair physical position (hg19)
	 * @param rate		Recombination rate in cM/Mb between the two positions
	 */
	public void put(int st, int end, double rate) {
		
		Range rng = new Range(st, end, st_pos.size());
		st_pos.add(new Range(st, end, st_pos.size()));
		gen_map.put(rng, new Double(rate));
	}
	
	/**
	 * This finds the genetic map recombination rate between to positions by 
	 * estimating a high-resolution genetic map
	 * 
	 * Assumes that recombination rate changes uniformly between data points
	 * in order to estimate higher resolution
	 * 
	 * @param Physical position of the start site
	 * @param Physical position of the target site
	 * @return The positive recombination rate between the two positions
	 */
public double getRecombRate(int up_pos, int dwn_pos) {
		
		//Case 1: positions are equal and by definition no recombination
		if(up_pos == dwn_pos)
			return 0.0;
		
		//Case 2: positions cannot be found in map (invalid or beyond range)
		Range dwn_rng = getRange(dwn_pos);
		Range up_rng = null;
		Range last_rng = st_pos.get(st_pos.size() - 1);
		
		if(dwn_rng != null)
			up_rng = getRangeFromDwnRng(up_pos, dwn_rng);
		else if(dwn_rng == null && dwn_pos > 0 && up_pos > 0) {
			
			double avg_rate = (double) rate_sum / last_rng.getEnd();
			return avg_rate * (double) (up_pos - dwn_pos);
		}
		else {
			System.out.println("Fatal Error: down position " + dwn_pos + " is unacceptable. "
					+ "Check genetic map for irregularities and api for more info");
			System.exit(0);
		}
		
		if(up_rng == null && dwn_pos > 0 && up_pos > 0) {
			
			double avg_rate = (double) rate_sum / last_rng.getEnd();
			return avg_rate * (double) (up_pos - dwn_pos);
		}
		else if(up_rng == null) {
			System.out.println("Fatal Error: up position " + up_pos + " is unacceptable. "
					+ "Check genetic map for irregularities and api for more info");
			System.exit(0);
		}
		
		//Case 3: recombination rate can be estimated within single GenMap entry
		if(up_rng.equals(dwn_rng))
			return getLocalizedRecombRate(up_rng, up_pos, dwn_pos);
		
		double dwn_recomb_rate = getLocalizedRecombRate(dwn_rng, dwn_rng.getEnd(), dwn_pos); 
		double up_recomb_rate = getLocalizedRecombRate(up_rng, up_pos, up_rng.getSt()); 
		
		//Case 4: recombination rate can be estimated with exactly 2 GenMap entries
		if(up_rng.getSt() == (dwn_rng.getEnd() + 1))
			return dwn_recomb_rate + up_recomb_rate;
		
		//Case 5: recombination rate is estimated with 2+ GenMap entries
		double tot_recomb_rate = up_recomb_rate + dwn_recomb_rate;
		int up_index = up_rng.getIndex();
		int dwn_index = dwn_rng.getIndex();
		
		for(int i = dwn_index + 1; i < up_index; i++) { 
			
			Range r = st_pos.get(i);
			tot_recomb_rate += gen_map.get(r);
		}
		
		return tot_recomb_rate;
	}
	
	
	public double get(int pos) {
		
		Range key = getRange(pos);
		
		if(key == null)
			return 0.0;
		
		return gen_map.get(key);
	}
	
	public double get(Range r) {
		return gen_map.get(r);
	}
	
	public Set<Range> getRangeSet() {
		return gen_map.keySet();
	}
	
	public Range getRange(int pos) {
		
		for(Range r : gen_map.keySet()) {
			if(r.getSt() <= pos && r.getEnd() >= pos)
				return r;
		}
		
		return null;
	}
	
	private Range getRangeFromDwnRng(int pos, Range rng) {
		
		for(int i = rng.getIndex(); i < st_pos.size(); i++) {
			
			Range r = st_pos.get(i);
			if(r.getSt() <= pos && r.getEnd() >= pos)
				return r;
		}
		
		return null;
	}
	
	private double getLocalizedRecombRate(Range rng, int up_pos, int dwn_pos) {
		
		int dist = up_pos - dwn_pos;
		double dist_ratio = (double) dist / (double) rng.getPhysRange();
		double local_recomb_rate = gen_map.get(rng) * dist_ratio;
		
		return local_recomb_rate;
	}
}

