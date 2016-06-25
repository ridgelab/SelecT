package tools;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * Class to represent the windows that regions are broken 
 * up into for statistical calculations.
 * 
 */
@SuppressWarnings("serial")
public class Window implements Serializable {
	
	private int st_pos;
	private int st_index;
	private int end_pos;
	private int end_index;
	
	private List<SNP> all_snps;
	
	/**
	 * Simple constructor
	 */
	public Window() {
		st_pos = 0;
		end_pos = 0;
		all_snps = new ArrayList<SNP>();
	}
	
	/**
	 * Constructor with fields for start and end of window
	 * 
	 * @param st_pos	start position for window
	 * @param end_pos	end position for window
	 */
	public Window(int st_pos, int end_pos) {

		this.st_pos = st_pos;
		this.end_pos = end_pos;
	
		st_index = -1;
		end_index = -1;
		
		all_snps = new ArrayList<SNP>();
	}
	
	/**
	 * Constructor with fields for start of window, end of window, and 
	 * a start index to show where to begin analysis
	 * 
	 * @param st_pos	start position for window
	 * @param end_pos	end position for window
	 * @param st_index	start index to begin analysis
	 */
	public Window (int st_pos, int end_pos, int st_index) {
		
		this.st_pos = st_pos;
		this.end_pos = end_pos;
		this.st_index = st_index;
		this.all_snps = new ArrayList<SNP>();
	}
	
	/**
	 * Builds a SNP and adds it to the window
	 * 
	 * @param pos		position to add the SNP
	 * @param a0		reference allele
	 * @param a1		alternate allele
	 * @param snp_id	ID of SNP
	 */
	public void addSNP(int pos, String a0, String a1, String snp_id) {
		all_snps.add(new SNP(pos, a0, a1, snp_id));
	}
	
	/**
	 * Adds a SNP to the Window
	 * 
	 * @param s		SNP to add
	 */
	public void addSNP(SNP s) {
		all_snps.add(s);
	}
	
	/**
	 * Retrieves a SNP from the window
	 * 
	 * @param index		position to retrieve SNP from
	 * @return			returns null if this index isn't within the window's boundaries, the SNP otherwise
	 */
	public SNP getSNP(int index) {
		if (containsIndex(index)) {
			return all_snps.get(index - st_index);
		}
		else {
			return null;
		}
	}
	
	/**
	 * Retrieves a SNP from the window if the position, 
	 * reference allele, and alternate allele match
	 * 
	 * @param pos	position of desired SNP
	 * @param a0	reference allele
	 * @param a1	alternate allele
	 * @return		the desired SNP if found, null otherwise
	 */
	public SNP getSNP(int pos, String a0, String a1) {
		
		for (SNP s : all_snps) {
			if (s.getPosition() == pos) {
				if (s.getAllele0().equals(a0) && s.getAllele1().equals(a1)) {
					return s;
				}
				if (s.getAllele0().equals(a1) && s.getAllele1().equals(a0)) {
					return s;
				}
			}
		}
		return null;
	}
	
	/**
	 * Gets the index for a SNP from the window's individual array
	 * 
	 * @param snp	SNP whose index is is desired
	 * @return		the index if found, -1 otherwise
	 */
	public int getSnpIndex(SNP snp) {
		
		for (int i = 0; i < all_snps.size(); i++) {
			SNP s = all_snps.get(i);
			if (s.sameAs(snp)) {
				return st_index + i;	
			}
		}		
		return -1;		
	}
	
	/**
	 * Tests whether a SNP is represented or not
	 * 
	 * @param snp	SNP of interest
	 * @return		True if found
	 */
	public boolean containsSNP(SNP snp) {
		
		for (SNP s : all_snps) {
			if (s.sameAs(snp)) {
				return true;
			}
		}
		
		return false;
	}
	
	/**
	 * Tests whether an index is represented of not
	 * 
	 * @param index		index of interest
	 * @return			true if within the bounds of the window, false if not
	 */
	public boolean containsIndex(int index) {
		if (index >= st_index && index <= end_index) {
			return true;
		}
		return false;
	}
	
	public void setEndIndex(int index) {
		end_index = index;
	}
	
	public void setStIndex(int index) {
		st_index = index;
	}
	
	public int getStPos() {
		return st_pos;
	}

	public int getStIndex() {
		return st_index;
	}

	public int getEndPos() {
		return end_pos;
	}

	public int getEndIndex() {
		return end_index;
	}
	
	public List<SNP> getSNPs() {
		return all_snps;
	}
	
	public void setSNPs(List<SNP> new_snps) {
		all_snps = new_snps;
	}
	
	public int getSnpListSize() {
		return all_snps.size();
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result
				+ ((all_snps == null) ? 0 : all_snps.hashCode());
		result = prime * result + end_pos;
		result = prime * result + st_pos;
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		Window other = (Window) obj;
		if (all_snps == null) {
			if (other.all_snps != null) {
				return false;
			}
		} else if (!all_snps.equals(other.all_snps)) {
			return false;
		}
		if (end_pos != other.end_pos) {
			return false;
		}
		if (st_pos != other.st_pos) {
			return false;
		}
		return true;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Window [st_pos=" + st_pos + ", end_pos=" + end_pos 
				+ ", st_index=" + st_index + ", end_index=" + end_index + "] " 
				+ "size=" + all_snps.size() + "\n");
		
		for (int i = 0; i < all_snps.size(); i++) {
			sb.append("\t" + all_snps.get(i).toString() + "\n");
		}
		
		return sb.toString();
	}

	
}
