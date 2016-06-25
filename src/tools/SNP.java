package tools;

import java.io.Serializable;

/**
 * Class to store SNP (Single Nucleotide Polymorphism) data
 *
 */
@SuppressWarnings("serial")
public class SNP implements Comparable<SNP>, Serializable {
	
	private int pos;
	
	private String snp_id;
	private String allele0;
	private String allele1;
	
	/**
	 * Simple constructor
	 */
	public SNP() {
		pos = -1;
		allele0 = "";
		allele1 = "";
		snp_id = "";
	}
	
	/**
	 * Constructor with fields for allele position and ID
	 * 
	 * @param pos		position of the SNP
	 * @param snp_id	ID of the SNP
	 */
	public SNP(int pos, String snp_id) {
		
		this.pos = pos;
		this.snp_id = snp_id;
		allele0 = "";
		allele1 = "";
	}
	
	/**
	 * Constructor with fields for allele position, reference allele, 
	 * alternate allele, and SNP ID
	 * 
	 * @param pos		allele position
	 * @param a0		reference allele
	 * @param a1		alternate allele
	 * @param snp_id	ID of the SNP
	 */
	public SNP(int pos, String a0, String a1, String snp_id) {

		this.pos = pos;
		this.allele0 = a0;
		this.allele1 = a1;
		this.snp_id = snp_id;
	}
	
	/**
	 * 
	 * @return	the reference allele
	 */
	public String getAllele0() {
		return allele0;
	}
	
	/**
	 * 
	 * @return the alternate allele
	 */
	public String getAllele1() {
		return allele1;
	}
	
	public int getPosition() {
		return pos;
	}
	
	public String getSnpID() {
		return snp_id;
	}
	
	/**
	 * Compares the given SNP to this one. Return true if
	 * given SNP is the same as this one in position, reference 
	 * allele, and alternate allele or if the positions match and 
	 * the reference alleles match the alternate allele (they've 
	 * been switched), but neither contains the flag "-ref".
	 * 
	 * @param s		SNP to compare
	 * @return		True if given SNP is the same as thing one in position, 
	 */	
	public boolean sameAs(SNP s) {
		
		//checks position and allele values
		if (s.getPosition() == pos
				&& s.getAllele0().equals(allele0)
				&& s.getAllele1().equals(allele1)) {
			return true;
		}
		//ensures the allele values aren't switched
		if (s.getPosition() == pos
				&& s.getAllele0().equals(allele1)
				&& s.getAllele1().equals(allele0)
				&& !getSnpID().contains("-ref")
				&& !s.getSnpID().contains("-ref")) {
			return true;
		}
		
		return false;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((allele0 == null) ? 0 : allele0.hashCode());
		result = prime * result + ((allele1 == null) ? 0 : allele1.hashCode());
		result = prime * result + pos;
		result = prime * result + ((snp_id == null) ? 0 : snp_id.hashCode());
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
		SNP other = (SNP) obj;
		if (allele0 == null) {
			if (other.allele0 != null) {
				return false;
			}
		} else if (!allele0.equals(other.allele0)) {
			return false;
		}
		if (allele1 == null) {
			if (other.allele1 != null) {
				return false;
			}
		} else if (!allele1.equals(other.allele1)) {
			return false;
		}
		if (pos != other.pos) {
			return false;
		}
		if (snp_id == null) {
			if (other.snp_id != null) {
				return false;
			}
		} else if (!snp_id.equals(other.snp_id)) {
			return false;
		}
		return true;
	}

	@Override
	public String toString() {
		return "SNP [pos=" + pos + ", snp_id=" + snp_id + ", a0=" + allele0 
				+ ", a1=" + allele1 + "]";
	}

	@Override
	public int compareTo(SNP s) {
		
		if (this.sameAs(s)) {
			return 0;
		}
		if (this.pos < s.getPosition()) {
			return -1;
		}
		if (this.pos > s.getPosition()) {
			return 1;
		}
		
		return this.snp_id.compareTo(s.getSnpID());
	}
}
