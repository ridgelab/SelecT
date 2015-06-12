package tools;

import java.io.Serializable;

@SuppressWarnings("serial")
public class SNP implements Comparable<SNP>, Serializable {
	
	private int pos;
	
	private String snp_id;
	private String a0;
	private String a1;
	
	
	public SNP() {
		pos = -1;
		a0 = "";
		a1 = "";
		snp_id = "";
	}
	
	public SNP(int pos, String snp_id) {
		
		this.pos = pos;
		this.snp_id = snp_id;
		a0 = "";
		a1 = "";
	}
	
	public SNP(int pos, String a0, String a1, String snp_id) {

		this.pos = pos;
		this.a0 = a0;
		this.a1 = a1;
		this.snp_id = snp_id;
	}
	
	public String getAllele0() {
		return a0;
	}
	
	public String getAllele1() {
		return a1;
	}
	
	public int getPosition() {
		return pos;
	}
	
	public String getSnpID() {
		return snp_id;
	}
	
	public boolean sameAs(SNP s) {
		
		//checks position and allele values
		if(s.getPosition() == pos
				&& s.getAllele0().equals(a0)
				&& s.getAllele1().equals(a1)) {
			return true;
		}
		//ensures the allele values aren't switched
		if(s.getPosition() == pos
				&& s.getAllele0().equals(a1)
				&& s.getAllele1().equals(a0)
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
		result = prime * result + ((a0 == null) ? 0 : a0.hashCode());
		result = prime * result + ((a1 == null) ? 0 : a1.hashCode());
		result = prime * result + pos;
		result = prime * result + ((snp_id == null) ? 0 : snp_id.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		SNP other = (SNP) obj;
		if (a0 == null) {
			if (other.a0 != null)
				return false;
		} else if (!a0.equals(other.a0))
			return false;
		if (a1 == null) {
			if (other.a1 != null)
				return false;
		} else if (!a1.equals(other.a1))
			return false;
		if (pos != other.pos)
			return false;
		if (snp_id == null) {
			if (other.snp_id != null)
				return false;
		} else if (!snp_id.equals(other.snp_id))
			return false;
		return true;
	}

	@Override
	public String toString() {
		return "SNP [pos=" + pos + ", snp_id=" + snp_id + ", a0=" + a0 
				+ ", a1=" + a1 + "]";
	}

	@Override
	public int compareTo(SNP s) {
		
		if(this.equals(s))
			return 0;
		if(this.pos < s.getPosition())
			return -1;
		if(this.pos > s.getPosition())
			return 1;
		
		return this.snp_id.compareTo(s.getSnpID());
	}
}
