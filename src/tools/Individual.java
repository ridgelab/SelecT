package tools;

import java.io.Serializable;
import java.util.Arrays;

/**
 * Stores information for individuals from the input data.
 */
@SuppressWarnings("serial")
public class Individual implements Serializable {
	
	private final int ARRAY_SIZE = 256; // 32 bytes at a time
	
	private int id;
	private int strand1_index;
	private int strand2_index;
	private int chr;
	private boolean[] strand1;
	private boolean[] strand2;
	
	/**
	 * Simple constructor
	 */
	public Individual() {
		
		id = 0;
		chr = 0;
		strand1_index = 0;
		strand2_index = 0;
		strand1 = new boolean[ARRAY_SIZE];
		strand2 = new boolean[ARRAY_SIZE];
	}
	
	/**
	 * Constructor with arguments. Takes in the individual's 
	 * ID and the chromosome number of interest.
	 * 
	 * @param id		individual ID
	 * @param chr		chromsosome number
	 */
	public Individual(int id, int chr) {
		
		this.id = id;
		this.chr = chr;
		strand1_index = 0;
		strand2_index = 0;
		strand1 = new boolean[ARRAY_SIZE];
		strand2 = new boolean[ARRAY_SIZE];
	}
	//**************************************************************************//
	public boolean addAlleleToStrand1(String a) {
		
		boolean valid = true;
		if (!isValidAllele(a)) {
			valid = false;
		}
		
		char allele = a.charAt(0);
		
		if ((strand1_index + 1) == strand1.length) {
			strand1 = increaseSize(strand1);
		}
		
		if (!valid) {
			strand1[strand1_index] = false;
		}
		if (allele == '0') {
			strand1[strand1_index] = false;
		}
		else if (allele == '1') {
			strand1[strand1_index] = true;
		}
		
		strand1_index++;
		return valid;
	}
	
	public boolean addAlleleToStrand2(String a) {
		
		boolean valid = true;
		if (!isValidAllele(a)) {
			valid = false;
		}
		
		char allele = a.charAt(0);
		
		if((strand2_index + 1) == strand2.length) {
			strand2 = increaseSize(strand2);
		}
		
		if (!valid) {
			strand2[strand2_index] = false;
		}
		else if (allele == '0') {
			strand2[strand2_index] = false;
		}
		else if (allele == '1') {
			strand2[strand2_index] = true;
		}
		
		strand2_index++;
		return valid;
	}
	
	public void addAlleleToStrand1(boolean allele) {
		
		if ((strand1_index + 1) == strand1.length) {
			strand1 = increaseSize(strand1);
		}
		
		strand1[strand1_index] = allele;
		strand1_index++;
	}
	
	public void addAlleleToStrand2(boolean allele) {
		
		if ((strand2_index + 1) == strand2.length) {
			strand2 = increaseSize(strand2);
		}
		
		strand2[strand2_index] = allele;
		strand2_index++;
	}
	

	private boolean isValidAllele(String a) {
	    
        if (a.length() != 1) {
            return false;
        }
    
        char ch = a.charAt(0);
        if (!Character.isDigit(ch)) {
            return false;
        }
    
        if (ch != '1' && ch != '0') {
            return false;
        }
    
        return true;
    }   

	//********************************************************************************//
	/**
	 * Checks validity of adding an allele to the specified strand
	 * 
	 * @param allele			allele character of allele to add
	 * @param addToLeading		addToLeading flag specifying which strand receives the allele. True if leading strand.
	 * @return 					True if allele was valid 
	 */
	public boolean addAlleleToStrand(char allele, boolean addToLeading ) {
		boolean[] strand =  null;
		int strand_index = -1;
		
		if (addToLeading) {
			strand = this.strand1;
			strand_index = this.strand1_index;
		}
		else {
			strand = this.strand2;
			strand_index = this.strand2_index;
		}
			
		
		boolean valid = true;
		if (!isValidAllele(allele)) {
			valid = false;
		}
		
		if ((strand_index + 1) == strand.length) {
			strand = increaseSize(strand);
		}
		
		if (!valid) {
			strand[strand_index] = false;
		}
		if (allele == '0') {
			strand[strand_index] = false;
		}
		else if (allele == '1') {
			strand[strand_index] = true;
		}
		
		strand_index++;
		
		if (addToLeading) {
			this.strand1 = strand;
			this.strand1_index = strand_index;
		}
		else {
			this.strand2 = strand;
			this.strand2_index = strand_index;
		}
		return valid;
	}
	
	/**
	 * Adds an allele to the specified strand
	 * 
	 * @param allele			True for the paternal allele, false for maternal
	 * @param addToLeading		True to add the leading strand
	 */
	public void addAlleleToStrand(boolean allele, boolean addToLeading) {
		
		boolean[] strand =  null;
		int strand_index = -1;
		
		if (addToLeading) {
			strand = this.strand1;
			strand_index = this.strand1_index;
		}
		else {
			strand = this.strand2;
			strand_index = this.strand2_index;
		}
		
		if ((strand_index + 1) == strand.length) {
			strand = increaseSize(strand);
		}
		
		strand[strand_index] = allele;
		strand_index++;
		
		if (addToLeading) {
			this.strand1 = strand;
			this.strand1_index = strand_index;
		}
		else {
			this.strand2 = strand;
			this.strand2_index = strand_index;
		}
	}
	
	/**
	 * Gets the allele from the specified location in the desired strand
	 * 
	 * @param index				index to get allele from
	 * @param getFromLeading	True if allele from leading strand is desired
	 * @return					allele as int
	 */	
	public int getAlleleFromStrand(int index, boolean getFromLeading) {
		
		boolean allele = false;
		if (getFromLeading) {
			allele = strand1[index];
		}
		else {
			allele = strand2[index];
		}
		
		if (allele) {
			return 1;
		}
		else {
			return 0;
		}
	}
	
	private boolean[] increaseSize(boolean[] str) {
		boolean[] new_str = new boolean[str.length * 2];
		
		for (int i = 0; i < str.length; i++) {
			new_str[i] = str[i];
		}
		
		return new_str;
	}
	
	private boolean isValidAllele(char allele) {

		if (!Character.isDigit(allele)) {
			return false;
		}
		
		if (allele != '1' && allele != '0') {
			return false;
		}
		
		return true;
	}
	
	public int getStrandSize() {
		return strand1.length;
	}
	
	public int getID() {
		return id;
	}
	
	public void setID(int id) {
		this.id = id;
	}
	
	public int getChr() {
		return chr;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Individual [id=" + id + ", chr=" + chr + "]\n");
		
		sb.append("\tStrand1 [");
		for (int i = 0; i < strand1_index; i++) {
			sb.append(getAlleleFromStrand(i, true));
		}
		sb.append("] size=" + strand1_index + "\n");
		
		sb.append("\tStrand2 [");
		for (int i = 0; i < strand2_index; i++) {
			sb.append(getAlleleFromStrand(i, false));
		}
		sb.append("] size=" + strand2_index + "\n");
		
		return sb.toString();
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + chr;
		result = prime * result + id;
		result = prime * result + Arrays.hashCode(strand1);
		result = prime * result + Arrays.hashCode(strand2);
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
		Individual other = (Individual) obj;
		if (id != other.id) {
			return false;
		}
		if (chr != other.chr) {
			return false;
		}
		if (!Arrays.equals(strand1, other.strand1)) {
			return false;
		}
		if (!Arrays.equals(strand2, other.strand2)) {
			return false;
		}
		return true;
	}
}
