package tools;

import java.util.LinkedList;

/**
 * ExtendedHaplotype stores a haplotype as a list of IDStrands.
 * The IDStrand is a wrapper containing an individual ID and an 
 * identifier showing the strand to be 1 or 2. 
 */
public class ExtendedHaplotype {

	/*
	 * use a FIFO method for organizing information
	 * 		use LinkedList .removeFirst()
	 */
	private LinkedList<IDStrand> haplo;
	
	/**
	 * Simple constructor
	 */
	public ExtendedHaplotype() {
		haplo = new LinkedList<IDStrand>();
	}	
	
	/**
	 * Adds a new individual and strand to the extended haplotype.
	 * 
	 * @param id		individual ID
	 * @param strand	strand identifier (1 or 2)
	 */
	public void add(int id, int strand) {
		haplo.add(new IDStrand(id, strand));
	}
	
	public int getNextID() {
		
		IDStrand top = haplo.getFirst();
		return top.getID();
	}
	
	public int getNextStrand() {
		
		IDStrand top = haplo.getFirst();
		return top.getStrand();
	}
	
	public void removeFirst() {
		haplo.removeFirst();
	}
	
	public boolean isEmpty() {
		return haplo.isEmpty();
	}
	
	public int size() {
		return haplo.size();
	}
	
	@Override
	public String toString() {
		
		StringBuilder sb = new StringBuilder();
		
		sb.append("ExtendedHaplotype [haplo.size = " + haplo.size() + "]\n");
		for(IDStrand ids : haplo) {
			sb.append("\t" + ids.toString() + "\n");
		}
		
		return sb.toString();
	}
	
	
	
	
//------------------------------------------------------------------------------
//---------------------------Package-only class---------------------------------	
//------------------------------------------------------------------------------
	
	//to access this data maybe I should just do two "get" commands...
	private class IDStrand {
		
		int id;
		int strand;
		
		public IDStrand(int id, int strand) {
			this.id = id;
			
			if (strand == 1 || strand == 2) {
				this.strand = strand;
			}
			else {
				strand = 0; //error
			}
		}
		
		public int getID() {
			return id;
		}
		
		public int getStrand() {
			return strand;
		}

		@Override
		public String toString() {
			return "IDStrand [id=" + id + ", strand=" + strand + "]";
		}
	}
//------------------------------------------------------------------------------	
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------	
}
