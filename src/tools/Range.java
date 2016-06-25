package tools;

import java.io.Serializable;

/**
 * A wrapper for a range of positions in a chromsome, used
 * to track genetic map information which is often sparse.
 *
 */
@SuppressWarnings("serial")
public class Range implements Serializable {
	
	private int start;
	private int end;
	private int index;
	
	/**
	 * Creates range object
	 * 
	 * @param st		start of the range
	 * @param end		end of the range
	 * @param index		an index in the range TODO What is this for?
	 */
	public Range(int st, int end, int index) {
		this.start = st;
		this.end = end;
		this.index = index;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}
	
	public int getIndex() {
		return index;
	}
	
	/**
	 * Discovers if the given position lies within the range
	 * 
	 * @param pos	position to test
	 * @return		True if pos lies within range
	 */
	public boolean contains(int pos) {
		if (pos >= start && pos <= end) {
			return true;
		}
		
		return false;
	}
	
	/**
	 * Gets the length of the range (end - start)
	 * 
	 * @return	integer length of the range
	 */
	public int getPhysRange() {
		return end - start;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + end;
		result = prime * result + start;
		result = prime * result + index;
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
		Range other = (Range) obj;
		if (end != other.end) {
			return false;
		}
		if (start != other.start) {
			return false;
		}
		if (index != other.index) {
			return false;
		}
		return true;
	}

	@Override
	public String toString() {
		return "Range [st=" + start + ", end=" + end + ", index=" + index + "]";
	}
}
