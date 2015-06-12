package tools;

import java.io.Serializable;

@SuppressWarnings("serial")
public class Range implements Serializable {
	
	private int st;
	private int end;
	private int index;
	
	public Range(int st, int end, int index) {
		this.st = st;
		this.end = end;
		this.index = index;
	}

	public int getSt() {
		return st;
	}

	public int getEnd() {
		return end;
	}
	
	public int getIndex() {
		return index;
	}
	
	public boolean contains(int pos) {
		if(pos >= st && pos <= end) 
			return true;
		
		return false;
	}
	
	public int getPhysRange() {
		return end - st;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + end;
		result = prime * result + st;
		result = prime * result + index;
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
		Range other = (Range) obj;
		if (end != other.end)
			return false;
		if (st != other.st)
			return false;
		if (index != other.index)
			return false;
		return true;
	}

	@Override
	public String toString() {
		return "Range [st=" + st + ", end=" + end + ", index=" + index + "]";
	}
}
