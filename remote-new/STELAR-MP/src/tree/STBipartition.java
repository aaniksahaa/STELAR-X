package tree;

import utils.BitSet;

public class STBipartition {
	public BitSet cluster1;
	public BitSet cluster2;
	private int _hash = 0;
	
	public STBipartition(BitSet c1, BitSet c2) {
		if (c1.nextSetBit(0) < c2.nextSetBit(0)) {
			cluster1 = (BitSet) c1.clone();
			cluster2 = (BitSet) c2.clone();
		} else {
			cluster1 = (BitSet) c2.clone();
			cluster2 = (BitSet) c1.clone();
		}
	}
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		if (!(obj instanceof STBipartition)) return false;
		STBipartition other = (STBipartition) obj;
		return cluster1.equals(other.cluster1) && cluster2.equals(other.cluster2);
	}
	
	@Override
	public int hashCode() {
		if (_hash == 0) {
			long[] words1 = cluster1.toLongArray();
			long[] words2 = cluster2.toLongArray();
			int result = 1;
			for (long word : words1) {
				result = 31 * result + (int)(word ^ (word >>> 32));
			}
			for (long word : words2) {
				result = 31 * result + (int)(word ^ (word >>> 32));
			}
			_hash = result;
		}
		return _hash;
	}
	
	@Override
	public String toString() {
		return cluster1.toString() + "|" + cluster2.toString();
	}
	
	public String print(String[] taxonIdToLabel) {
		StringBuilder sb1 = new StringBuilder();
		StringBuilder sb2 = new StringBuilder();
		
		sb1.append("(");
		boolean first1 = true;
		for (int i = cluster1.nextSetBit(0); i >= 0; i = cluster1.nextSetBit(i + 1)) {
			if (!first1) sb1.append(",");
			sb1.append(taxonIdToLabel[i]);
			first1 = false;
		}
		sb1.append(")");
		
		sb2.append("(");
		boolean first2 = true;
		for (int i = cluster2.nextSetBit(0); i >= 0; i = cluster2.nextSetBit(i + 1)) {
			if (!first2) sb2.append(",");
			sb2.append(taxonIdToLabel[i]);
			first2 = false;
		}
		sb2.append(")");
		
		return sb1.toString() + " | " + sb2.toString();
	}
	
	private boolean contains(BitSet container, BitSet contained) {
		BitSet temp = (BitSet) contained.clone();
		temp.andNot(container);
		return temp.isEmpty();
	}
	
	public boolean isDominatedBy(STBipartition dominant) {
		return (contains(dominant.cluster1, this.cluster1) && contains(dominant.cluster2, this.cluster2)) ||
				(contains(dominant.cluster2, this.cluster1) && contains(dominant.cluster1, this.cluster2));
	}
	
}
