package de.unijena.bioinf.decoy.model;

import java.util.Comparator;

public class VertexComparator implements Comparator<DecoyTreeVertex>{

	@Override
	public int compare(DecoyTreeVertex o1, DecoyTreeVertex o2) {
		return Double.compare(o1.mf.getMass(),o2.mf.getMass());
	}
}
