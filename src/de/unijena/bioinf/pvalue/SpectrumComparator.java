package de.unijena.bioinf.pvalue;

import java.util.Comparator;

class SpectrumComparator implements Comparator<double[]>{

	@Override
	public int compare(double[] o1, double[] o2) {
		return Double.compare(o1[0], o2[0]);
	}
	
}