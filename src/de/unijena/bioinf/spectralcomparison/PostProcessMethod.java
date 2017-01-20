package de.unijena.bioinf.spectralcomparison;

import de.unijena.bioinf.statistics.Result;

public interface PostProcessMethod {
	public Result[][] getResultList(ComparisonMethod cm);
	public String getMethodNameLong();
	
}
