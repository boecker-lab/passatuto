package de.unijena.bioinf.spectralcomparison;

import java.util.List;

import de.unijena.bioinf.decoy.model.MassBank;
import de.unijena.bioinf.statistics.Result;

public interface ComparisonMethod {

	public String getMethodNameLong();
	public void setSpectra(List<MassBank> s1, List<MassBank> s2);
//	public void setTrees(List<Tree> s1, List<Tree> s2);
	public List<MassBank>[] getSpectra();
	public Result getResult(List<double[]> peaksLeft, List<double[]> peaksRight);
	public AlignmentMatrix getAlignmentMatrix();
}
