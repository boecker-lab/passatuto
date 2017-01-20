package de.unijena.bioinf.em;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public interface InterfaceMLDistribution {
	
	public Map<String, Double> MStep(List<Double> expectation);
	public List<Double> EStep(Map<String, Double> parameters) throws Exception;
	double getDensity(double x, Map<String, Double> parameters) throws Exception;
	public double getDensityUnprocessed(double x, Map<String, Double> parameters) throws Exception;
	public String getDescription();
	public void preprocessSample(List<Double> sample);	
	double getProbability(double x, Map<String, Double> parameters) throws Exception;
	double getCumulativeProbabilityUnprocessed(double x, Map<String, Double> parameters) throws Exception;
	double getNumericalMean(Map<String, Double> parameters) throws Exception;
	
}
