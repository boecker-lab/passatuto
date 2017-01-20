package de.unijena.bioinf.em;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

public class MLExponentialDistribution extends MLDistribution implements InterfaceMLDistribution {

	public MLExponentialDistribution(boolean inverse, double shift) {
		super(inverse, shift);
		// TODO Auto-generated constructor stub
	}

	public String getDescription(){
		return "ExponentialDistribution"+(inverse?"Inverse":"");
	}

	public Map<String, Double> MStep(List<Double> expectation){
		Map<String, Double> parameters=new HashMap<String, Double>();

		double pi=EMUtils.getPi(expectation);
		parameters.put("pi", pi);

		double mu=EMUtils.getMu(sample, expectation);
		parameters.put("mu", mu);
		
		return parameters;
	}

	@Override
	public double getNumericalMean(Map<String, Double> parameters) throws Exception {
		return new ExponentialDistribution(parameters.get("mu")).getNumericalMean();
	}
	
	@Override
	public double getDensityUnprocessed(double x, Map<String, Double> parameters) throws Exception {
		return new ExponentialDistribution(parameters.get("mu")).density(x);
	}
	
	@Override
	public double getCumulativeProbabilityUnprocessed(double x1, Map<String, Double> parameters) throws Exception {
		return new ExponentialDistribution(parameters.get("mu")).cumulativeProbability(x1);
	}
}
