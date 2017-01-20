package de.unijena.bioinf.em;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.NormalDistribution;

public class MLNormalDistribution extends MLDistribution implements InterfaceMLDistribution {

	public MLNormalDistribution(boolean inverse, double shift) {
		super(inverse, shift);
		// TODO Auto-generated constructor stub
	}

	public String getDescription(){
		return "NormalDistribution"+(inverse?"Inverse":"");
	}

	public Map<String, Double> MStep(List<Double> expectation){
		Map<String, Double> parameters=new HashMap<String, Double>();

		double pi=EMUtils.getPi(expectation);
		parameters.put("pi", pi);
		
		double mu=EMUtils.getMu(sample, expectation);
		parameters.put("mu",mu);

		double sd=EMUtils.getSD(sample, expectation, mu);
		parameters.put("sd", Math.max(1E-10,sd));
		
		return parameters;
	}
	
	@Override
	public double getNumericalMean(Map<String, Double> parameters) throws Exception {
		return new NormalDistribution(parameters.get("mu"), parameters.get("sd")).getNumericalMean();
	}

	@Override
	public double getDensityUnprocessed(double x, Map<String, Double> parameters) {
		return new NormalDistribution(parameters.get("mu"), parameters.get("sd")).density(x);
	}

	@Override
	public double getCumulativeProbabilityUnprocessed(double x1, Map<String, Double> parameters) {
		return new NormalDistribution(parameters.get("mu"), parameters.get("sd")).cumulativeProbability(x1);
	}
	
}
