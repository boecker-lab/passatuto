package de.unijena.bioinf.em;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.LogisticDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

public class MLLogisticDistribution extends MLDistribution implements InterfaceMLDistribution {

	public MLLogisticDistribution(boolean inverse, double shift) {
		super(inverse, shift);
		// TODO Auto-generated constructor stub
	}

	public String getDescription(){
		return "LogisticDistribution"+(inverse?"Inverse":"");
	}

	public Map<String, Double> MStep(List<Double> expectation){
		Map<String, Double> parameters=new HashMap<String, Double>();

		double pi=EMUtils.getPi(expectation);
		parameters.put("pi", pi);
		
		double mu=EMUtils.getMu(sample, expectation);
		parameters.put("mu",mu);

		double sd=EMUtils.getSD(sample, expectation, mu);
		double s=Math.pow(Math.pow(sd, 2)*3/Math.pow(Math.PI, 2),0.5);
		parameters.put("s", s);
		
		return parameters;
	}
	
	@Override
	public double getNumericalMean(Map<String, Double> parameters) throws Exception {
		return new LogisticDistribution(parameters.get("mu"), parameters.get("s")).getNumericalMean();
	}

	@Override
	public double getDensityUnprocessed(double x, Map<String, Double> parameters) {
		return new LogisticDistribution(parameters.get("mu"), parameters.get("s")).density(x);
	}
	
	@Override
	public double getCumulativeProbabilityUnprocessed(double x1, Map<String, Double> parameters) {
		return new LogisticDistribution(parameters.get("mu"), parameters.get("s")).cumulativeProbability(x1);
	}
}
