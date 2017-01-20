package de.unijena.bioinf.em;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.special.Gamma;

public class MLBetaDistribution extends MLDistribution implements InterfaceMLDistribution {
	
	public MLBetaDistribution(boolean inverse, double shift) {
		super(inverse, shift);
		// TODO Auto-generated constructor stub
	}

	public String getDescription(){
		return "BetaDistribution"+(inverse?"Inverse":"");
	}

	public Map<String, Double> MStep(List<Double> expectation){
		Map<String, Double> parameters=new HashMap<String, Double>();

		double pi=EMUtils.getPi(expectation);
		parameters.put("pi", pi);

		double mu=EMUtils.getMu(sample, expectation);
		double sd=EMUtils.getSD(sample, expectation,mu);

		double alpha=Math.pow(mu,2)*(1-mu)/Math.pow(sd,2)-mu;
		double beta=Math.pow(1-mu,2)*mu/Math.pow(sd,2)-(1-mu);
		
		parameters.put("alpha", alpha);
		parameters.put("beta",beta);
		
		return parameters;
	}

	@Override
	public double getNumericalMean(Map<String, Double> parameters) throws Exception {
		return new BetaDistribution(parameters.get("alpha"), parameters.get("beta")).getNumericalMean();
	}
	
	@Override
	public double getDensityUnprocessed(double x, Map<String, Double> parameters) throws Exception {
		return new BetaDistribution(parameters.get("alpha"), parameters.get("beta")).density(x);
	}
	
	@Override
	public double getCumulativeProbabilityUnprocessed(double x1, Map<String, Double> parameters) throws Exception {
		return new BetaDistribution(parameters.get("alpha"), parameters.get("beta")).cumulativeProbability(x1);
	}
}
