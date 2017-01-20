package de.unijena.bioinf.em;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.special.Gamma;

public class MLGammaDistribution extends MLDistribution implements InterfaceMLDistribution {
	
	public MLGammaDistribution(boolean inverse, double shift ) {
		super(inverse, shift);
		// TODO Auto-generated constructor stub
	}

	public String getDescription(){
		return "GammaDistribution"+(inverse?"Inverse":"");
	}

	public Map<String, Double> MStep(List<Double> expectation){
		Map<String, Double> parameters=new HashMap<String, Double>();

		double pi=EMUtils.getPi(expectation);
		parameters.put("pi", pi);

		double mu=EMUtils.getMu(sample, expectation);		
		double logmu=EMUtils.getLogMu(sample, expectation);
//		double sd=EMUtils.getSD(sample, expectation,mu);
//		double n=EMUtils.getN(expectation);
		double s=Math.log(mu)-logmu;

//		double scale=Math.pow(sd,2)/mu;
//		double shape=mu/scale;
		
//		double shape=1/Math.pow(sd/mu,2)-1d/n;
//		double scale=mu/shape;
		
		double shape=Double.POSITIVE_INFINITY;
		double shape2=(3-s+Math.pow(Math.pow(s-3,2)+24*s,0.5))/(12*s);
		int it=0;
		while(it<500&&Math.abs(shape-shape2)>1E-13){
			shape=shape2;
			shape2=shape-(Math.log(shape)-Gamma.digamma(shape)-s)/(1/shape-Gamma.trigamma(shape));
			it++;
		}
		shape=shape2;
		double scale=mu/shape;//wikipedia
		if(Double.isInfinite(shape)){
			System.out.println();
		}
		
		parameters.put("scale", scale);
		parameters.put("shape",shape);
		
		return parameters;
	}

	@Override
	public double getNumericalMean(Map<String, Double> parameters) throws Exception {
		return new GammaDistribution(parameters.get("shape"), parameters.get("scale")).getNumericalMean();
	}
	
	@Override
	public double getDensityUnprocessed(double x, Map<String, Double> parameters) throws Exception {
		return new GammaDistribution(parameters.get("shape"), parameters.get("scale")).density(x);
	}
	
	@Override
	public double getCumulativeProbabilityUnprocessed(double x1, Map<String, Double> parameters) throws Exception {
		return new GammaDistribution(parameters.get("shape"), parameters.get("scale")).cumulativeProbability(x1);
	}
}
