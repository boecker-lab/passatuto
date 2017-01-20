package de.unijena.bioinf.em;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.GumbelDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

public class MLGumbelDistribution extends MLDistribution implements InterfaceMLDistribution {

	public MLGumbelDistribution(boolean inverse, double shift) {
		super(inverse, shift);
		// TODO Auto-generated constructor stub
	}

	public String getDescription(){
		return "GumbelDistribution"+(inverse?"Inverse":"");
	}

	public Map<String, Double> MStep(List<Double> expectation){
		Map<String, Double> parameters=new HashMap<String, Double>();

		double pi=EMUtils.getPi(expectation);
		parameters.put("pi", pi);

		double mu=EMUtils.getMu(sample, expectation);
//		double sd=EMUtils.getSD(sample, expectation, mu);
		double musqr=getMuSqr(sample, expectation);
		
//		double beta=Math.pow(6*Math.pow(sd, 2)/Math.pow(Math.PI,2),0.5);
		double beta=Math.pow((musqr-Math.pow(mu,2))/(1.978-Math.pow(0.5772156649,2)), 0.5); //Revista de Matematica: Teoria y Aplicaciones 2005 12(1 & 2) : 151–156
		parameters.put("beta", beta);
				
		double loc=mu-beta*0.5772156649;
		parameters.put("loc",loc);
		
		return parameters;
	}

	@Override
	public double getNumericalMean(Map<String, Double> parameters) throws Exception {
		return new GumbelDistribution(parameters.get("loc"), parameters.get("beta")).getNumericalMean();
	}

	@Override
	public double getDensityUnprocessed(double x, Map<String, Double> parameters) {
		return new GumbelDistribution(parameters.get("loc"), parameters.get("beta")).density(x);
	}
	
	@Override
	public double getCumulativeProbabilityUnprocessed(double x1, Map<String, Double> parameters) {
		return new GumbelDistribution(parameters.get("loc"), parameters.get("beta")).cumulativeProbability(x1);
	}
	
	private static double getMuSqr(List<Double> sample, List<Double> expectation){
		double muDividend=0;
		double muDivisor=0;
		for(int j=0;j<expectation.size();j++){
			double e=expectation.get(j);
			double x=sample.get(j);
			muDividend+=e*Math.pow(x,2);
			muDivisor+=e;
		}
		double mu=muDividend/muDivisor;
		return mu;
	}	
}
