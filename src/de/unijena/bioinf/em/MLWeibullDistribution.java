package de.unijena.bioinf.em;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.WeibullDistribution;
import org.apache.commons.math3.special.Gamma;

import umontreal.iro.lecuyer.probdist.WeibullDist;

public class MLWeibullDistribution extends MLDistribution implements InterfaceMLDistribution {

	public MLWeibullDistribution(boolean inverse, double shift) {
		super(inverse, shift);
		// TODO Auto-generated constructor stub
	}

	public String getDescription(){
		return "WeibullDistribution"+(inverse?"Inverse":"");
	}

	public Map<String, Double> MStep(List<Double> expectation){
		Map<String, Double> parameters=new HashMap<String, Double>();

		double pi=EMUtils.getPi(expectation);
		parameters.put("pi", pi);
		
		double ZLogX=getZLogX(sample, expectation);
		double Z=getZ(expectation);
		
		double alpha2=parameters.containsKey("alpha")?parameters.get("alpha"):1;
		double alpha=Double.POSITIVE_INFINITY;
		int it=0;
		while(it<2000&&Math.abs(alpha-alpha2)>1E-7){
			alpha=alpha2;
			double ZXAlphaLogX=getZXAlphaLogX(sample, expectation, alpha);
			double ZXAlpha=getZXAlpha(sample, expectation, alpha);
			double ZXAlphaLogX2=getZXAlphaLogX2(sample, expectation, alpha);
			double fn=ZLogX/Z+1/alpha-ZXAlphaLogX/ZXAlpha;
			double dfn=(1+Math.pow(alpha,  2))+(ZXAlpha*ZXAlphaLogX2-Math.pow(ZXAlphaLogX, 2))/Math.pow(ZXAlpha, 2);
			alpha2=alpha+fn/dfn;
			it++;
		}
		alpha=alpha2;
		
		parameters.put("alpha", alpha);
		
		double ZXAlpha=getZXAlpha(sample, expectation, alpha);
		double lambda=Math.pow(ZXAlpha/Z,1/alpha);
		parameters.put("lambda", lambda);
		
		return parameters;
	}
	
	public static double getZXAlpha(List<Double> sample, List<Double> expectation, double alpha){
		double mu=0;
		for(int j=0;j<expectation.size();j++){
			double e=expectation.get(j);
			double x=sample.get(j);
			mu+=e*Math.pow(x, alpha);
		}
		return mu;
	}	
	
	public static double getZXAlphaLogX2(List<Double> sample, List<Double> expectation, double alpha){
		double mu=0;
		for(int j=0;j<expectation.size();j++){
			double e=expectation.get(j);
			double x=sample.get(j);
			mu+=e*Math.pow(x, alpha)*Math.pow(Math.log(x),2);
		}
		return mu;
	}	
	
	public static double getZXAlphaLogX(List<Double> sample, List<Double> expectation, double alpha){
		double mu=0;
		for(int j=0;j<expectation.size();j++){
			double e=expectation.get(j);
			double x=sample.get(j);
			mu+=e*Math.pow(x, alpha)*Math.log(x);
		}
		return mu;
	}	
	
	public static double getZLogX(List<Double> sample, List<Double> expectation){
		double mu=0;
		for(int j=0;j<expectation.size();j++){
			double e=expectation.get(j);
			double x=sample.get(j);
			mu+=e*Math.log(x);
		}
		return mu;
	}	
	
	public static double getZ(List<Double> expectation){
		double mu=0;
		for(int j=0;j<expectation.size();j++){
			double e=expectation.get(j);
			mu+=e;
		}
		return mu;
	}

	@Override
	public double getNumericalMean(Map<String, Double> parameters) throws Exception {
		return new WeibullDistribution(parameters.get("alpha"), parameters.get("lambda")).getNumericalMean();
	}
	
	@Override
	public double getDensityUnprocessed(double x, Map<String, Double> parameters) throws Exception {
		return new WeibullDistribution(parameters.get("alpha"), parameters.get("lambda")).density(x);
	}
	
	@Override
	public double getCumulativeProbabilityUnprocessed(double x1, Map<String, Double> parameters) throws Exception {
		return new WeibullDistribution(parameters.get("alpha"), parameters.get("lambda")).cumulativeProbability(x1);
	}
}
