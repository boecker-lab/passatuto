package de.unijena.bioinf.em;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public abstract class MLDistribution implements InterfaceMLDistribution{
	public double shift=+1E-5;
	static double borders=+1E-4;
	public boolean inverse=false;
	public Double inverse_x=null;
	public Double max_sample_value=null;
	public Double min_sample_value=null;
	public List<Double> sample_original=null;
	public List<Double> sample;

	public MLDistribution(boolean inverse, double shift){
		this.inverse=inverse;
		this.shift=shift;
	}

	@Override
	public List<Double> EStep(Map<String,Double> parameters) throws Exception{
		List<Double> expectation=new ArrayList<Double>();
		for(int i=0;i<sample.size();i++){
			double x=sample.get(i);
			double value=parameters.get("pi")*getDensityUnprocessed(x,parameters);
			expectation.add(value);
		}
		return expectation;
	}

	@Override
	public void preprocessSample(List<Double> sample_final){	
		sample_original=sample_final;

		this.sample=new ArrayList<Double>();
		List<Double> sample_tmp=new ArrayList<Double>();
		if(!inverse){
			sample_tmp.addAll(sample_final);
		}else{
			double max=Double.NEGATIVE_INFINITY;
			for(Double d:sample_final){
				max=Math.max(max, d);
			}
			inverse_x=max;
			for(Double d:sample_final){
				sample_tmp.add(max-d);
			}
		}

		max_sample_value=Double.NEGATIVE_INFINITY;
		min_sample_value=Double.POSITIVE_INFINITY;
		for(double s:sample_tmp){
			double s_tmp=s+shift;
			sample.add(s_tmp);
			max_sample_value=Math.max(max_sample_value, s_tmp);
			min_sample_value=Math.min(min_sample_value, s_tmp);
		}
		
	}

	@Override
	public double getDensity(double x, Map<String, Double> parameters) throws Exception {
		if(inverse) x=inverse_x-x;
		x+=shift;
		return getDensityUnprocessed(x,parameters);
	}

	@Override
	public double getProbability(double x, Map<String, Double> parameters) throws Exception {
		
		if(inverse) x=inverse_x-x;
		x+=shift;
		double y=max_sample_value+borders;
		if(inverse){
			y=x;
			x=min_sample_value-borders;
		}
		double result=getCumulativeProbabilityUnprocessed(y,parameters)-getCumulativeProbabilityUnprocessed(x,parameters);
		return result;
		
	}




}
