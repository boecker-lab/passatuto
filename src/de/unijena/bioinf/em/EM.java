package de.unijena.bioinf.em;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;

import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.GumbelDistribution;
import org.apache.commons.math3.distribution.LogisticDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomGenerator;

public class EM {

	public final Sample sample;
	public List<InterfaceMLDistribution> distributions;
	int maxIterations=10000;
	double minChangeLikelihood=1E-20;

	public EM(Sample sample, List<InterfaceMLDistribution> distributions){
		for(int i=0;i<distributions.size();i++)distributions.get(i).preprocessSample(sample.sample);
		this.sample=sample;
		this.distributions=distributions;		
	}

	public EM(Sample sample, InterfaceMLDistribution distribution){						
		this(sample, Arrays.asList(new InterfaceMLDistribution[]{distribution}));
	}

	public EM(Sample sample, List<InterfaceMLDistribution> distributions, int maxIterations, double minChangeLikelihood){
		this(sample, distributions);
		this.maxIterations=maxIterations;
		this.minChangeLikelihood=minChangeLikelihood;
	}

	public EM(Sample sample, InterfaceMLDistribution distribution, int maxIterations, double minChangeLikelihood){
		this(sample, distribution);
		this.maxIterations=maxIterations;
		this.minChangeLikelihood=minChangeLikelihood;
	}

	public EMResult doEM(boolean verb) throws Exception{
		List<List<Double>> expectations=EMUtils.initializeRandomExpectations(sample.sample.size(), distributions.size());

		List<Map<String, Double>> parameters=null;

		int iteration=0;
		double likelihood=Double.NEGATIVE_INFINITY;
		double changeLikelihood=Double.POSITIVE_INFINITY;

		while(iteration++<maxIterations&&changeLikelihood>minChangeLikelihood){		
			List<Map<String, Double>> parameters_new=EMUtils.MStep(distributions, expectations);

			if(verb){
				System.out.println("iteration "+iteration);
				for(Map<String, Double> p:parameters_new){
					for(Entry<String, Double> e:p.entrySet()){
						System.out.print(e.getKey()+": "+e.getValue()+"\t");
					}
					System.out.println();
				}
			}

			List<List<Double>> expectations_new=EMUtils.EStep(distributions, parameters_new);

			double currLikelihood=EMUtils.getLogLikelihood(distributions, parameters_new);
			changeLikelihood=currLikelihood-likelihood;
			if(verb){
				System.out.println("LogLik "+currLikelihood);
				System.out.println("Delta LogLik "+changeLikelihood+"\n");
			}
			if(changeLikelihood<0)
				System.out.println("delta likelihood decreased?!" );
			likelihood=currLikelihood;

			expectations=expectations_new;
			parameters=parameters_new;
		}

		parameters=EMUtils.MStep(distributions, expectations);

		EMResult result=new EMResult(sample.sample, expectations, distributions, parameters, likelihood);
		return result;

	}

	private static void getStandardOutput(Process proc) throws IOException {
		BufferedReader stdInput = new BufferedReader(new 
				InputStreamReader(proc.getInputStream()));

		BufferedReader stdError = new BufferedReader(new 
				InputStreamReader(proc.getErrorStream()));

		// read the output from the command
		System.out.println("Here is the standard output of the command:\n");
		String s = null;
		while ((s = stdInput.readLine()) != null) {
			System.out.println(s);
		}

		// read any errors from the attempted command
		System.out.println("Here is the standard error of the command (if any):\n");
		while ((s = stdError.readLine()) != null) {
			System.out.println(s);
		}
	}

	private boolean checkExpectations(List<List<Double>> expectations){
		for(int i=0;i<expectations.get(0).size();i++){
			double sum=0;
			for(int j=0;j<expectations.size();j++){
				sum+=expectations.get(j).get(i);
			}
			if(Math.abs(1-sum)>1E-7)return false;
		}
		return true;	
	}

	private boolean checkParameters(List<Map<String, Double>> parameters){
		double sum=0;
		for(Map<String, Double> p:parameters){
			sum+=p.get("pi");				
		}
		if(Math.abs(1-sum)>1E-7)return false;
		return true;


	}
}
