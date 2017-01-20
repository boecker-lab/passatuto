package de.unijena.bioinf.em;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Random;
import java.util.TreeSet;

import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.RealDistribution;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.xy.IntervalXYDataset;
import org.jfree.data.xy.XYDataset;

import de.unijena.bioinf.em.Sample.Type;

public class EMUtils {

	public static void normalizeExpectations(List<List<Double>> expectations){
		for(int i=0;i<expectations.get(0).size();i++){
			double sum=0;
			for(int j=0;j<expectations.size();j++){
				sum+=expectations.get(j).get(i);
			}
			for(int j=0;j<expectations.size();j++){
				if(sum==0)expectations.get(j).set(i,0d);
				else expectations.get(j).set(i,expectations.get(j).get(i)/sum);
			}
		}			
	}

	public static double getLogLikelihood(List<InterfaceMLDistribution> distributions, List<Map<String, Double>> parameters) throws Exception{
		List<List<Double>> expectations=EStepWithoutNormalization(distributions,parameters);
		double sum=0;		
		for(int i=0;i<expectations.get(0).size();i++){			
			double sumTmp=0;
			for(int j=0;j<expectations.size();j++){
				sumTmp+=expectations.get(j).get(i);
			}
			if(sumTmp!=0)sum+=Math.log(sumTmp);
		}
		return sum;
	}

	public static List<Map<String, Double>> MStep(List<InterfaceMLDistribution> distributions, List<List<Double>> expectations){
		List<Map<String, Double>> parameters=new ArrayList<Map<String, Double>>();		
		for(int i=0;i<expectations.size();i++){
			parameters.add(distributions.get(i).MStep(expectations.get(i)));
		}
		
		return parameters;
	}
	
	private static List<List<Double>> EStepWithoutNormalization(List<InterfaceMLDistribution> distributions, List<Map<String,Double>> parameters) throws Exception{
		List<List<Double>> expectations=new ArrayList<List<Double>>();
		for(int i=0;i<parameters.size();i++){
			expectations.add(distributions.get(i).EStep(parameters.get(i)));
		}
		return expectations;
	}

	public static List<List<Double>> EStep(List<InterfaceMLDistribution> distributions, List<Map<String,Double>> parameters) throws Exception{
		List<List<Double>> expectations=EStepWithoutNormalization(distributions,parameters);
		normalizeExpectations(expectations);
		return expectations;
	}

	public static List<List<Double>> initializeRandomExpectations(int numberSamples, int numberDistributions){
		List<List<Double>> expectations=new ArrayList<List<Double>>();
		Random r=new Random();
		for(int n=0;n<numberDistributions;n++){
			List<Double> expectation=new ArrayList<Double>();
			expectations.add(expectation);
			for(int i=0;i<numberSamples;i++){
				expectation.add(r.nextDouble());
			}
		}
		normalizeExpectations(expectations);
		return expectations;
	}

	public static double getPi(List<Double> expectation){
		double piDividend=0;
		double piDivisor=0;
		for(double e:expectation){
			piDividend+=e;
			piDivisor++;
		}
		double pi=piDividend/piDivisor;
		return pi;
	}
	
	static public double[] getArrayOfValues(List<Double> sample){
		double[] tmp=new double[sample.size()];
		for(int i=0;i<tmp.length;i++){
			tmp[i]=sample.get(i);
		}
		return tmp;
	}

	public static double getMu(List<Double> sample, List<Double> expectation){
		double muDividend=0;
		double muDivisor=0;
		for(int j=0;j<expectation.size();j++){
			double e=expectation.get(j);
			double x=sample.get(j);
			muDividend+=e*x;
			muDivisor+=e;
		}
		double mu=muDividend/muDivisor;
		return mu;
	}			
	
	
	public static double deviation(List<Double> d1, List<Double> d2){
		double sum=0;
		for(int i=0;i<d1.size();i++){
			sum+=Math.abs(Math.pow(d1.get(i)-d2.get(i),2));
		}
		return sum/d1.size();
	}
	
	public static double getLogMu(List<Double> sample, List<Double> expectation){
		double muDividend=0;
		double muDivisor=0;
		for(int j=0;j<expectation.size();j++){
			double e=expectation.get(j);
			double x=sample.get(j);
			muDividend+=e*Math.log(x);
			muDivisor+=e;
		}
		double mu=muDividend/muDivisor;
		return mu;
	}	

	public static double getSD(List<Double> sample, List<Double> expectation, Double mu){
		double sdDividend=0;
		double sdDivisor=0;
		if(mu==null)mu=getMu(sample, expectation);

		for(int j=0;j<expectation.size();j++){
			double e=expectation.get(j);
			double x=sample.get(j);
			sdDividend+=e*Math.pow(x-mu,2);
			sdDivisor+=e;
		}
		double sd=Math.pow(sdDividend/sdDivisor,0.5);
		return sd;
	}
	
	public static double getLogSD(List<Double> sample, List<Double> expectation, Double mu){
		double sdDividend=0;
		double sdDivisor=0;
		if(mu==null)mu=getLogMu(sample, expectation);

		for(int j=0;j<expectation.size();j++){
			double e=expectation.get(j);
			double x=sample.get(j);
			sdDividend+=e*Math.pow(Math.log(x)-mu,2);
			sdDivisor+=e;
		}
		double sd=Math.pow(sdDividend/sdDivisor,0.5);
		return sd;
	}

	public static double getN(List<Double> expectation){
		double n=0;
		for(int j=0;j<expectation.size();j++){
			double e=expectation.get(j);			
			n+=e;
		}
		return n;
	}

	public static List<Sample> getScores(File inputFile, boolean log, boolean discardZeroEntries) throws Exception{
		Sample TPs=new Sample();
		Sample FPs=new Sample();
		List<Sample> all=new ArrayList<Sample>();
		all.add(FPs);
		all.add(TPs);
		BufferedReader br=new BufferedReader(new FileReader(inputFile));
		br.readLine();
		br.readLine();
		String line;
		int i=0;
		double minValue=Double.POSITIVE_INFINITY;
		while((line=br.readLine())!=null){
			String[] l=line.split("\t");
			String type=l[1];
			double score=Double.parseDouble(l[2]);
			if(discardZeroEntries&&score<=1E-5){
				continue;
			}
			if(log)score=-Math.log(1-score);
			if(type.equals("TruePositiveMatch"))
				TPs.add(score,type,i);
			if(type.equals("FalsePositiveMatch"))
				FPs.add(score,type,i);
			minValue=Math.min(minValue, score);
			i++;
		}
//		if(log){
//			for(int k=0;k<TPs.size();k++){
//				TPs.sample.set(k,TPs.sample.get(k)-minValue);
//			}
//			for(int k=0;k<FPs.size();k++){
//				FPs.sample.set(k,FPs.sample.get(k)-minValue);
//			}
//		}
		br.close();
		return all;
	}

	public static List<List<Double>> getSimulatedScores(int[] n, List<RealDistribution> distributions){
		List<List<Double>> sample_sep=new ArrayList<List<Double>>();
		for(int i=0;i<n.length;i++){
			RealDistribution nd=distributions.get(i);
			List<Double> sample_curr=new ArrayList<Double>();
			for(double d:nd.sample(n[i])){
				sample_curr.add(d);
			}
			sample_sep.add(sample_curr);
		}
		return sample_sep;
	}

	public static void writeData(File outputFile, List<Double> data, List<Type> types) throws Exception{
		if(!outputFile.getParentFile().exists())outputFile.getParentFile().mkdirs();
		BufferedWriter bw=new BufferedWriter(new FileWriter(outputFile));
		int i=0;
		for(Double d:data){
			bw.write(Double.toString(d));
			if(types!=null)bw.write("\t"+types.get(i++));
			bw.newLine();
		}
		bw.close();

	}

	public static void writeDistribution(File outputFile, List<Double> data, List<InterfaceMLDistribution> distributions, List<Map<String, Double>> parameters) throws Exception{
		if(!outputFile.getParentFile().exists())outputFile.getParentFile().mkdirs();
		BufferedWriter bw=new BufferedWriter(new FileWriter(outputFile));

		double minValue=Double.POSITIVE_INFINITY;
		double maxValue=Double.NEGATIVE_INFINITY;
		for(double d:data){
			minValue=Math.min(minValue, d);
			maxValue=Math.max(maxValue, d);
		}

		double step=(maxValue-minValue)/1000;
		for(double x=minValue;x<=maxValue;x+=step){
			bw.write(Double.toString(x));
			for(int i=0;i<distributions.size();i++){
				InterfaceMLDistribution d=distributions.get(i);
				Map<String, Double> p=parameters.get(i);
				Double v=p.get("pi")*d.getDensity(x, p);
				bw.write("\t"+v);				
			}
			bw.newLine();
		}
		bw.close();


	}

	public static void writeQValues(File outputFile, String description, List<Double> qValuesReal, List<Double> qValuesPred)throws Exception{
		if(!outputFile.getParentFile().exists())outputFile.getParentFile().mkdirs();
		BufferedWriter bw=new BufferedWriter(new FileWriter(outputFile));
		bw.write(description);
		bw.newLine();
		bw.write("calculated qValue\testimated qValue");
		bw.newLine();
		double sum=0;
		int n=0;
		for(int i=0;i<qValuesReal.size();i++){
			sum+=qValuesPred.get(i);
			n+=1;
			if(i==qValuesReal.size()-1||Double.compare(qValuesReal.get(i+1),qValuesReal.get(i))!=0){
				bw.write(qValuesReal.get(i)+"\t"+sum/n);
				bw.newLine();
				n=0;
				sum=0;
			}
		}


		bw.close();
	}
	
	public static void writePValues(File outputFile, String description, List<Double> pValuesPred, int steps)throws Exception{
		Locale.setDefault(Locale.US);
		DecimalFormat df=new DecimalFormat("0.000");
		int[] freq=new int[steps];
		int all=0;
		for(double d:pValuesPred){
			freq[(int)Math.floor(d*steps)]++;
			all++;
		}
		if(!outputFile.getParentFile().exists())outputFile.getParentFile().mkdirs();
		BufferedWriter bw=new BufferedWriter(new FileWriter(outputFile));
		bw.write(description);
		bw.newLine();
		bw.write("estimated pValue\tFrequency");
		bw.newLine();
		for(int i=0;i<freq.length;i++){
			bw.write(df.format(1.0*i/steps+0.5/steps)+"\t"+1.0*freq[i]/all);
			bw.newLine();			
		}
		bw.close();
	}
}
