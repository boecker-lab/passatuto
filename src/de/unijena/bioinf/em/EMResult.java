package de.unijena.bioinf.em;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

import de.unijena.bioinf.em.Sample.Type;

class EMResult implements Comparable<EMResult>{
	List<Map<String, Double>> parameters;
	final List<List<Double>> expectations;
	List<InterfaceMLDistribution> distributions;
	final List<Double> values;
	double likelihood;
	public Double deviation_QValuesByEFDR=null;
	public Double deviation_QValuesByPEP=null;


	public EMResult(List<Double> values, List<List<Double>> expectations, List<InterfaceMLDistribution> distributions, List<Map<String, Double>> parameters, double likelihood){
		this.distributions=distributions;
		this.parameters=parameters;
		this.likelihood=likelihood;
		this.expectations=expectations;
		this.values=new ArrayList<Double>();
		this.values.addAll(values);
	}

	@Override
	public int compareTo(EMResult o) {
		return Double.compare(o.likelihood,this.likelihood);
	}

	public String toString(){
		String description="";
		description+=likelihood;
		description+="\t"+deviation_QValuesByPEP;
		description+="\t"+deviation_QValuesByEFDR;

		TreeMap<Double, List<Integer>> means=new TreeMap<Double,List<Integer>>();
		for(int i=0;i<distributions.size();i++){
			InterfaceMLDistribution d=distributions.get(i);
			Map<String, Double> p=parameters.get(i);
			Double mean=Double.NaN;
			try{
				mean=d.getNumericalMean(p);				
			}catch(Exception e){}
			if(!means.containsKey(mean))means.put(mean,new ArrayList<Integer>());
			means.get(mean).add(i);
		}

		for(List<Integer> iList:means.values()){
			for(int i:iList){
				InterfaceMLDistribution d=distributions.get(i);
				Map<String, Double> p=parameters.get(i);
				description+="\t"+d.getDescription()+":";
				for(Entry<String, Double> e:p.entrySet()){
					description+=" "+e.getKey()+"="+e.getValue();
				}			
			}
		}
		return description;
	}

	public List<Double> getEFDRs() throws Exception{
		List<Double> fdrs=new ArrayList<Double>();
		boolean[] FPs = getFPIndizes();

		List<List<Double>> sortedExpectations=getSortedExpectationsByValue(values, this.expectations);
		double sum=0;
		double n=0;
		for(int j=0;j<sortedExpectations.get(0).size();j++){
			for(int i=0;i<FPs.length;i++){
				boolean FP=FPs[i];	
				double e=sortedExpectations.get(i).get(j);
				if(FP)sum+=e;
				n+=e;
			}
			fdrs.add(sum/n);
		}
		return fdrs;
	}

	public List<Double> getFDRsByPEP() throws Exception{
		List<Double> fdrs=new ArrayList<Double>();
		boolean[] FPs = getFPIndizes();
		List<Double> sortedValues=new ArrayList<Double>();
		sortedValues.addAll(values);

		Collections.sort(sortedValues,Collections.reverseOrder());

		for(int j=0;j<sortedValues.size();j++){
			double fpProb=0;
			double tpProb=0;
			double v=sortedValues.get(j);
			for(int i=0;i<FPs.length;i++){
				boolean FP=FPs[i];
				double e=distributions.get(i).getProbability(v, parameters.get(i));
				e*=parameters.get(i).get("pi");
				if(FP)fpProb+=e;
				else tpProb+=e;	
			}
			fdrs.add(fpProb/(fpProb+tpProb));
		}

		return fdrs;
	}
	
	public Map<Type,List<Double>> getPValues(Sample sample) throws Exception{
		Map<Type,List<Double>> result=new HashMap<Type,List<Double>>();
		boolean[] FPs = getFPIndizes();

		for(int j=0;j<sample.sample.size();j++){
			double fpSum=0;
			double v=sample.sample.get(j);
			Type t=sample.types.get(j);
			for(int i=0;i<FPs.length;i++){
				boolean FP=FPs[i];
				if(!FP)continue;
				double e=distributions.get(i).getProbability(v, parameters.get(i));
				e/=distributions.get(i).getProbability(Double.NEGATIVE_INFINITY, parameters.get(i));
				fpSum+=e;
			}
			if(!result.containsKey(t))result.put(t,new ArrayList<Double>());
			result.get(t).add(fpSum);
		}
		return result;
	}

	public List<Double> getQValuesByEDFR() throws Exception{
		List<Double> fdrs=getEFDRs();
		for(int i=fdrs.size()-2;i>=0;i--){
			fdrs.set(i,Math.min(fdrs.get(i+1), fdrs.get(i)));
		}
		return fdrs;
	}

	public List<Double> getQValuesByPEP() throws Exception{
		List<Double> fdrs=getFDRsByPEP();
		for(int i=fdrs.size()-2;i>=0;i--){
			fdrs.set(i,Math.min(fdrs.get(i+1), fdrs.get(i)));
		}
		return fdrs;
	}

	public List<List<Double>> getSortedExpectationsByValue(List<Double> values, List<List<Double>> expectations){
		List<Double> values_tmp=new ArrayList<Double>();
		values_tmp.addAll(values);
		List<List<Double>> expectations_tmp=new ArrayList<List<Double>>();
		for(List<Double> e:expectations){
			List<Double> e_tmp=new ArrayList<Double>();
			e_tmp.addAll(e);
			expectations_tmp.add(e_tmp);
		}

		for(int i=0;i<values_tmp.size()-1;i++){
			double maxValue=Integer.MIN_VALUE;
			int maxIndex=Integer.MAX_VALUE;
			for(int j=i;j<values_tmp.size();j++){
				if(values_tmp.get(j)>maxValue){
					maxValue=values_tmp.get(j);
					maxIndex=j;
				}
			}
			change(i,maxIndex, values_tmp);
			for(List<Double> e:expectations_tmp){				
				change(i,maxIndex, e);
			}
		}
		return expectations_tmp;
	}

	public boolean[] getFPIndizes() throws Exception{
		double minMean=Double.POSITIVE_INFINITY;
		int minIndex=-1;
		double maxMean=Double.NEGATIVE_INFINITY;
		int maxIndex=-1;
		boolean[] FPs=new boolean[parameters.size()];
		for(int i=0;i<parameters.size();i++){
			double mean=distributions.get(i).getNumericalMean(parameters.get(i));			
			if(mean<0.5)FPs[i]=true;
			else FPs[i]=false;
			if(minMean>mean){
				minMean=mean;
				minIndex=i;
			}
			if(maxMean<mean){
				maxMean=mean;
				maxIndex=i;
			}
		}
		FPs[minIndex]=true;
		FPs[maxIndex]=false;
		return FPs;
	}

	public static void change(int i, int j, List list){
		Object o=list.get(i);
		list.set(i,list.get(j));
		list.set(j, o);
	}


}