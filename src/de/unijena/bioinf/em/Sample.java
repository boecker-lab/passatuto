package de.unijena.bioinf.em;

import java.util.ArrayList;
import java.util.List;

public class Sample {
	
	public enum Type{FalsePositiveMatch, TruePositiveMatch, DecoyMatch}
	
	public List<Double> sample=new ArrayList<Double>();
	public List<Type> types=new ArrayList<Type>();
	public List<Integer> ranks=new ArrayList<Integer>();
	
	public Sample(){
		
	}
	
	public Sample(Sample s){
		for(Double d:s.sample)this.sample.add(d);
		for(Type d:s.types)this.types.add(d);
		for(Integer d:s.ranks)this.ranks.add(d);
	}
	
	public Sample(List<Sample> sample){
		for(Sample s:sample)this.add(s);
	}
	
	public Sample(List<Double> sample, List<Type> types, List<Integer> ranks){
		this.sample=sample;
		this.types=types;
		this.ranks=ranks;
	}
	
	public void add(Double score, String type, Integer rank){		
		Type res=null;
		for(Type s:Type.values()){
			if(s.toString().equals(type))res=s;
		}
		add(score, res, rank);
		
	}
	
	public void add(Double score, Type type, int rank){
		sample.add(score);	
		types.add(type);
		ranks.add(rank);
	}

	public void add(Sample s) {
		for(int i=0;i<s.size();i++){
			add(s.sample.get(i), s.types.get(i), s.ranks.get(i));
		}
	}
	
	public int size(){
		return sample.size();
	}
	
	public Sample sortByRank(){
		Sample s=new Sample(this);
		for(int i=0;i<s.ranks.size()-1;i++){
			int minValue=Integer.MAX_VALUE;
			int minIndex=Integer.MAX_VALUE;
			for(int j=i;j<s.ranks.size();j++){
				if(s.ranks.get(j)<minValue){
					minValue=s.ranks.get(j);
					minIndex=j;
				}
			}
			change(i,minIndex, s.sample);
			change(i,minIndex, s.types);
			change(i,minIndex, s.ranks);
		}
		return s;
	}
	
	public static void change(int i, int j, List list){
		Object o=list.get(i);
		list.set(i,list.get(j));
		list.set(j, o);
	}
	
	public List<Double> getFDRs(){
		Sample s=sortByRank();
		List<Double> fdrs=new ArrayList<Double>();
		int numberFPs=0;;
		int numberTPs=0;
		for(int i=0;i<s.types.size();i++){
			if(s.types.get(i).equals(Type.FalsePositiveMatch))numberFPs++;
			if(s.types.get(i).equals(Type.TruePositiveMatch))numberTPs++;
			fdrs.add(1d*numberFPs/(numberFPs+numberTPs));
		}
		return fdrs;
	}
	
	public List<Double> getQValues(){
		List<Double> fdrs=getFDRs();
		for(int i=fdrs.size()-2;i>=0;i--){
			fdrs.set(i,Math.min(fdrs.get(i+1), fdrs.get(i)));
		}
		return fdrs;
	}
	
}
