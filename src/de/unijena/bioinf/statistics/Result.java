package de.unijena.bioinf.statistics;

import de.unijena.bioinf.decoy.model.MassBank;

public class Result implements Comparable<Result>{
	String resultID;
	public int matchedPeaks;
	public double score;
	boolean isTrueMatch;
	MassBank massbank;
	
	public Result(String[] line){
		resultID=line[1];
		matchedPeaks=Integer.parseInt(line[3]);
		score=Double.parseDouble(line[4]);
	}
	
	public Result(int matchedPeaks, double score){
		this.matchedPeaks=matchedPeaks;
		this.score=score;
	}
	
	public Result(String resultID, int matchedPeaks, double score){
		this(matchedPeaks, score);
		this.resultID=resultID;
	}

	@Override
	public int compareTo(Result r) {
		return Double.compare(r.score,score);
	}
	
	public String getDB(){
		return resultID.replaceAll("\\d","");
	}
	
	public String getNumber(){
		return resultID.replaceAll("\\D","");
	}
	
	public String getDataset(){
		return resultID.substring(0,1);
	}
}