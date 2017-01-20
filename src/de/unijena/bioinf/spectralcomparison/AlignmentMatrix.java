package de.unijena.bioinf.spectralcomparison;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.List;

import de.unijena.bioinf.statistics.Result;

public class AlignmentMatrix {
	
	String queryFile;
	String resultFile;
	List<String> left;
	List<String> right;
	Result[][] results;
	
	public AlignmentMatrix(String queryFile, String resultFile, List<String> leftEntries, List<String> rightEntries, Result[][] results){
		this.queryFile=queryFile;
		this.resultFile=resultFile;
		this.left=leftEntries;
		this.right=rightEntries;
		this.results=results;
	}
	
	public void writeToCSV(File outputFile) throws Exception{
		if(!outputFile.getParentFile().exists())outputFile.getParentFile().mkdirs();
		BufferedWriter bw=new BufferedWriter(new FileWriter(outputFile));
		bw.write("Query: "+queryFile+";");
		bw.newLine();
		bw.write("DB: "+resultFile+";");
		bw.newLine();
		bw.newLine();
		
		for(String s:right)bw.write(","+s);
		bw.newLine();
		for(int i=0;i<left.size();i++){
			bw.write(left.get(i));
			for(Result r:results[i]){
				bw.write(","+r.score+" "+r.matchedPeaks);
			}
			bw.newLine();
		}
		bw.close();
	}

}
