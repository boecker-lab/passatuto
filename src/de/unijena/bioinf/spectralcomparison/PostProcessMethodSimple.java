package de.unijena.bioinf.spectralcomparison;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import de.unijena.bioinf.decoy.model.MassBank;
import de.unijena.bioinf.statistics.Result;

public class PostProcessMethodSimple implements PostProcessMethod{
	static public final String methodNameLong="";
	ParametersSpectralComparison p;
	
	public String getMethodNameLong(){
		return methodNameLong;
	}
	
	public PostProcessMethodSimple(ParametersSpectralComparison p){
		this.p=p;
	}

	@Override
	public Result[][] getResultList(ComparisonMethod cm) {
		BufferedWriter bw=null;
		if(p.rt){
			File f=p.getOutputFileRunningTime();
			if(!f.getParentFile().exists())f.getParentFile().mkdirs();
			try{
				bw=new BufferedWriter(new FileWriter(f));
			}catch(IOException e){
				System.err.println("could not create output file for running times: "+f.getAbsolutePath());
			}
		}
		
		List<MassBank>[] spectra=cm.getSpectra();
		Result[][] results=new Result[spectra[0].size()][spectra[1].size()];
		long allTimes=0;
		for(int i=0;i<results.length;i++){
			long current=System.currentTimeMillis();
			for(int j=0;j<results[i].length;j++){
				results[i][j]=cm.getResult(spectra[0].get(i).peaks, spectra[1].get(j).peaks);			
			}
			long time=System.currentTimeMillis()-current;
			allTimes+=time;
			if(bw!=null){
				try{
					bw.write("time for processing "+spectra[0].get(i).massbankID+" (in milliseconds): "+time);
					bw.newLine();
				}catch(IOException e){
					System.err.println("could not write to output file for running times");
				}	
			}
		}
		if(bw!=null){
			try{
				bw.write("time for processing all compounds (in milliseconds): "+ allTimes);
				bw.newLine();
				bw.close();				
			}catch(IOException e){
				System.err.println("could not close output file for running times");
			}
		}
		return results;
	}
}
