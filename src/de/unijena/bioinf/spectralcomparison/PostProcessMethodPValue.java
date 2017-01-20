package de.unijena.bioinf.spectralcomparison;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.TreeMap;

import de.unijena.bioinf.babelms.dot.Graph;
import de.unijena.bioinf.decoy.model.MassBank;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction;
import de.unijena.bioinf.pvalue.SisterSpectrumGeneratorConditionalPeaks;
import de.unijena.bioinf.pvalue.SisterSpectrumGeneratorRandomPeaks;
import de.unijena.bioinf.statistics.Result;

public class PostProcessMethodPValue implements PostProcessMethod{
	static public final String methodNameLong="PValue";
	ParametersSpectralComparison p;
	Random r=new Random();
	
	SisterSpectrumGeneratorRandomPeaks g;
	TreeMap<Double, Double> oversampling;
	TreeMap<Double, Integer> resultingDistribution;
	int z;
	int N=100000;
	double step=0.1;
	
	public String getMethodNameLong(){
		return methodNameLong;
	}
	
	public PostProcessMethodPValue(ParametersSpectralComparison p){
		this.p=p;
	}

	@Override
	public Result[][] getResultList(ComparisonMethod cm) {
		List<MassBank>[] spectra=cm.getSpectra();
		Result[][] results=new Result[spectra[0].size()][spectra[1].size()];
		for(int i=0;i<results.length;i++){
			List<Double> tmp=new ArrayList<Double>();
			for(int j=0;j<results[i].length;j++){
				results[i][j]=cm.getResult(spectra[0].get(i).peaks, spectra[1].get(j).peaks);
				tmp.add(results[i][j].score);
			}
						
			if(spectra[0].get(i).massbankID.equals("MPOriginal001292")){
				Collections.sort(tmp, Collections.reverseOrder());
				for(double d:tmp)System.out.println(d);
				getPValues(spectra[0].get(i).peaks, p.p2, cm);
			}
						
		}
		return results;
	}
	
	public void getPValues(List<double[]> originalSpectrum,ParametersDecoySpectrumConstruction  p, ComparisonMethod cm){
		oversampling=new TreeMap<Double,Double>();
		resultingDistribution=new TreeMap<Double,Integer>();
		for(double i=0;i<=1;i+=step)oversampling.put(i,1.0);
		g=new SisterSpectrumGeneratorRandomPeaks(originalSpectrum.size(), p);
		List<double[]> startSpectrum=g.getStartSpectrum();
		System.out.println(cm.getResult(originalSpectrum, startSpectrum).score);
		try{
			BufferedWriter bwP=new BufferedWriter(new FileWriter(new File("U:\\MSBlast\\TrajectorySplittingP.txt")));
			BufferedWriter bwMy=new BufferedWriter(new FileWriter(new File("U:\\MSBlast\\TrajectorySplittingMy.txt")));
			BufferedWriter bwN=new BufferedWriter(new FileWriter(new File("U:\\MSBlast\\TrajectorySplittingN.txt")));
			for(int test=0;test<100;test++){
				z=0;
				double minOversampling=Double.MAX_VALUE;
				for(double d:oversampling.values())minOversampling=Math.min(minOversampling, d);
				simulateDPRTrajectory(startSpectrum, originalSpectrum, cm, minOversampling, 1);

				
				int n=0;
				double sum=0;
				
				for(Entry<Double, Double> e:oversampling.entrySet()){
					if(!resultingDistribution.containsKey(e.getKey()))resultingDistribution.put(e.getKey(), 0);
					resultingDistribution.put(e.getKey(), resultingDistribution.get(e.getKey())+1);
					bwN.write(e.getKey()+" "+resultingDistribution.get(e.getKey())+"\n");
				}
				
				for(Entry<Double, Integer> e:resultingDistribution.entrySet()){					
					n+=e.getValue();
					sum+=e.getValue()/oversampling.get(e.getKey());
				}
			
				for(double i:oversampling.keySet()){
					bwP.write(i+" "+resultingDistribution.get(i)/oversampling.get(i)/sum+"\n");
					bwMy.write(i+" "+sum*oversampling.get(i)/resultingDistribution.get(i)+"\n");
					oversampling.put(i,sum*oversampling.get(i)/resultingDistribution.get(i));
				}
				bwN.write(n+"\n");
				bwN.write("\n");
				bwN.flush();
				bwP.write(sum+"\n");
				bwP.write("\n");
				bwP.flush();
				bwMy.write("\n");
				bwMy.flush();
				resultingDistribution.clear();
				System.out.println(test+". done!");
			}
			bwN.close();
			bwP.close();
			bwMy.close();
		}catch(Exception e){
			System.err.println(e);
		}
	}
	
	public void simulateDPRTrajectory(List<double[]> decoySpectrum, List<double[]> originalSpectrum, ComparisonMethod cm, double omega, int level){
		while(z<N){
			List<double[]> sisterSpectrum=g.getSisterSpectrum(decoySpectrum);
			double oversamplingFactorOriginal=getOversamplingFactor(cm.getResult(originalSpectrum, decoySpectrum).score);
			double scoreSisterSpectrum=cm.getResult(originalSpectrum, sisterSpectrum).score;
			double oversamplingFactorSister=getOversamplingFactor(scoreSisterSpectrum);
			if(oversamplingFactorSister<omega)return;
			if(oversamplingFactorSister>oversamplingFactorOriginal){
				double y=oversamplingFactorSister/oversamplingFactorOriginal;
				double remaining=oversamplingFactorSister%oversamplingFactorOriginal;
				double rand=r.nextDouble();
				double yInt=Math.floor(y);
				if(rand<remaining)yInt=Math.floor(y)+1;
				for(long i=0;i<yInt;i++){
					rand=r.nextDouble();
					double newOmega=oversamplingFactorOriginal+rand*(oversamplingFactorSister-oversamplingFactorOriginal);
					simulateDPRTrajectory(sisterSpectrum, originalSpectrum, cm, newOmega, level++);
				}
			}
			double group=getOversamplingKey(scoreSisterSpectrum);
			if(!resultingDistribution.containsKey(group))resultingDistribution.put(group, 0);
			resultingDistribution.put(group, resultingDistribution.get(group)+1);
			System.out.println(z);
			z++;
		}
	}
	
	public double getOversamplingFactor(double score){
		double oversamplingKey=getOversamplingKey(score);
		return oversampling.get(oversamplingKey);
	}
	
	public double getOversamplingKey(double score){
		double overSamplingFactorCeiling=oversampling.ceilingKey(score);
		double overSamplingFactorFloor=oversampling.floorKey(score);
		return Math.abs(score-overSamplingFactorCeiling)<Math.abs(score-overSamplingFactorFloor)?overSamplingFactorCeiling:overSamplingFactorFloor;
	}
}
