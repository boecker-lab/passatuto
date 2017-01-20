package de.unijena.bioinf.spectralcomparison;

import java.util.ArrayList;
import java.util.List;

import de.unijena.bioinf.decoy.model.MassBank;
import de.unijena.bioinf.statistics.Result;

public class PostProcessMethodFingerprint implements PostProcessMethod{
	static public final String methodNameLong="FP";
	ParametersSpectralComparison p;
	
	public String getMethodNameLong(){
		return methodNameLong;
	}
	
	public PostProcessMethodFingerprint(ParametersSpectralComparison p){
		this.p=p;
	}

	@Override
	public Result[][] getResultList(ComparisonMethod cm) {
		List<MassBank>[] spectra=cm.getSpectra();
		Result[][] results=new Result[spectra[0].size()][spectra[1].size()];
		double[][] scoresQueryVsQuery=new double[spectra[0].size()][spectra[0].size()];
		double[][] scoresDBVsQuery=new double[spectra[1].size()][spectra[0].size()];
		 
		for(int j=0;j<spectra[0].size();j++){
			for(int i=0;i<spectra[0].size();i++){
				Result r=cm.getResult(spectra[0].get(i).peaks, spectra[0].get(j).peaks);
				scoresQueryVsQuery[i][j]=r.score;
			}
			for(int i=0;i<spectra[1].size();i++){
				Result r=cm.getResult(spectra[1].get(i).peaks, spectra[0].get(j).peaks);
				results[j][i]=r;
				scoresDBVsQuery[i][j]=r.score;
			}			
		}
		double[][] scoresFP=calculateFP(scoresQueryVsQuery, scoresDBVsQuery);
		for(int i=0;i<scoresFP.length;i++){
			for(int j=0; j<scoresFP[i].length;j++){
				results[i][j].score=scoresFP[i][j];
			}
		}
		return results;
	}
	
	public double[][] calculateFP(double [][] scoresQueryVsQuery, double[][] scoresDBVsQuery) {

		double[] meansQueryVsQuery = getMeans(scoresQueryVsQuery);
		double[] meansDBVsQuery = getMeans(scoresDBVsQuery);
		
		double[][] correlationCoefficient = new double[meansQueryVsQuery.length][meansDBVsQuery.length];
		
		for (int i=0; i<meansQueryVsQuery.length; i++){
			for (int j=0; j<meansDBVsQuery.length; j++){

				double combinedSum =0;
				double isum =0;
				double jsum=0;

				for (int k=0; k<scoresQueryVsQuery[i].length; k++){
					if(!Double.isInfinite(scoresQueryVsQuery[i][k])&&!Double.isInfinite(scoresDBVsQuery[j][k])){
						combinedSum += (scoresQueryVsQuery[i][k]-meansQueryVsQuery[i])*(scoresDBVsQuery[j][k]-meansDBVsQuery[j]);
						isum += Math.pow(scoresQueryVsQuery[i][k]-meansQueryVsQuery[i],2);
						jsum += Math.pow(scoresDBVsQuery[j][k]-meansDBVsQuery[j],2);
					}
				}

				correlationCoefficient[i][j]=combinedSum/(Math.pow(isum,0.5)*Math.pow(jsum,0.5));
			}
		}
		return correlationCoefficient;
	}
	
	public double[] getMeans(double[][] scores){
		double[] means = new double[scores.length];
		for (int i=0; i<scores.length; i++){
			double sum =0;
			int n=0;
			for (int j=0; j<scores[i].length; j++){                    	
				if(!Double.isInfinite(scores[i][j])){
					sum += scores[i][j];
					n++;
				}
			}
			means[i] = sum/n;
		}
		return means;
	}
}
