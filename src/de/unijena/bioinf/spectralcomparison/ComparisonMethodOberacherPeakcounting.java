package de.unijena.bioinf.spectralcomparison;

import java.util.ArrayList;
import java.util.List;

import de.unijena.bioinf.decoy.model.MassBank;
import de.unijena.bioinf.deocy.Utils;
import de.unijena.bioinf.spectralcomparison.ParametersSpectralComparison.COMPARISONPOSTPROCESS;
import de.unijena.bioinf.statistics.Result;

public class ComparisonMethodOberacherPeakcounting implements ComparisonMethod{

	static public final String methodNameLong="OberacherPeakcounting";
	ParametersSpectralComparison p;
	public List<MassBank> s1;
	public List<MassBank> s2;


	public ComparisonMethodOberacherPeakcounting(ParametersSpectralComparison p) {
		this.p=p;
	}

	public String getMethodNameLong(){
		return methodNameLong;
	}

	public void setSpectra(List<MassBank> s1, List<MassBank> s2) {
		this.s1=s1;
		this.s2=s2;
	}
	
	public List<MassBank>[] getSpectra() {
		return new List[]{s1,s2};
	}

	@Override
	public AlignmentMatrix getAlignmentMatrix() {
		Result[][] results=p.getPostProcess().getResultList(this);
		List<String> left=new ArrayList<String>();
		for(MassBank m:s1)left.add(m.massbankID);
		List<String> right=new ArrayList<String>();
		for(MassBank m:s2)right.add(m.massbankID);
		String queryFile=p.p1.getOutputFolderMassbank().getAbsolutePath();
		String resultFile=p.p2.getOutputFolderMassbank().getAbsolutePath();
		AlignmentMatrix matrix=new AlignmentMatrix(queryFile, resultFile, left, right, results);
		return matrix;
	}


	public Result getResult(List<double[]> peaksLeft, List<double[]> peaksRight){

		int numberMatches=0;
		for(double[] peaks1:peaksLeft){
			double[] bestMatchingPeak=null;
			double error=Utils.getAbsoluteErrorForMass(peaks1[0], p.ppm, p.ae);
			for(double[] peaks2:peaksRight){
				double diff=Math.abs(peaks1[0]-peaks2[0]);
				if(diff<error&&(bestMatchingPeak==null||Math.abs(peaks1[0]-bestMatchingPeak[0])>diff)){
					bestMatchingPeak=peaks2;
				}
			}
			if(bestMatchingPeak!=null){
				numberMatches++;
			}
		}

		Result r=new Result(numberMatches, 2.0*numberMatches/(peaksLeft.size()+peaksRight.size()));
		return r;
	}

	public double[][] calculateFP(double [][] scores) {

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

		double[][] correlationCoefficient = new double[scores.length][scores.length];
		for (int i=0; i<correlationCoefficient.length; i++){
			for (int j=0; j<correlationCoefficient.length; j++){

				double combinedSum =0;
				double isum =0;
				double jsum=0;

				for (int k=0; k<correlationCoefficient[i].length; k++){
					if(!Double.isInfinite(scores[i][k])&&!Double.isInfinite(scores[j][k])){
						combinedSum += (scores[i][k]-means[i])*(scores[j][k]-means[j]);
						isum += Math.pow(scores[i][k]-means[i],2);
						jsum += Math.pow(scores[j][k]-means[j],2);
					}
				}

				correlationCoefficient[i][j]=combinedSum/(Math.pow(isum,0.5)*Math.pow(jsum,0.5));
			}
		}
		return correlationCoefficient;
	}


}
