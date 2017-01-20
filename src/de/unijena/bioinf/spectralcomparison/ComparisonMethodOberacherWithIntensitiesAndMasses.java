package de.unijena.bioinf.spectralcomparison;

import java.util.ArrayList;
import java.util.List;

import de.unijena.bioinf.decoy.model.MassBank;
import de.unijena.bioinf.deocy.Utils;
import de.unijena.bioinf.spectralcomparison.ParametersSpectralComparison.COMPARISONPOSTPROCESS;
import de.unijena.bioinf.statistics.Result;

public class ComparisonMethodOberacherWithIntensitiesAndMasses implements ComparisonMethod{

	static public final String methodNameLong="OberacherWithIntensitiesAndMasses";
	ParametersSpectralComparison p;
	public List<MassBank> s1;
	public List<MassBank> s2;


	public ComparisonMethodOberacherWithIntensitiesAndMasses(ParametersSpectralComparison p) {
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
		double intensityDiff=0;
		double massDiff=0;
		double intensitySum1=0;
		double intensitySum2=0;

		double maxLeft=0;
		for(double[] p:peaksLeft)maxLeft=Math.max(maxLeft, p[1]);
		double maxRight=0;
		for(double[] p:peaksRight)maxRight=Math.max(maxRight, p[1]);

		List<double[]> p1=new ArrayList<double[]>();
		for(double[] p:peaksLeft)if(p[1]/maxLeft*100>=5)p1.add(new double[]{p[0], p[1]/maxLeft*100});
		List<double[]> p2=new ArrayList<double[]>();
		for(double[] p:peaksRight)if(p[1]/maxRight*100>=5)p2.add(new double[]{p[0], p[1]/maxRight*100});
		
		for(double[] p:p1)intensitySum1+=p[1];
		for(double[] p:p2)intensitySum2+=p[1];
		
		for(double[] peaks1:p1){
			double[] bestMatchingPeak=null;
			double error=Utils.getAbsoluteErrorForMass(peaks1[0], p.ppm, p.ae);
			for(double[] peaks2:p2){
				double diff=Math.abs(peaks1[0]-peaks2[0]);
				if(diff<error&&(bestMatchingPeak==null||Math.abs(peaks1[0]-bestMatchingPeak[0])>diff)){
					bestMatchingPeak=peaks2;
				}
			}
			if(bestMatchingPeak!=null){
				numberMatches++;
				intensityDiff+=Math.abs(peaks1[1]-bestMatchingPeak[1]);
				massDiff+=Math.abs(peaks1[0]-bestMatchingPeak[0]);
			}
		}

		Result r=new Result(numberMatches, numberMatches!=0?1.0*Math.pow(numberMatches,4)*Math.pow(intensitySum1+2*intensitySum2,1.25)/(Math.pow(p1.size()+2*p2.size(),2)+intensityDiff+massDiff):0);
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
