package de.unijena.bioinf.spectralcomparison;

import java.util.ArrayList;
import java.util.List;

import de.unijena.bioinf.decoy.model.MassBank;
import de.unijena.bioinf.deocy.Utils;
import de.unijena.bioinf.spectralcomparison.ParametersSpectralComparison.COMPARISONPOSTPROCESS;
import de.unijena.bioinf.statistics.Result;

public class ComparisonMethodEquality implements ComparisonMethod{

	static public final String methodNameLong="Equality";
	ParametersSpectralComparison p;
	public List<MassBank> s1;
	public List<MassBank> s2;


	public ComparisonMethodEquality(ParametersSpectralComparison p) {
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

		Result resultNotEqual=new Result(0, 0);
		Result resultEqual=new Result(1, 1);

		List<double[]> p1=new ArrayList<double[]>();		
		for(double[] p:peaksLeft){
				p1.add(new double[]{p[0],p[1]});
		}
		List<double[]> p2=new ArrayList<double[]>();
		for(double[] p:peaksRight){
				p2.add(new double[]{p[0],p[1]});
		}
		
		if(p1.size()!=p2.size())return resultNotEqual;
		
		for(double[] peaks1:p1){
			boolean found=false;
			for(double[] peaks2:p2){
				if(peaks1[0]==peaks2[0]&&peaks1[1]==peaks2[1]){
					found=true;
					break;
				}
			}
			if(!found)return resultNotEqual;
		}		
		
		return resultEqual;
	}

//	public double[][] calculateFP(double [][] scores) {
//
//		double[] means = new double[scores.length];
//		for (int i=0; i<scores.length; i++){
//			double sum =0;
//			int n=0;
//			for (int j=0; j<scores[i].length; j++){                    	
//				if(!Double.isInfinite(scores[i][j])){
//					sum += scores[i][j];
//					n++;
//				}
//			}
//			means[i] = sum/n;
//		}
//
//		double[][] correlationCoefficient = new double[scores.length][scores.length];
//		for (int i=0; i<correlationCoefficient.length; i++){
//			for (int j=0; j<correlationCoefficient.length; j++){
//
//				double combinedSum =0;
//				double isum =0;
//				double jsum=0;
//
//				for (int k=0; k<correlationCoefficient[i].length; k++){
//					if(!Double.isInfinite(scores[i][k])&&!Double.isInfinite(scores[j][k])){
//						combinedSum += (scores[i][k]-means[i])*(scores[j][k]-means[j]);
//						isum += Math.pow(scores[i][k]-means[i],2);
//						jsum += Math.pow(scores[j][k]-means[j],2);
//					}
//				}
//
//				correlationCoefficient[i][j]=combinedSum/(Math.pow(isum,0.5)*Math.pow(jsum,0.5));
//			}
//		}
//		return correlationCoefficient;
//	}


}
