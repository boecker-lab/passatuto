package de.unijena.bioinf.spectralcomparison;

import java.io.File;
import java.util.Arrays;

import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction;

public class ParametersSpectralComparison {
	
	public static enum COMPARISONMETHOD{OBERACHER, OBERACHERWITHINTENSITIES, OBERACHERWITHINTENSITIESANDMASSES, OBERACHERPEAKCOUNTING, MASSBANK, MASSBANKEFFECTIVEPEAKS, COSINEDISTANCE, TREEALIGNMENT, EQUALITY};
	public static enum COMPARISONPOSTPROCESS{FP, SIMPLE};
	
	static String sep=System.getProperty("file.separator");
	
	String base;
	String baseRT="RunningTimes";
	ParametersDecoySpectrumConstruction p1;
	ParametersDecoySpectrumConstruction p2;
	COMPARISONMETHOD m;
	COMPARISONPOSTPROCESS pp;
	String addFolder="";
	String addFolderRunningTimes="";
	public boolean rt=false;
	public boolean removePP=false;
	public double ppm;
	public double ae;
	
	public ParametersSpectralComparison(COMPARISONMETHOD m,  COMPARISONPOSTPROCESS pp, double ppm, double ae, String addFolder, ParametersDecoySpectrumConstruction p1, ParametersDecoySpectrumConstruction p2){
		this(m, pp, ppm, ae, p1, p2);
		this.addFolder=addFolder;
	}
	
	public ParametersSpectralComparison(COMPARISONMETHOD m, COMPARISONPOSTPROCESS pp, double ppm, double ae, ParametersDecoySpectrumConstruction p1, ParametersDecoySpectrumConstruction p2){
		this.p1=p1;
		this.p2=p2;
		this.m=m;
		this.pp=pp;
		this.ppm=ppm;
		this.ae=ae;	
	}
	
	public ParametersSpectralComparison(String[] argsAll){
		int indexOfFirstDataset=Integer.MIN_VALUE;
		int indexOfSecondDataset=Integer.MIN_VALUE;
		for(int j=0;j<argsAll.length;j++){
			if(argsAll[j].equals("-ds1"))indexOfFirstDataset=j;
			if(argsAll[j].equals("-ds2"))indexOfSecondDataset=j;
		}
		String args[]=Arrays.copyOfRange(argsAll, 0, indexOfFirstDataset);
		int i=0;
		while(i<args.length){
			if(args[i].equals("-base"))this.base=args[++i];
			if(args[i].equals("-method")){
				String method=args[++i];
				for(COMPARISONMETHOD m:COMPARISONMETHOD.values()){
					if(method.equals(getComparisonMethodName(m)))this.m=m;
				}
			}
			if(args[i].equals("-addFolder"))this.addFolder=args[++i];
			if(args[i].equals("-rt")){
				this.rt=true;
				if(!args[i+1].startsWith("-"))addFolderRunningTimes=args[++i];
			}
			if(args[i].equals("-rpp")){
				this.removePP=true;
			}
			if(args[i].equals("-ppm"))this.ppm=Double.parseDouble(args[++i]);
			if(args[i].equals("-ae"))this.ae=Double.parseDouble(args[++i]);
			if(args[i].equals("-postprocess")){				
					String post=args[++i];
					if(post.equals("FP"))this.pp=COMPARISONPOSTPROCESS.FP;
					if(post.equals("Simple"))this.pp=COMPARISONPOSTPROCESS.SIMPLE;
			}
			i++;
		}
		p1=new ParametersDecoySpectrumConstruction(Arrays.copyOfRange(argsAll, indexOfFirstDataset+1, indexOfSecondDataset));
		p2=new ParametersDecoySpectrumConstruction(Arrays.copyOfRange(argsAll, indexOfSecondDataset+1, argsAll.length));
	}
	
	public ComparisonMethod getComparisonMethod(){
		switch(m){
			case OBERACHER: return new ComparisonMethodOberacher(this);
			case OBERACHERWITHINTENSITIES: return new ComparisonMethodOberacherWithIntensities(this);
			case OBERACHERWITHINTENSITIESANDMASSES: return new ComparisonMethodOberacherWithIntensitiesAndMasses(this);
			case OBERACHERPEAKCOUNTING: return new ComparisonMethodOberacherPeakcounting(this);
			case MASSBANK: return new ComparisonMethodMassBank(this);
			case MASSBANKEFFECTIVEPEAKS: return new ComparisonMethodMassBankEffectivePeaks(this);
			case COSINEDISTANCE: return new ComparisonMethodCosineDistance(this);
			case TREEALIGNMENT: return new ComparisonMethodTreeAlignment(this);
			case EQUALITY: return new ComparisonMethodEquality(this);
			default: return null;
		}
	}
	
	public static String getComparisonMethodName(COMPARISONMETHOD m){
		switch(m){
			case OBERACHER: return ComparisonMethodOberacher.methodNameLong;
			case OBERACHERWITHINTENSITIES: return ComparisonMethodOberacherWithIntensities.methodNameLong;
			case OBERACHERWITHINTENSITIESANDMASSES: return ComparisonMethodOberacherWithIntensitiesAndMasses.methodNameLong;
			case OBERACHERPEAKCOUNTING: return ComparisonMethodOberacherPeakcounting.methodNameLong;
			case MASSBANK: return ComparisonMethodMassBank.methodNameLong;
			case MASSBANKEFFECTIVEPEAKS: return ComparisonMethodMassBankEffectivePeaks.methodNameLong;
			case COSINEDISTANCE: return ComparisonMethodCosineDistance.methodNameLong;
			case TREEALIGNMENT: return ComparisonMethodTreeAlignment.methodNameLong;
			case EQUALITY: return ComparisonMethodEquality.methodNameLong;
			default: return null;
		}
	}
	
	public String getSearchMethod(){//TODO: SearchMethodClass
		return getComparisonMethod().getMethodNameLong();
	}
	
	public PostProcessMethod getPostProcess(){//TODO: SearchMethodClass
		switch(pp){
			case FP: return new PostProcessMethodFingerprint(this);
			case SIMPLE: return new PostProcessMethodSimple(this);
			default: return null;
		}
	}
	
	public File getOutputFileComparison(){
		return new File(base+sep+getSearchMethod()+sep+p1.getChargePath()+sep+p1.getMassbankFirstID()+"-"+p2.getMassbankFirstID()+(!addFolder.isEmpty()?(sep+p1.getMassbankFirstID()+p1.addFolder+"-"+p2.getMassbankFirstID()+addFolder):"")+getPostProcess().getMethodNameLong()+".csv");
	}

	public boolean isCorrectCombination() {
		// TODO Auto-generated method stub
		return true;
	}
	
	public File getOutputFileRunningTime(){
		return new File(base+sep+baseRT+sep+getSearchMethod()+sep+p1.getChargePath()+sep+p1.getMassbankFirstID()+"-"+p2.getMassbankFirstID()+(!addFolder.isEmpty()?(sep+p1.getMassbankFirstID()+p1.addFolder+"-"+p2.getMassbankFirstID()+addFolder):"")+getPostProcess().getMethodNameLong()+addFolderRunningTimes+".txt");
	}
	
}
