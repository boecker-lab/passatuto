package de.unijena.bioinf.msblast;

import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import de.unijena.bioinf.babelms.dot.Graph;
import de.unijena.bioinf.decoy.decoytrees.DefaultDecoyTree;
import de.unijena.bioinf.decoy.model.DecoyTreeVertex;
import de.unijena.bioinf.decoy.model.decoytreeconstructors.ConditionalFastConstructor;
import de.unijena.bioinf.decoy.model.decoytreeconstructors.ConditionalPeaksConstructor;
import de.unijena.bioinf.decoy.model.decoytreeconstructors.DefaultDecoyFileConstructor;
import de.unijena.bioinf.decoy.model.decoytreeconstructors.DefaultDecoyTreeConstructor;
import de.unijena.bioinf.decoy.model.decoytreeconstructors.MixedSpectrumConstructor;
import de.unijena.bioinf.decoy.model.decoytreeconstructors.RandomPeaksConstructor;
import de.unijena.bioinf.decoy.model.decoytreeconstructors.RandomTreeDecoyTreeConstructor;
import de.unijena.bioinf.decoy.model.decoytreeconstructors.RerootDecoyTreeConstructor;
import de.unijena.bioinf.deocy.Utils;

public class ParametersDecoySpectrumConstruction{
	public static enum CHARGE{POS,NEG};
	public static enum DATASET{AGILENT,METLIN,MASSBANKORBI,MASSBANKQTOF,MASSBANK,GNPS};
	public static enum DATA{FILES,TREES};
	public static enum METHOD{ORIGINAL, REROOT, RANDOMPEAKS, CONDITIONALPEAKS, CONDITIONALFAST, MIXEDSPECTRUM, RANDOMTREE};
	public static enum POSTPROCESS{PPDiff};
	
	static String sep=System.getProperty("file.separator");
	static String baseOriginalFiles="OriginalFiles";
	static String baseOriginalTrees="OriginalTrees";
	static String baseDB="DataBases";
	static String baseRT="RunningTimes";
	static String dot="dot";
	static String massbank="massbank";
	static String ms="ms";
		
	public String base;
	CHARGE c;
	DATASET ds;
	DATA d;
	METHOD m;
	List<POSTPROCESS> pp;
	public String addFolder="";
	public boolean rt;
	public double ppm;
	public double ae;
	
	Map<File,Graph> allGraphsInOriginal;
	Map<File,List<double[]>> allPeaksInOriginal;
	
	DecimalFormat df=new DecimalFormat("000000");
	
	
	public ParametersDecoySpectrumConstruction(String base, CHARGE c, DATASET ds, DATA d, METHOD m, List<POSTPROCESS> pp, double ppm, double ae){
		this.base=base;
		this.c=c;		
		this.ds=ds;
		this.d=d;
		this.m=m;
		this.pp=pp;
		this.ppm=ppm;
		this.ae=ae;
	}
	
	public ParametersDecoySpectrumConstruction(String[] args){
		this.pp=new ArrayList<POSTPROCESS>();
		int i=0;
		while(i<args.length){
			if(args[i].equals("-base"))this.base=args[++i];
			if(args[i].equals("-c")){
				String charge=args[++i];
				if(charge.equals("pos"))this.c=CHARGE.POS;
				if(charge.equals("neg"))this.c=CHARGE.NEG;
			}
			if(args[i].equals("-db")){
				String db=args[++i];
				if(db.equals("agilent"))this.ds=DATASET.AGILENT;
				if(db.equals("metlin"))this.ds=DATASET.METLIN;
				if(db.equals("massbankorbi"))this.ds=DATASET.MASSBANKORBI;
				if(db.equals("massbankqtof"))this.ds=DATASET.MASSBANKQTOF;
				if(db.equals("massbank"))this.ds=DATASET.MASSBANK;
				if(db.equals("gnps"))this.ds=DATASET.GNPS;
			}
			if(args[i].equals("-d")){
				String data=args[++i];
				for(DATA d:DATA.values()){
					if(data.equals(getDataName(d)))this.d=d;
				}
			}
			if(args[i].equals("-method")){
				String db=args[++i];
				for(METHOD m:METHOD.values()){
					if(db.equals(getDecoyTreeConstructorName(m)))this.m=m;
				}
			}
			if(args[i].equals("-addFolder"))this.addFolder=args[++i];
			if(args[i].equals("-rt"))this.rt=true;
			if(args[i].equals("-ppm"))this.ppm=Double.parseDouble(args[++i]);
			if(args[i].equals("-ae"))this.ae=Double.parseDouble(args[++i]);
			if(args[i].equals("-postprocess")){
				while(!args[i+1].startsWith("-")){
					String post=args[++i];
					if(post.equals("PPDiff"))this.pp.add(POSTPROCESS.PPDiff);
				}
			}
			i++;
		}
		
	}
	
	public boolean isCorrectCombination(){
		if(d.equals(DATA.FILES)&&(m.equals(METHOD.REROOT)||m.equals(METHOD.RANDOMTREE)))return false;
		return true;
	}
	
	public ParametersDecoySpectrumConstruction(String base, CHARGE c, DATASET ds, DATA d, METHOD m, List<POSTPROCESS> pp, String addFolder, double ppm, double ae){
		this(base, c, ds, d, m, pp, ppm, ae);
		this.addFolder=addFolder;
	}
	
	public Map<File,Graph> getAllGraphsInOriginal(){
		if (allGraphsInOriginal==null){
			allGraphsInOriginal=Utils.readGraphs(getFolderOriginalFiles());
		}
		return allGraphsInOriginal;		
	}
	
	public Map<File,List<double[]>> getAllPeaksInOriginal(){
		if (allPeaksInOriginal==null){
			allPeaksInOriginal=new HashMap<File,List<double[]>>();
			Map<File,Graph> allGraphs=getAllGraphsInOriginal();
			for(Entry<File,Graph> e:allGraphs.entrySet()){
				DefaultDecoyTree dt=new DefaultDecoyTree(e.getValue());
				List<double[]> r=new ArrayList<double[]>();
				for(DecoyTreeVertex dv:dt.getVertices()){
					r.add(new double[]{dv.mz,dv.intensity});
				}
				allPeaksInOriginal.put(e.getKey(), r);
			}
		}
		return allPeaksInOriginal;		
	}
	
	public void deleteAllGraphsInOriginalDB(){
		if(allGraphsInOriginal!=null)allGraphsInOriginal.clear();
	}
	
	public File getInchiKeysFile(){
		return new File(base+sep+getDatasetLong()+".csv");
	}
	
	public DefaultDecoyTreeConstructor getDecoyTreeConstructor(){
		switch(m){
//			case ORIGINALTREES: return new DefaultDecoyTreeConstructor(this);
			case ORIGINAL: return new DefaultDecoyTreeConstructor(this);
			case REROOT: return new RerootDecoyTreeConstructor(this);
			case RANDOMPEAKS: return new RandomPeaksConstructor(this);
			case CONDITIONALPEAKS: return new ConditionalPeaksConstructor(this);
			case CONDITIONALFAST: return new ConditionalFastConstructor(this);
			case MIXEDSPECTRUM: return new MixedSpectrumConstructor(this);
			case RANDOMTREE:return new RandomTreeDecoyTreeConstructor(this);
			default: return null;
		}
	}
	
	public String getDataName(){
		return getDataName(d);
	}
	
	public static String getDataName(DATA d){
		switch(d){
			case FILES: return "Files";
			case TREES: return "Trees";
			default: return null;
		}
	}
	
	public static String getDecoyTreeConstructorName(METHOD m){		
		switch(m){
//			case ORIGINALTREES: return DefaultDecoyTreeConstructor.methodNameLong;
			case ORIGINAL: return DefaultDecoyFileConstructor.methodNameLong;
			case REROOT: return RerootDecoyTreeConstructor.methodNameLong;
			case RANDOMPEAKS: return RandomPeaksConstructor.methodNameLong;
			case CONDITIONALPEAKS: return ConditionalPeaksConstructor.methodNameLong;
			case CONDITIONALFAST: return ConditionalFastConstructor.methodNameLong;
			case MIXEDSPECTRUM: return MixedSpectrumConstructor.methodNameLong;
			case RANDOMTREE: return RandomTreeDecoyTreeConstructor.methodNameLong;
			default: return null;
		}
	}
	
	public String getPostProcessNameLong(){
		String result="";
		if(pp.contains(POSTPROCESS.PPDiff))result+="PPDiff";
		return result;
	}
	
	public String getChargePath(){
		switch(c){
			case POS:return "pos";
			case NEG:return "neg";
			default: return null;
		}
	}
	
	public String getChargeShort(){
		switch(c){
			case POS:return "P";
			case NEG:return "N";
			default: return null;
		}
	}
	
	public String getChargeLong(){
		switch(c){
			case POS:return "Positive";
			case NEG:return "Negative";
			default: return null;
		}
	}
	
	public String getDatasetLong(){
		switch(ds){
			case AGILENT:return "Agilent";
			case METLIN:return "Metlin";
			case MASSBANKORBI:return "MassbankOrbi";
			case MASSBANKQTOF:return "MassbankQTof";
			case MASSBANK:return "Massbank";
			case GNPS:return "Gnps";
			default: return null;
		}
	}
	
	public String getDatasetShort(){
		switch(ds){
			case AGILENT:return "A";
			case METLIN:return "M";
			case MASSBANKORBI:return "O";
			case MASSBANKQTOF:return "Q";
			case MASSBANK:return "B";
			case GNPS:return "G";
			default: return null;
		}
	}
	
	public File getOutputFolder(){
		return new File(base+sep+baseDB+sep+getDatasetLong()+sep+getChargePath()+sep+getInstrumentName());
	}
	
	public File getOutputFolderRunningTime(){
		return new File(base+sep+baseRT+sep+getDatasetLong()+sep+getChargePath()+sep+getInstrumentName());
	}
	
	public File getOutputFolderDot(){
		return new File(getOutputFolder().getAbsolutePath()+sep+dot+addFolder);
	}
	
	public File getOutputFileRunningTime(){
		return new File(getOutputFolderRunningTime().getAbsolutePath()+sep+"runningTime_"+getInstrumentName()+addFolder+".txt");
	}
	
	public File getOutputFileDot(String fileName){
		return new File(getOutputFolderDot().getAbsolutePath()+sep+fileName+".dot");
	}
	
	public File getOutputFileMassbank(String fileName){
		return new File(getOutputFolderMassbank().getAbsolutePath()+sep+fileName+".txt");
	}
	
	public File getOutputFolderMassbank(){
		return new File(getOutputFolder().getAbsolutePath()+sep+massbank+addFolder);
	}
	
	public File getFolderOriginalFiles(){
		String baseOriginal=d.equals(DATA.TREES)?baseOriginalTrees:baseOriginalFiles;
		String folderOriginal=d.equals(DATA.TREES)?dot:ms;
		return new File(base+sep+baseOriginal+sep+getDatasetLong()+sep+getChargePath()+sep+folderOriginal);
	}
	
	public String getMethodNameLong(){
		return getDataName()+getDecoyTreeConstructor().getMethodNameLong();
	}
	
	public String getMassbankFirstID(){
		return getDatasetShort()+getChargeShort()+getMethodNameLong()+getPostProcessNameLong();
	}
	
	public String getMassbankID(int number){
		return getMassbankFirstID()+df.format(number);
	}
	
	public String getInstrumentName(){
		return getDatasetLong()+getChargeLong()+getMethodNameLong()+getPostProcessNameLong();
	}
	
	public Boolean isPositive(){
		switch(c){
			case POS:return true;
			case NEG:return false;
			default: return null;
		}
	}
	
	public Double getMinMass(){
		switch(ds){
			case AGILENT:return 50.0;
			case METLIN:return 29.0;
			case MASSBANKORBI:return 29.0;
			case MASSBANKQTOF:return 40.0;
			case MASSBANK: return 29.0;
			case GNPS: return 50.0;
			default: return null;
		}
	}
}
