package de.unijena.bioinf.decoy.model;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.deocy.Utils;

public class MassBank implements java.io.Serializable, Comparable<MassBank>{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -4955544804945648760L;

	/**
	 * 
	 */

	static String sep=System.getProperty("file.separator");
		
	public String massbankID;
	public String recordTitle;
	public String instrument;
	public MolecularFormula mf=null;
	public Double massParentIon=null;
	public String inchi=null;
	public String inchiBackup=null;
	public boolean isPositive;
	public List<double[]> peaks;
	Double idealParentIon=null;
	
	public MassBank(File f) throws Exception{
		this(f, true);
	}
	
	public MassBank(File f, boolean withPeaks) throws Exception{
		BufferedReader br=new BufferedReader(new FileReader(f));
		String line;
		while((line=br.readLine())!=null){
			if(line.startsWith("ACCESSION"))this.massbankID=line.substring(11);
			if(line.startsWith("RECORD_TITLE"))this.recordTitle=line.substring(14);
			if(line.startsWith("CH$FORMULA"))if(this.mf==null)this.mf=MolecularFormula.parse(line.substring(12));
			if(line.startsWith("AC$INSTRUMENT"))this.instrument=line.substring(15);
			if(line.startsWith("MS$FOCUSED_ION: PRECURSOR_TYPE")){
				if(line.substring(31).equals("[M+H]+"))isPositive=true;
				else isPositive=false;
				if(this.mf!=null){
					double hAdduct=isPositive?1.00728:-1.00728;		
					idealParentIon=this.mf.getMass()+hAdduct;
				}
			}
				
			if(line.startsWith("CH$IUPAC")){
				String l=line.substring(10);
				this.inchiBackup=l;
				if(!l.equals("N/A")){
					String tmp[]=l.split("/");
					String first=tmp[1];
					int i=first.indexOf('.');
					if(i>=0)first=first.substring(0,i);
//					this.mf=MolecularFormula.parse(first);
					String second=tmp[2];
					i=second.indexOf(';');
					if(i>=0)second=second.substring(0,i);
					this.inchi=first+"/"+second;
				}
			}			
			if(line.startsWith("PK$PEAK")){
				peaks=new ArrayList<double[]>();
				while(!(line=br.readLine()).equals("//")){
					String[] p=line.trim().split(" ");
					double mz=Double.parseDouble(p[0]);
					if(massParentIon==null||Math.abs(massParentIon-idealParentIon)>Math.abs(mz-idealParentIon))massParentIon=mz;
					if(withPeaks)peaks.add(new double[]{mz,Double.parseDouble(p[1]),Double.parseDouble(p[2])});
				}
			}
		}
		br.close();
	}
	
	public void removeParentPeaks(double ppm, double ae){
		Iterator<double[]> peakIt=peaks.iterator();
		while(peakIt.hasNext()){
			double[] p=peakIt.next();
			double error=Utils.getAbsoluteErrorForMass(p[0], ppm, ae);
			if(Math.abs(p[0]-idealParentIon)<error)peakIt.remove();
		}
		
	}
	
	public MassBank(List<double[]> p, MolecularFormula mf, double massParentIon, String massbankID, String recordTitle, String instrument, String inchi, boolean isPositive, boolean useParentPeakDiffs) throws IOException {		
		this.massbankID=massbankID;
		this.recordTitle=recordTitle;
		this.mf=mf;
		this.massParentIon=massParentIon;
		this.instrument=instrument;
		this.inchi=inchi;
		this.isPositive=isPositive;
		
		peaks=new ArrayList<double[]>(p);
		
		if(useParentPeakDiffs){
			List<double[]> newPeaks=new ArrayList<double[]>();
			double maxIntensity=0;
			for(int i=0;i<p.size();i++){
				maxIntensity=Math.max(p.get(i)[1], maxIntensity);
			}
			for (double[] v : p) {
				double diffMass=mf.getMass()-v[0];
				if(diffMass>0){
					newPeaks.add(new double[]{mf.getMass()-v[0], maxIntensity-v[1]});
				}
			}
			peaks.addAll(newPeaks);
		}

	}
	
	public void writeSpectrumAsMassbank(File folder) throws IOException {
		if(!folder.exists())folder.mkdirs();
		BufferedWriter writer=new BufferedWriter(new FileWriter(folder+sep+massbankID+".txt"));
		double hAdduct=isPositive?1.00728:-1.00728;
		final BufferedWriter bw = new BufferedWriter(writer);
		bw.write("ACCESSION: "+massbankID+"\n");
		bw.write("RECORD_TITLE: "+recordTitle+"\n");
		SimpleDateFormat sdf = new SimpleDateFormat("yyyy.MM.dd");
		bw.write("DATE: "+sdf.format(Calendar.getInstance().getTime())+"\n");
		bw.write("AUTHORS: Administrator\n");
		bw.write("LICENSE: CC BY-SA\n");
		bw.write("CH$NAME: "+recordTitle+"\n");
		bw.write("CH$COMPOUND_CLASS: Natural Product\n");
		bw.write("CH$FORMULA: "+this.mf+"\n");
		bw.write("CH$EXACT_MASS: "+this.mf.getMass()+"\n");
		bw.write("CH$SMILES: N/A\n");
		bw.write("CH$IUPAC: "+(inchi==null?"N/A":inchi)+"\n");
		bw.write("CH$LINK: N/A\n");
		bw.write("AC$INSTRUMENT: "+instrument+"\n");
		bw.write("AC$INSTRUMENT_TYPE: "+instrument+"\n");
		bw.write("AC$MASS_SPECTROMETRY: MS_TYPE MS2\n");
		bw.write("AC$MASS_SPECTROMETRY: ION_MODE POSITIVE\n");
		bw.write("AC$MASS_SPECTROMETRY: COLLISION_ENERGY 10 V\n");
		bw.write("MS$FOCUSED_ION: PRECURSOR_M/Z "+Math.round(this.mf.getMass()+hAdduct)+"\n");
		bw.write("MS$FOCUSED_ION: PRECURSOR_TYPE [M+H]+\n");
		bw.write("PK$NUM_PEAK: "+peaks.size()+"\n");
		bw.write("PK$PEAK: m/z int. rel.int.\n");
		
		
		double maxIntensity=0;
		for(int i=0;i<peaks.size();i++){
			maxIntensity=Math.max(peaks.get(i)[1], maxIntensity);
		}
		
		Map<Double,Double> intensityMap=new TreeMap<Double,Double>();
		for(int i=0;i<peaks.size();i++){
			double mass=Math.round(peaks.get(i)[0]*1000000)/1000000.0;
			double intensity=peaks.get(i)[1];
			if(!intensityMap.containsKey(mass)){
				intensityMap.put(mass,intensity);
			}else{
				if(intensityMap.get(mass)<intensity){
					intensityMap.put(mass,intensity);
				}				
			}
		}

		List<double[]> p =new ArrayList<double[]>();
		for (Entry<Double,Double> e : intensityMap.entrySet()) {
			p.add(new double[]{e.getKey()+hAdduct, e.getValue(), e.getValue()/maxIntensity*999});			
		}

		for (double[] v : p) {
			bw.write("  "+v[0]+" "+v[1]+" "+Math.round(v[2])+"\n");
		}
		
		bw.write("//\n");
		bw.flush();
		bw.close();
	}
	
	public boolean isDecoy(){
		return (!instrument.contains("Original"));
	}
	
	public boolean hasEqualInChiKey(MassBank mb){
		if(this.inchi==null||mb.inchi==null)return false;
		return (this.mf.equals(mb.mf)&&this.inchi.equals(mb.inchi));
	}
	
	public boolean hasEqualMass(MassBank mb, double ppm, double ae){
		double error=Utils.getAbsoluteErrorForMass(getMassParentIon(), ppm, ae);
//		if(Math.abs(getMassParentIon()-mb.getMassParentIon())>error){
//		if(Math.abs(getMassParentIon()-mb.idealParentIon)>error){
		if(Math.abs(idealParentIon-mb.idealParentIon)>error){
			return false;
		}
		return true;
	}
	
	public double getMassParentIon(){
		return massParentIon;
	}

	@Override
	public int compareTo(MassBank o) {
		return this.massbankID.compareTo(o.massbankID);
	}
}
