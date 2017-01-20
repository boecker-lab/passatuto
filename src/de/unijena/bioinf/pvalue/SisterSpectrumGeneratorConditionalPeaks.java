package de.unijena.bioinf.pvalue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.Map.Entry;
import java.util.TreeMap;

import de.unijena.bioinf.babelms.dot.Graph;
import de.unijena.bioinf.babelms.dot.Vertex;
import de.unijena.bioinf.decoy.decoytrees.DefaultDecoySpectrum;
import de.unijena.bioinf.decoy.decoytrees.DefaultDecoyTree;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction;
import de.unijena.bioinf.decoy.model.DecoyTreeVertex;
import de.unijena.bioinf.decoy.model.decoytreeconstructors.ConditionalPeaksConstructor;
import de.unijena.bioinf.deocy.Utils;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction.CHARGE;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction.DATASET;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction.METHOD;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction.POSTPROCESS;

public class SisterSpectrumGeneratorConditionalPeaks{
	
	Random r=new Random();
	double[] parentPeak;
	int numberPeaks;
	Map<Double,List<double[]>> conditionalPeaksFromDB;
	List<double[]> allPeaksFromDB;
	
	public SisterSpectrumGeneratorConditionalPeaks(double[] parentPeak, int numberPeaks, ParametersDecoySpectrumConstruction p){
		this.parentPeak=parentPeak;
		this.numberPeaks=numberPeaks;
		this.conditionalPeaksFromDB=ConditionalPeaksConstructor.getAllConditionalPeaksFromDB(parentPeak[0],p.getAllGraphsInOriginal().values());
		this.allPeaksFromDB=ConditionalPeaksConstructor.getAllPeaksFromDB(parentPeak[0],p.getAllGraphsInOriginal().values());
	}
	
//	public static void main(String[] args){
//		
//		Locale.setDefault(Locale.US);
//		
//		Random r=new Random();
//		int numberPeaks=10;
//		
//		ParametersDecoySpectrumConstruction p=new ParametersDecoySpectrumConstruction(CHARGE.POS, DATASET.AGILENT, METHOD.ORIGINAL, Arrays.asList(new POSTPROCESS[]{}));		
//		
//		List<double[]> allPeaks=new ArrayList<double[]>();
//		Collection<Graph> graphs=p.getAllGraphsInOriginalDB().values();
//		for(Graph g:graphs){
//			for(Vertex v:g.getVertices()){
//				DecoyTreeVertex dtv=new DecoyTreeVertex(v);
//				allPeaks.add(new double[]{dtv.mf.getMass(),dtv.intensity});
//			}
//		}
//		double[] parentPeak=allPeaks.get(r.nextInt(allPeaks.size()));
//	
//		
//		SisterSpectrumGenerator g=new SisterSpectrumGenerator(parentPeak, numberPeaks, p);
//		List<double[]> startSpectrum=g.getStartSpectrum();
//		List<double[]> sisterSpectrum=g.getSisterSpectrum(startSpectrum);
//	}
	
	public List<double[]> getStartSpectrum(){
		
//		Map<Double,List<double[]>> conditionalPeaksFromDBTmp=new HashMap<Double, List<double[]>>();
//		for(Entry<Double, List<double[]>> e:conditionalPeaksFromDB.entrySet()){
//			List<double[]> newList=new ArrayList<double[]>(e.getValue());
//			conditionalPeaksFromDBTmp.put(e.getKey(), newList);
//		}		
//		
//		List<double[]> allPeaksDecoy=new ArrayList<double[]>();
//		allPeaksDecoy.add(parentPeak);
//		
//		for(int i=0;i<numberPeaks-1;i++){
//			
//			List<double[]> currentPeaks=new ArrayList<double[]>();
//			for(double[] p:allPeaksDecoy){
//				if(conditionalPeaksFromDBTmp.get(p[0])!=null){
//					for(double[] p2:conditionalPeaksFromDBTmp.get(p[0])){
//						if(Double.compare(parentPeak[0],p2[0])!=0)currentPeaks.add(p2);
//					}
//				}
//			}
//			
//			if(currentPeaks.isEmpty()){
//				int index=r.nextInt(conditionalPeaksFromDBTmp.values().size());
//				List<double[]> tmp=(List<double[]>)conditionalPeaksFromDBTmp.values().toArray()[index];
//				index=r.nextInt(tmp.size());			
//				currentPeaks.add(tmp.get(index));
//			}
//			
//			int index=r.nextInt(currentPeaks.size());			
//			allPeaksDecoy.add(currentPeaks.get(index));
//			double addedMass=currentPeaks.get(index)[0];
//			for(List<double[]> peaks:conditionalPeaksFromDBTmp.values()){
//				Iterator<double[]> it=peaks.iterator();
//				while(it.hasNext()){
//					if(Double.compare(it.next()[0],addedMass)==0)it.remove();
//				}
//			}
//		}
//		Collections.sort(allPeaksDecoy,new SpectrumComparator());
//		return allPeaksDecoy;
		
		List<double[]> allPeaksFromDBTmp=new ArrayList<double[]>(allPeaksFromDB);
		
		List<double[]> allPeaksDecoy=new ArrayList<double[]>();
		allPeaksDecoy.add(parentPeak);
		Iterator<double[]> it=allPeaksFromDBTmp.iterator();
		while(it.hasNext())if(Double.compare(it.next()[0],parentPeak[0])==0)it.remove();
		for(int i=0;i<numberPeaks-1;i++){
			int rand=r.nextInt(allPeaksFromDBTmp.size());
			double[] nextPeak=allPeaksFromDBTmp.get(rand);
			allPeaksDecoy.add(nextPeak);
			it=allPeaksFromDBTmp.iterator();
			while(it.hasNext())if(Double.compare(it.next()[0],nextPeak[0])==0)it.remove();
		}
		
		Collections.sort(allPeaksDecoy,new SpectrumComparator());
		return allPeaksDecoy;
	}
	
	public List<double[]> getSisterSpectrum(List<double[]> originalSpectrum){
		List<double[]> spectrum=new ArrayList<double[]>(originalSpectrum);
		Collections.sort(spectrum, new SpectrumComparator());
		Map<double[],Integer> deletingOnePeak=new HashMap<double[],Integer>();
		int allN=0;		
		for(int i=0;i<spectrum.size()-1;i++){
			double[] deletedPeak=spectrum.get(i);
			int n=0;
			for(int j=0;j<spectrum.size();j++){				
				if(i!=j){
					double[] currentPeak=spectrum.get(j);
					List<double[]> condPeaks=conditionalPeaksFromDB.get(currentPeak[0]);
					if(condPeaks!=null){
						for(double[] p:condPeaks)
							if(Double.compare(p[0],deletedPeak[0])==0)
								n++;
					}
				}
			}
			deletingOnePeak.put(deletedPeak,n);
			allN+=n;
		}
		
		int k=r.nextInt(allN);
		int i=0;
		double[] deletedPeak=null;
		for(Entry<double[], Integer> e:deletingOnePeak.entrySet()){	
			if(i+e.getValue()>k){
				deletedPeak=e.getKey();
				break;
			}else{
				i+=e.getValue();
			}
		}
		
		Iterator<double[]> iteratorSpectrum=spectrum.iterator();
		while(iteratorSpectrum.hasNext()){
			if(iteratorSpectrum.next()==deletedPeak)iteratorSpectrum.remove();
		}	
		
		List<double[]> currentPeaks=new ArrayList<double[]>();
		for(double[] p:spectrum){
			if(conditionalPeaksFromDB.get(p[0])!=null)
				currentPeaks.addAll(conditionalPeaksFromDB.get(p[0]));
		}
		currentPeaks.addAll(allPeaksFromDB);
		
		boolean added=false;
		while(!added){
			int index=r.nextInt(currentPeaks.size());
			double[] peak=currentPeaks.get(index);
			boolean alreadyExists=false;
			for(double[] s:spectrum)if(Double.compare(s[0],peak[0])==0)alreadyExists=true;
			if(!alreadyExists){
				spectrum.add(peak);
				added=true;
			}
		}
		
		Collections.sort(spectrum, new SpectrumComparator());
		return spectrum;
	}
	
}