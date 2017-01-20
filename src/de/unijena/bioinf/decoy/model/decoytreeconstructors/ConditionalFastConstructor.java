package de.unijena.bioinf.decoy.model.decoytreeconstructors;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import de.unijena.bioinf.babelms.dot.Graph;
import de.unijena.bioinf.decoy.decoytrees.DecoySpectrum;
import de.unijena.bioinf.decoy.decoytrees.DefaultDecoySpectrum;
import de.unijena.bioinf.decoy.decoytrees.DefaultDecoyTree;
import de.unijena.bioinf.decoy.model.DecoyTreeVertex;
import de.unijena.bioinf.deocy.Utils;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction;

public class ConditionalFastConstructor extends DefaultDecoyTreeConstructor implements DecoySpectrumConstructor{
//	static private TreeMap<Double, List<Graph>> m=null;
	static private TreeMap<Double, List<List<double[]>>> m=null;
	static public final String methodNameLong="ConditionalFast";
	Random r = new Random();

	public ConditionalFastConstructor(ParametersDecoySpectrumConstruction p) {
		super(p);
	}

	public String getMethodNameLong(){
		return methodNameLong;
	}

//	public DecoySpectrum getDecoySpectrum(){		
//		if(m==null){
//			m=new TreeMap<Double, List<Graph>>();
//			for(Graph g:p.getAllGraphsInOriginal().values()){
//				List<double[]> peaks=DefaultDecoyTree.getPeaksIncludingMax(g,Double.POSITIVE_INFINITY);					
//				for(double[] p:peaks){
//					if(!m.containsKey(p[0]))m.put(p[0],new ArrayList<Graph>());
//					m.get(p[0]).add(g);
//				}
//			}
//		}
//
//		List<double[]> allPeaksDecoy=new ArrayList<double[]>();
//		DecoyTreeVertex root=inputGraph.getRoot();
//
//		List<double[]> currentPeaks=new ArrayList<double[]>();
//		currentPeaks.add(new double[]{root.mz,root.intensity});
//		int size=inputGraph.getVertices().size();
//		for(int i=0;i<size;i++){
//			System.out.println(i+" of "+size);
//			int index=r.nextInt(currentPeaks.size());
//			allPeaksDecoy.add(currentPeaks.get(index));
//			double addedMass=currentPeaks.get(index)[0];
//			currentPeaks.remove(index);
//
//			List<double[]> peaksTmp=new ArrayList<double[]>();
//			double error=Utils.getAbsoluteErrorForMass(addedMass, this.p.ppm, this.p.ae);
//
//			SortedMap<Double, List<Graph>> mSub=m.subMap(m.ceilingKey(addedMass-error), true, m.floorKey(addedMass+error), true);
//			List<Graph> graphsSub=new ArrayList<Graph>();
//			for(List<Graph> x:mSub.values())graphsSub.addAll(x);
//			boolean found=false;
//			List<double[]> peaks=new ArrayList<double[]>();
//			while(!found&&graphsSub.size()>0){
//				index=r.nextInt(graphsSub.size());
//				Graph g=graphsSub.get(index);
//				peaks=DefaultDecoyTree.getPeaksIncludingMax(g,root.mz);
//				if(peaks.size()>0)found=true;
//				graphsSub.remove(index);
//			}			
//			peaksTmp.addAll(peaks);
//			for(int t=0;t<Math.min(5,peaksTmp.size());t++){
//				index=r.nextInt(peaksTmp.size());
//				currentPeaks.add(peaksTmp.get(index));
//				peaksTmp.remove(index);
//			}
//		}
//		
//		Set<Double> masses=new HashSet<Double>();
//		Iterator<double[]> it=allPeaksDecoy.iterator();
//		while(it.hasNext()){
//			double mass=Math.round(it.next()[0]*1000000)/1000000.0;
//			boolean added=masses.add(mass);
//			if(!added)it.remove();		
//		}
//		
//		List<Graph> allGraphs=new ArrayList<Graph>();
//		allGraphs.addAll(p.getAllGraphsInOriginal().values());
//		while(allPeaksDecoy.size()<inputGraph.getVertices().size()){
//			int i=r.nextInt(allGraphs.size());
//			Graph currGraph=allGraphs.get(i);
//			DefaultDecoyTree dt=new DefaultDecoyTree(currGraph);
//			int i2=r.nextInt(dt.getVertices().size());			
//			double mass=dt.getVertices().get(i2).mz;
//			if(mass<root.mz){
//				boolean alreadyAdded=false;
//				double error=Utils.getAbsoluteErrorForMass(mass, p.ppm, p.ae);
//				for(double d[]:allPeaksDecoy){
//					if(Math.abs(mass-d[0])<error)alreadyAdded=true;
//				}
//				if(!alreadyAdded)allPeaksDecoy.add(new double[]{mass, dt.getVertices().get(i2).intensity});
//			}
//		}
//		
//		DefaultDecoySpectrum s=new DefaultDecoySpectrum(inputGraph.getParent(), allPeaksDecoy);
//		return s;
//	}
	
	public static void setM(ParametersDecoySpectrumConstruction p){
		if(m==null){
			m=new TreeMap<Double, List<List<double[]>>>();
			for(List<double[]> peaks:p.getAllPeaksInOriginal().values()){					
				for(double[] peak:peaks){
					if(!m.containsKey(peak[0]))m.put(peak[0],new ArrayList<List<double[]>>());
					m.get(peak[0]).add(peaks);
				}
			}
		}
	}
	
	
	public DecoySpectrum getDecoySpectrum(){		
		if(m==null){
			setM(p);
		}

		List<double[]> allPeaksDecoy=new ArrayList<double[]>();
		DecoyTreeVertex root=inputGraph.getRoot();

		List<double[]> currentPeaks=new ArrayList<double[]>();
		currentPeaks.add(new double[]{root.mz,root.intensity});
		int size=inputGraph.getVertices().size();
		for(int i=0;i<size;i++){
			System.out.println(i+" of "+size);
			int index=r.nextInt(currentPeaks.size());
			allPeaksDecoy.add(currentPeaks.get(index));
			double addedMass=currentPeaks.get(index)[0];
			currentPeaks.remove(index);

			List<double[]> peaksTmp=new ArrayList<double[]>();
			double error=Utils.getAbsoluteErrorForMass(addedMass, this.p.ppm, this.p.ae);

			SortedMap<Double, List<List<double[]>>> mSub=m.subMap(m.ceilingKey(addedMass-error), true, m.floorKey(addedMass+error), true);
			List<List<double[]>> peaksSub=new ArrayList<List<double[]>>();
			for(List<List<double[]>> x:mSub.values())peaksSub.addAll(x);
			boolean found=false;
			List<double[]> peaks=new ArrayList<double[]>();
			while(!found&&peaksSub.size()>0){
				index=r.nextInt(peaksSub.size());
				List<double[]> g=peaksSub.get(index);
				for(double[] d:g)if(d[0]<=root.mz)peaks.add(d);
				if(peaks.size()>0)found=true;
				peaksSub.remove(index);
			}			
			peaksTmp.addAll(peaks);
			for(int t=0;t<Math.min(5,peaksTmp.size());t++){
				index=r.nextInt(peaksTmp.size());
				currentPeaks.add(peaksTmp.get(index));
				peaksTmp.remove(index);
			}
		}
		
		Set<Double> masses=new HashSet<Double>();
		Iterator<double[]> it=allPeaksDecoy.iterator();
		while(it.hasNext()){
			double mass=Math.round(it.next()[0]*1000000)/1000000.0;
			boolean added=masses.add(mass);
			if(!added)it.remove();		
		}
		
		List<List<double[]>> allPeaks=new ArrayList<List<double[]>>();
		allPeaks.addAll(p.getAllPeaksInOriginal().values());
		while(allPeaksDecoy.size()<inputGraph.getVertices().size()){
			int i=r.nextInt(allPeaks.size());
			List<double[]> currPeaks=allPeaks.get(i);
			int i2=r.nextInt(currPeaks.size());			
			double mass=currPeaks.get(i2)[0];
			if(mass<root.mz){
				boolean alreadyAdded=false;
				double error=Utils.getAbsoluteErrorForMass(mass, p.ppm, p.ae);
				for(double d[]:allPeaksDecoy){
					if(Math.abs(mass-d[0])<error)alreadyAdded=true;
				}
				if(!alreadyAdded)allPeaksDecoy.add(new double[]{mass, currPeaks.get(i2)[1]});
			}
		}
		
		DefaultDecoySpectrum s=new DefaultDecoySpectrum(inputGraph.getParent(), allPeaksDecoy);
		return s;
	}


	public static List<double[]> getAllPeaksFromDB(double maxMass, Collection<Graph> graphsInOriginalDB){
		List<double[]> result=new ArrayList<double[]>();
		for(Graph g:graphsInOriginalDB){
			DefaultDecoyTree dt=new DefaultDecoyTree(g);
			List<double[]> peaks=dt.getPeaksIncludingMax(maxMass);
			result.addAll(peaks);		
		}
		return result;
	}

}
