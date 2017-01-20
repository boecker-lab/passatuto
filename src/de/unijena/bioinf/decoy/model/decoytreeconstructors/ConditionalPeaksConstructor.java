package de.unijena.bioinf.decoy.model.decoytreeconstructors;

import java.util.ArrayList;
import java.util.Collection;
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

public class ConditionalPeaksConstructor extends DefaultDecoyTreeConstructor implements DecoySpectrumConstructor{
	static private TreeMap<Double, List<Graph>> m=null;
	static public final String methodNameLong="ConditionalPeaks";
	Random r = new Random();

	public ConditionalPeaksConstructor(ParametersDecoySpectrumConstruction p) {
		super(p);
	}
	
	public String getMethodNameLong(){
		return methodNameLong;
	}
	
	public DecoySpectrum getDecoySpectrum(){		
		List<double[]> allPeaksDecoy=new ArrayList<double[]>();
		DecoyTreeVertex root=inputGraph.getRoot();
		allPeaksDecoy.add(new double[]{root.mz,root.intensity});
		
		Map<Double,List<double[]>> allPeaksDB=getAllConditionalPeaksFromDB(inputGraph.getRoot().mz);

		for(int i=0;i<inputGraph.getVertices().size()-1;i++){
			
			List<double[]> currentPeaks=new ArrayList<double[]>();
//			for(double[] p:allPeaksDecoy){
			double[] p=allPeaksDecoy.get(allPeaksDecoy.size()-1);
			double error=Utils.getAbsoluteErrorForMass(p[0], this.p.ppm, this.p.ae);
			for(Entry<Double, List<double[]>> peaksDB:allPeaksDB.entrySet()){
				if(Math.abs(p[0]-peaksDB.getKey())<error)currentPeaks.addAll(peaksDB.getValue());
			}
			
			int index=r.nextInt(currentPeaks.size());			
			allPeaksDecoy.add(currentPeaks.get(index));
			double addedMass=currentPeaks.get(index)[0];
			error=Utils.getAbsoluteErrorForMass(addedMass, this.p.ppm, this.p.ae);
			for(List<double[]> peaks:allPeaksDB.values()){
				Iterator<double[]> it=peaks.iterator();
				while(it.hasNext()){
					if(Math.abs(it.next()[0]-addedMass)<error)it.remove();
				}
			}
		}
		DefaultDecoySpectrum s=new DefaultDecoySpectrum(inputGraph.getParent(), allPeaksDecoy);
		return s;
	}
	
	Map<Double,List<double[]>> getAllConditionalPeaksFromDB(double maxMass){
		return getAllConditionalPeaksFromDB(maxMass, p.getAllGraphsInOriginal().values());
	}
	
	public static Map<Double,List<double[]>> getAllConditionalPeaksFromDB(double maxMass, Collection<Graph> graphsInOriginalDB){
		Map<Double, List<double[]>> result=new HashMap<Double, List<double[]>>();
		for(Graph g:graphsInOriginalDB){
			DefaultDecoyTree dt=new DefaultDecoyTree(g);
			List<double[]> peaks=dt.getPeaksIncludingMax(maxMass);
			for(int i=0;i<peaks.size();i++){
				for(int j=0;j<peaks.size();j++){					
					if(i!=j){
						if(!result.containsKey(peaks.get(i)[0]))result.put(peaks.get(i)[0],new ArrayList<double[]>());
						result.get(peaks.get(i)[0]).add(peaks.get(j));
					}
				}
			}
		}
		return result;
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
