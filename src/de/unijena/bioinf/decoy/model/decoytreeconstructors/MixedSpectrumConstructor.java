package de.unijena.bioinf.decoy.model.decoytreeconstructors;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import de.unijena.bioinf.babelms.dot.Graph;
import de.unijena.bioinf.babelms.dot.Vertex;
import de.unijena.bioinf.decoy.decoytrees.DecoySpectrum;
import de.unijena.bioinf.decoy.decoytrees.DefaultDecoySpectrum;
import de.unijena.bioinf.decoy.decoytrees.DefaultDecoyTree;
import de.unijena.bioinf.decoy.model.DecoyTreeVertex;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction;

public class MixedSpectrumConstructor extends DefaultDecoyTreeConstructor implements DecoySpectrumConstructor{
	
	static public final String methodNameLong="MixedSpectrum";
	Random r = new Random();

	public MixedSpectrumConstructor(ParametersDecoySpectrumConstruction p) {
		super(p);
	}
	
	public String getMethodNameLong(){
		return methodNameLong;
	}
	
	
//	public DecoySpectrum getDecoySpectrum(){		
//		List<double[]> allPeaksDecoy=new ArrayList<double[]>();		
//		DecoyTreeVertex root=inputGraph.getRoot();
//		allPeaksDecoy.add(new double[]{root.mz,root.intensity});
//				
//				
//		List<Graph> allGraphs=new ArrayList<Graph>();
//		allGraphs.addAll(p.getAllGraphsInOriginal().values());
//		allGraphs.remove(inputGraph);
//		int numberGraphs=0;
//		List<double[]> allPeaksDecoyTmp=new ArrayList<double[]>();	
//		while(numberGraphs<2||allPeaksDecoyTmp.size()<3*inputGraph.getVertices().size()){
//			int i=r.nextInt(allGraphs.size());
//			Graph currGraph=allGraphs.get(i);
//			DefaultDecoyTree dt=new DefaultDecoyTree(currGraph);
//			allGraphs.remove(i);
//			boolean GraphUsed=false;
//			for(DecoyTreeVertex v:dt.getVertices()){
//				double mass=v.mz;
//				if(mass<root.mz){
//					allPeaksDecoyTmp.add(new double[]{mass, v.intensity});
//					GraphUsed=true;
//				}
//			}			
//			if(GraphUsed)numberGraphs++;
//		}
//		
//		for (int j=0;j<inputGraph.getVertices().size()-1;j++){
//			int i=r.nextInt(allPeaksDecoyTmp.size());
//			allPeaksDecoy.add(allPeaksDecoyTmp.get(i));
//			allPeaksDecoyTmp.remove(i);
//		}	
//		
//		DefaultDecoySpectrum s=new DefaultDecoySpectrum(inputGraph.getParent(), allPeaksDecoy);
//		return s;
//	}
	
	public DecoySpectrum getDecoySpectrum(){		
		List<double[]> allPeaksDecoy=new ArrayList<double[]>();		
		DecoyTreeVertex root=inputGraph.getRoot();
		allPeaksDecoy.add(new double[]{root.mz,root.intensity});
				
				
		List<List<double[]>> allPeaks=new ArrayList<List<double[]>>();
		allPeaks.addAll(p.getAllPeaksInOriginal().values());
		int numberGraphs=0;
		List<double[]> allPeaksDecoyTmp=new ArrayList<double[]>();	
		while(numberGraphs<2||allPeaksDecoyTmp.size()<3*inputGraph.getVertices().size()){
			int i=r.nextInt(allPeaks.size());
			List<double[]> currPeaks=allPeaks.get(i);
			allPeaks.remove(i);
			boolean GraphUsed=false;
			for(double[] d:currPeaks){
				double mass=d[0];
				if(mass<root.mz){
					allPeaksDecoyTmp.add(new double[]{mass, d[1]});
					GraphUsed=true;
				}
			}			
			if(GraphUsed)numberGraphs++;
		}
		
		for (int j=0;j<inputGraph.getVertices().size()-1;j++){
			int i=r.nextInt(allPeaksDecoyTmp.size());
			allPeaksDecoy.add(allPeaksDecoyTmp.get(i));
			allPeaksDecoyTmp.remove(i);
		}	
		
		DefaultDecoySpectrum s=new DefaultDecoySpectrum(inputGraph.getParent(), allPeaksDecoy);
		return s;
	}


}
