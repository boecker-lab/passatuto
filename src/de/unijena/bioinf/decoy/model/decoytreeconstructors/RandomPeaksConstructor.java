package de.unijena.bioinf.decoy.model.decoytreeconstructors;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import de.unijena.bioinf.babelms.dot.Graph;
import de.unijena.bioinf.decoy.decoytrees.DecoySpectrum;
import de.unijena.bioinf.decoy.decoytrees.DefaultDecoySpectrum;
import de.unijena.bioinf.decoy.decoytrees.DefaultDecoyTree;
import de.unijena.bioinf.decoy.model.DecoyTreeVertex;
import de.unijena.bioinf.deocy.Utils;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction;

public class RandomPeaksConstructor extends DefaultDecoyTreeConstructor implements DecoySpectrumConstructor{
	
	static public final String methodNameLong="RandomPeaks";
	Random r = new Random();

	public RandomPeaksConstructor(ParametersDecoySpectrumConstruction p) {
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
		
	
	public DecoySpectrum getDecoySpectrum(){
		List<double[]> allPeaksDecoy=new ArrayList<double[]>();
		DecoyTreeVertex root=inputGraph.getRoot();
		allPeaksDecoy.add(new double[]{root.mz,root.intensity});
		
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
	
				
//		List<double[]> allPeaksDB=getAllPeaksFromDB(inputGraph.getRoot().mz);
//		for(int i=0;i<inputGraph.getVertices().size()-1;i++){
//			int index=r.nextInt(allPeaksDB.size());			
//			allPeaksDecoy.add(allPeaksDB.get(index));
//			double addedMass=allPeaksDB.get(index)[0];
//			Iterator<double[]> it=allPeaksDB.iterator();
//			while(it.hasNext()){
//				if(Double.compare(it.next()[0],addedMass)==0)it.remove();
//			}
//		}
		
		DefaultDecoySpectrum s=new DefaultDecoySpectrum(inputGraph.getParent(), allPeaksDecoy);
		return s;
	}
	
	List<double[]> getAllPeaksFromDB(double maxMass){
		List<double[]> peaks=new ArrayList<double[]>();
		for(Graph g:p.getAllGraphsInOriginal().values()){
			DefaultDecoyTree dt=new DefaultDecoyTree(g);
			peaks.addAll(dt.getPeaks(maxMass));
		}
		return peaks;
	}

}
