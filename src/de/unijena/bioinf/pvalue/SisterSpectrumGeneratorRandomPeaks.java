package de.unijena.bioinf.pvalue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
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

public class SisterSpectrumGeneratorRandomPeaks{
	
	Random r=new Random();
	int numberPeaks;
	Map<Double,List<double[]>> conditionalPeaksFromDB;
	List<double[]> allPeaksFromDB;
	
	public SisterSpectrumGeneratorRandomPeaks(int numberPeaks, ParametersDecoySpectrumConstruction p){
		this.numberPeaks=numberPeaks;
		this.conditionalPeaksFromDB=ConditionalPeaksConstructor.getAllConditionalPeaksFromDB(Double.MAX_VALUE,p.getAllGraphsInOriginal().values());
		this.allPeaksFromDB=ConditionalPeaksConstructor.getAllPeaksFromDB(Double.MAX_VALUE,p.getAllGraphsInOriginal().values());
	}
	
	public List<double[]> getStartSpectrum(){
		
		List<double[]> allPeaksFromDBTmp=new ArrayList<double[]>(allPeaksFromDB);
		
		List<double[]> allPeaksDecoy=new ArrayList<double[]>();
		for(int i=0;i<numberPeaks;i++){
			int rand=r.nextInt(allPeaksFromDBTmp.size());
			double[] nextPeak=allPeaksFromDBTmp.get(rand);
			allPeaksDecoy.add(nextPeak);
			Iterator<double[]> it=allPeaksFromDBTmp.iterator();
			while(it.hasNext())if(Double.compare(it.next()[0],nextPeak[0])==0)it.remove();
		}
		
		Collections.sort(allPeaksDecoy,new SpectrumComparator());
		return allPeaksDecoy;
	}
	
	public List<double[]> getSisterSpectrum(List<double[]> originalSpectrum){
		List<double[]> spectrum=new ArrayList<double[]>(originalSpectrum);
		int rand=r.nextInt(spectrum.size());
		spectrum.remove(rand);
		
		List<double[]> allPeaksFromDBTmp=new ArrayList<double[]>(allPeaksFromDB);		
		
		boolean added=false;
		while(!added){
			rand=r.nextInt(allPeaksFromDBTmp.size());
			double[] nextPeak=allPeaksFromDBTmp.get(rand);
			boolean alreadyExists=false;
			for(double[] d:spectrum) if(Double.compare(d[0],nextPeak[0])==0)alreadyExists=true;
			if(!alreadyExists){
				spectrum.add(nextPeak);
				added=true;
			}	
		}				
		
		Collections.sort(spectrum, new SpectrumComparator());
		return spectrum;
	}
	
}