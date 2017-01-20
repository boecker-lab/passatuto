package de.unijena.bioinf.deocy;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;

import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.babelms.dot.Edge;
import de.unijena.bioinf.babelms.dot.Graph;
import de.unijena.bioinf.babelms.dot.Parser;
import de.unijena.bioinf.babelms.dot.Vertex;
import de.unijena.bioinf.decoy.decoytrees.DefaultDecoyTree;
import de.unijena.bioinf.decoy.model.DecoyTreeEdge;
import de.unijena.bioinf.decoy.model.DecoyTreeVertex;

public class Utils {

	public static Map<File,Graph> readGraphs(File input){
		return readGraphs(input,true);
	}

	public static Map<File,Graph> readGraphs(File input,boolean storeGraphs){
		Map<File,Graph> result=new HashMap<File,Graph>();

		for(File f:input.listFiles()){
			if(f.getName().endsWith(".dot")){
				try{
					Graph g=Parser.parse(new FileReader(f));
					if(g!=null&&g.getVertices().size()>1){
						result.put(f,g);
						System.out.println("Graph for "+f.getName()+" is read.");
					}else{
						System.err.println("Graph for "+f.getName()+" is empty.");
					}
				}catch(IOException e){
					System.err.println("Graph "+f+" was not read!");
				}				
			}else if(f.getName().endsWith(".ms")||f.getName().endsWith(".txt")){
				try{					
					BufferedReader br=new BufferedReader(new FileReader(f));
					
					Graph g=new Graph();				
					List<Peak> peaks=new ArrayList<Peak>();
					String line="";
					MolecularFormula mf=null;
					double ionization=0;
					boolean intrinsicCharged=false;
					Double parentMass=null;
					double maxIntensity=0;
					while((line=br.readLine())!=null){
						if(line.startsWith(">formula")){
							mf=MolecularFormula.parse(line.substring(9));
							if(mf.maybeCharged()){
								ionization=-0.0019666;
								intrinsicCharged=true;
							}
						}
						if(line.startsWith(">parentmass")){
							parentMass=Double.parseDouble(line.substring(12));
						}
						if(line.startsWith(">ionization")&&!intrinsicCharged){
							String ion=line.substring(12);
							if(ion.equals("[M+H]+"))ionization=1.00728;
							else if(ion.equals("[M+H-H2O]+"))ionization=-17.00219;
							else if(ion.equals("[M+Na]+"))ionization=22.98922;
							else if(ion.equals("[M+NH4]+"))ionization=18.03437;
							else if(ion.equals("[M+K]+"))ionization=39.098301;
							else{
								System.err.println("Ionization unknown");
								continue;
							}
						}
						if(line.startsWith(">charge")&&!intrinsicCharged){
							String charge=line.substring(8);
							if(charge.equals("1"))ionization=1.00728;				
						}
						if(line.matches("\\d*[\\.\\d*]+\\s\\d*[\\.\\d*]+")){
							String l[]=line.split("\\s");
							double mz=Double.parseDouble(l[0]);
							double intensity=Double.parseDouble(l[1]);
							if(mz-0.1<mf.getMass()+ionization)peaks.add(new Peak(mz, intensity));
							maxIntensity=Math.max(maxIntensity,intensity);
						}
					}
					
					if(peaks.isEmpty())continue;
					
					Collections.sort(peaks);
					List<List<Peak>> peakGroups=new ArrayList<List<Peak>>();
					List<Peak> currentListPeaks=new ArrayList<Peak>();
					peakGroups.add(currentListPeaks);
					
					for(int i=0;i<peaks.size();i++){
						if(i==0){
							currentListPeaks.add(peaks.get(i));
						}else{
							double error=Utils.getAbsoluteErrorForMass(peaks.get(i).mz, 10, 2);
							if(Math.abs(peaks.get(i).mz-peaks.get(i-1).mz)>error){
								currentListPeaks=new ArrayList<Peak>();
								peakGroups.add(currentListPeaks);	
							}
							currentListPeaks.add(peaks.get(i));
						}
					}
					
					int i=0;
					boolean parentPeakFound=false;
					for(List<Peak> peakGroup:peakGroups){
						Peak peakWithMaxInt=null;
						for(Peak p:peakGroup){
							if(peakWithMaxInt==null||peakWithMaxInt.intensity<p.intensity){
								peakWithMaxInt=p;
							}
						}
						Vertex currentVertex = new Vertex("v"+i++);
						if(Math.abs(peakWithMaxInt.mz-ionization-mf.getMass())<0.2){
							currentVertex.getProperties().put("label", mf.toString()+"\\nexact mass "+(peakWithMaxInt.mz-ionization)+" Da, "+peakWithMaxInt.intensity+"%");
							parentPeakFound=true;
						}else{
							currentVertex.getProperties().put("label", "exact mass "+(peakWithMaxInt.mz-ionization)+" Da, "+peakWithMaxInt.intensity+"%");
						}
						g.getVertices().add(currentVertex);
					}
					if(!parentPeakFound){
						if(parentMass==null){
							parentMass=mf.getMass()+ionization;
						}
						System.out.println("parent peak added");
						Vertex currentVertex = new Vertex("v"+i++);
						currentVertex.getProperties().put("label", mf.toString()+"\\nexact mass "+(parentMass-ionization)+" Da, 1%");
						g.getVertices().add(currentVertex);
					}
//					if(!parentPeakFound){
//						System.err.println("no parent peak found for "+f.getName());
//						continue;
//					}
					if(g.getVertices().size()<=1){
						System.err.println("no peaks found for "+f.getName());
						continue;
					}
					if(storeGraphs)result.put(f,g);
					else result.put(f,null);
					br.close();
					System.out.println("Graph for "+f.getName()+" is read.");
				}catch(IOException e){
					System.err.println("Could not read file "+f.getAbsolutePath());
				}
			}
		}

		return result;
	}
	
	public static double getIntensityFromVertex(Vertex v){
		String label=v.getProperties().get("label");
		int start=label.indexOf("Da,")+4;
		int end=label.indexOf("%");
		return Double.parseDouble(label.substring(start,end).trim());
		
	}
	
	public static void getVerticesOfSpecialDepth(DefaultDecoyTree g, DecoyTreeVertex root, int depth, Set<DecoyTreeVertex> vertices){
		if(depth==0){
			vertices.add(root);
		}else{
			for(Edge e:g.getOutgoingEdgesFor(root)){
				getVerticesOfSpecialDepth(g, g.getVertex(e.getTail()), depth-1, vertices);
			}
		}
		
	}
	
	public static double getAbsoluteErrorForMass(double mass, double ppm, double ae){
		return Math.max(ppm*1e-6*mass, ae*1e-3);
	}

}

class Peak implements Comparable<Peak>{
	double mz;
	double intensity;
	
	public Peak(double mz, double intensity){
		this.mz=mz;
		this.intensity=intensity;
	}

	@Override
	public int compareTo(Peak o) {
		return Double.compare(mz, o.mz);
	}

}