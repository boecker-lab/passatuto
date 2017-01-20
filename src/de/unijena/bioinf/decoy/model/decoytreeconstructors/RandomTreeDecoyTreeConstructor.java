package de.unijena.bioinf.decoy.model.decoytreeconstructors;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;

import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.babelms.dot.Edge;
import de.unijena.bioinf.babelms.dot.Graph;
import de.unijena.bioinf.babelms.dot.Vertex;
import de.unijena.bioinf.decoy.decoytrees.DefaultDecoyTree;
import de.unijena.bioinf.decoy.model.DecoyTreeEdge;
import de.unijena.bioinf.decoy.model.DecoyTreeVertex;
import de.unijena.bioinf.deocy.Utils;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction;

public class RandomTreeDecoyTreeConstructor extends DefaultDecoyTreeConstructor implements DecoySpectrumConstructor{

	static public final String methodNameLong="RandomTree";	
	Random r = new Random();
	double minMass;

	public RandomTreeDecoyTreeConstructor(ParametersDecoySpectrumConstruction p){
		super(p);
		this.minMass=p.getMinMass();
	}

	public String getMethodNameLong(){
		return methodNameLong;
	}

//	public DefaultDecoyTree getDecoySpectrum(){
//		List<Graph> allGraphs=new ArrayList<Graph>();
//		allGraphs.addAll(p.getAllGraphsInOriginal().values());
//		List<MolecularFormula> mfs=new ArrayList<MolecularFormula>();
//
//		for(Graph t:allGraphs){
//			DefaultDecoyTree g=new DefaultDecoyTree(t);
//			for(DecoyTreeEdge e:g.getEdges())mfs.add(e.mf);
//		}
//
//		DefaultDecoyTree decoyTree=new DefaultDecoyTree(inputGraph);
//		Random r=new Random();
//		List<DecoyTreeEdge> currentEdges=decoyTree.getOutgoingEdgesFor(decoyTree.getRoot());
//		while(!currentEdges.isEmpty()){
//			DecoyTreeEdge e=currentEdges.get(0);
//			DecoyTreeVertex vH=null;
//			List<MolecularFormula> mfsTmp=new ArrayList<MolecularFormula>();
//			while(mfsTmp.isEmpty()){
//				if(vH==null){
//					vH=decoyTree.getVertex(e.getHead());
//				}else{
//					System.out.print("rearrangement "+e.toString());
//					vH=decoyTree.getVertex(decoyTree.getIncommingEdgesFor(vH.getName()).get(0).getHead());
//					System.out.println(" to "+vH.getName());
//				}
//				MolecularFormula mfH=vH.mf;
//				for(MolecularFormula mf:mfs)if(mfH.isSubtractable(mf))mfsTmp.add(mf);
//
//			}				
//			MolecularFormula mf=mfsTmp.get(r.nextInt(mfsTmp.size()));
//			e.mf=mf;
//			e.setHead(vH.getName());
//			decoyTree.getVertex(e.getTail()).mf=vH.mf.subtract(mf);					
//			currentEdges.remove(e);
//			for(DecoyTreeEdge newE:decoyTree.getOutgoingEdgesFor(decoyTree.getVertex(e.getTail())))currentEdges.add(newE);
//		}
//
//		if(!currentEdges.isEmpty())return null;
//
//		return decoyTree;
//	}
	
	public DefaultDecoyTree getDecoySpectrum(){
		List<Graph> allGraphs=new ArrayList<Graph>();
		allGraphs.addAll(p.getAllGraphsInOriginal().values());
		List<MolecularFormula> mfs=new ArrayList<MolecularFormula>();

		for(Graph t:allGraphs){
			for(Edge e:t.getEdges()){
				MolecularFormula mf=MolecularFormula.parse(e.getProperties().get("label").replaceAll("\\\\n.*",""));
				mfs.add(mf);
			}
		}

		DefaultDecoyTree decoyTree=new DefaultDecoyTree(inputGraph);
		Random r=new Random();
		List<DecoyTreeEdge> currentEdges=decoyTree.getOutgoingEdgesFor(decoyTree.getRoot());
		while(!currentEdges.isEmpty()){
			DecoyTreeEdge e=currentEdges.get(0);
			DecoyTreeVertex vH=null;
			List<MolecularFormula> mfsTmp=new ArrayList<MolecularFormula>();
			while(mfsTmp.isEmpty()){
				if(vH==null){
					vH=decoyTree.getVertex(e.getHead());
				}else{
					System.out.print("rearrangement "+e.toString());
					vH=decoyTree.getVertex(decoyTree.getIncommingEdgesFor(vH.getName()).get(0).getHead());
					System.out.println(" to "+vH.getName());
				}
				MolecularFormula mfH=vH.mf;
				for(MolecularFormula mf:mfs)if(mfH.isSubtractable(mf))mfsTmp.add(mf);

			}				
			MolecularFormula mf=mfsTmp.get(r.nextInt(mfsTmp.size()));
			e.mf=mf;
			e.setHead(vH.getName());
			decoyTree.getVertex(e.getTail()).mf=vH.mf.subtract(mf);					
			currentEdges.remove(e);
			for(DecoyTreeEdge newE:decoyTree.getOutgoingEdgesFor(decoyTree.getVertex(e.getTail())))currentEdges.add(newE);
		}

		if(!currentEdges.isEmpty())return null;

		return decoyTree;
	}
}
