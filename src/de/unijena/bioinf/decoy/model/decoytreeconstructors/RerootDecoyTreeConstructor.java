package de.unijena.bioinf.decoy.model.decoytreeconstructors;

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
import de.unijena.bioinf.babelms.dot.Graph;
import de.unijena.bioinf.babelms.dot.Vertex;
import de.unijena.bioinf.decoy.decoytrees.DefaultDecoyTree;
import de.unijena.bioinf.decoy.model.DecoyTreeEdge;
import de.unijena.bioinf.decoy.model.DecoyTreeVertex;
import de.unijena.bioinf.deocy.Utils;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction;

public class RerootDecoyTreeConstructor extends DefaultDecoyTreeConstructor implements DecoySpectrumConstructor{
	
	static public final String methodNameLong="Reroot";	
	Random r = new Random();
	double minMass;
	
	public RerootDecoyTreeConstructor(ParametersDecoySpectrumConstruction p){
		super(p);
		this.minMass=p.getMinMass();
	}
	
	public String getMethodNameLong(){
		return methodNameLong;
	}
	
	public DefaultDecoyTree getDecoySpectrum(){
		List<DefaultDecoyTree> decoyTrees=new ArrayList<DefaultDecoyTree>();
		MolecularFormula parentMF=inputGraph.getRoot().getMolecularFormulaFromLabel();
		for(Vertex v:inputGraph.getVertices()){
			if(v!=inputGraph.getRoot()){
				DefaultDecoyTree decoyTree=new DefaultDecoyTree(inputGraph);
				DecoyTreeVertex newRoot=decoyTree.getVertex(v.getName());
				double rootIntensity=decoyTree.getRoot().intensity;
//				decoyTree.getRoot().intensity=newRoot.intensity;
				newRoot.intensity=rootIntensity;
				labelRecursive(decoyTree, newRoot, parentMF, null);
				decoyTrees.add(decoyTree);
			}
		}

		DefaultDecoyTree dt=getRandomDecoyTree(decoyTrees,r, minMass);
		
		rearrangeImpossibleNodes(dt, r, minMass);
		
		return dt;
	}
	
	public static double getProbabilityOfRearrangements(DefaultDecoyTree dt, double minMass){
		return 1.0/(getNumberOfRearrangements(dt, minMass)+1);
	}

	public static int getNumberOfRearrangements(DefaultDecoyTree dt, double minMass){
		int rearrangements=0;
		for(DecoyTreeVertex fragment:dt.getVertices()){
			if(fragment.mf.shouldBeRearranged(minMass))rearrangements++;
		}
		return rearrangements;
	}
	
	public static void labelRecursive(DefaultDecoyTree dt, DecoyTreeVertex v, MolecularFormula mf, DecoyTreeEdge excludeEdge){
		v.mf=mf;
		v.getProperties().put("label",mf.toString());
		for(DecoyTreeEdge e:dt.getIncommingEdgesFor(v)){
			if(e!=excludeEdge){
				e.redirect();
			}
		}
		for(DecoyTreeEdge e:dt.getOutgoingEdgesFor(v)){
			if(e!=excludeEdge){
				MolecularFormula mfChild=MolecularFormula.parse(mf.toString());
				MolecularFormula mfLoss=e.getMolecularFormulaFromLabel();
				mfChild=mfChild.subtract(mfLoss);
				labelRecursive(dt, dt.getVertex(e.getTail()), mfChild, e);
			}
		}
	}

	public static void rearrangeImpossibleNodes(DefaultDecoyTree dt, Random r, double minMass) {
		List<DecoyTreeVertex> verticesToRearrange=new ArrayList<DecoyTreeVertex>();
		List<DecoyTreeVertex> allVertices=new ArrayList<DecoyTreeVertex>();
		for(DecoyTreeVertex v:dt.getVertices()){
			allVertices.add(v);
			if(v.mf.shouldBeRearranged(minMass)){
				verticesToRearrange.add(v);
			}
		}


		for(DecoyTreeVertex vtr:verticesToRearrange){
			List<DecoyTreeVertex> possibleParentVertices=new ArrayList<DecoyTreeVertex>();
			List<MolecularFormula> newChildMolecularFormula=new ArrayList<MolecularFormula>();

			List<DecoyTreeEdge> incommingEdges=dt.getIncommingEdgesFor(vtr);
			if(incommingEdges.size()!=0){
				DecoyTreeEdge childE=incommingEdges.get(0);
				for(DecoyTreeVertex v:allVertices){					
					MolecularFormula newVertexMF=v.mf.subtract(childE.mf);								
					if(!newVertexMF.shouldBeRearranged(minMass)){
						//TODO: consider doubled peaks 
						//boolean mfExists=false;
						//for(DecoyTreeVertex equalMF:allVertices){
						//if(equalMF.mf.equals(newVertexMF))mfExists=true;
						//}
						//if(!mfExists){
						possibleParentVertices.add(v);
						newChildMolecularFormula.add(newVertexMF);
						//}
					}
				}
				int choice=r.nextInt(possibleParentVertices.size());
				DecoyTreeVertex parentV=possibleParentVertices.get(choice);
				childE.setHead(parentV.getName());
				vtr.mf=newChildMolecularFormula.get(choice);				
				vtr.getProperties().put("label",vtr.mf.toString());
			}


		}
	}
	
	public static DefaultDecoyTree getRandomDecoyTree(List<DefaultDecoyTree> decoyTrees, Random r, double minMass){
		
		double rand=r.nextDouble();
		
		double sum=0;
		for(DefaultDecoyTree dt:decoyTrees){
			sum+=getProbabilityOfRearrangements(dt, minMass);
		}
		double add=0;
		DefaultDecoyTree result=null;
		for(DefaultDecoyTree dt:decoyTrees){
			add+=getProbabilityOfRearrangements(dt, minMass)/sum;
			if(add>=rand){
				result=dt;
				break;
			}
		}
		return result;
		
	}	


	
}
