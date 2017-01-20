package de.unijena.bioinf.decoy.model.decoytreeconstructors;

import java.util.ArrayList;
import java.util.List;

import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.babelms.dot.Graph;
import de.unijena.bioinf.decoy.decoytrees.DecoySpectrum;
import de.unijena.bioinf.decoy.decoytrees.DefaultDecoySpectrum;
import de.unijena.bioinf.decoy.decoytrees.DefaultDecoyTree;
import de.unijena.bioinf.decoy.model.DecoyTreeVertex;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction;

public class DefaultDecoyFileConstructor extends DefaultDecoyTreeConstructor implements DecoySpectrumConstructor{
	
	static public final String methodNameLong="Original";	
	
	public DefaultDecoyFileConstructor(ParametersDecoySpectrumConstruction p){
		super(p);
	}
	
	public String getMethodNameLong(){
		return methodNameLong;
	}
	
	public void setOriginalTree(Graph g){
		this.inputGraph=new DefaultDecoyTree(g);
	}
	
	public DecoySpectrum getDecoySpectrum(){
		MolecularFormula root=null;
		List<double[]> peaks=new ArrayList<double[]>();
		for(DecoyTreeVertex v:inputGraph.getVertices()){
			peaks.add(new double[]{v.mz,v.intensity});
			if(v.mf!=null&&(root==null||root.getMass()<v.mf.getMass()))root=v.mf;
		}
		return new DefaultDecoySpectrum(root, peaks);		
	}
	
}
