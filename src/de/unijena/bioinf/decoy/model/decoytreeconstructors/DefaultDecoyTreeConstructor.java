package de.unijena.bioinf.decoy.model.decoytreeconstructors;

import de.unijena.bioinf.babelms.dot.Graph;
import de.unijena.bioinf.decoy.decoytrees.DecoySpectrum;
import de.unijena.bioinf.decoy.decoytrees.DefaultDecoySpectrum;
import de.unijena.bioinf.decoy.decoytrees.DefaultDecoyTree;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction;

public class DefaultDecoyTreeConstructor implements DecoySpectrumConstructor{
	
	static public final String methodNameLong="Original";	
	public final ParametersDecoySpectrumConstruction p;
	public Graph originalGraph=null;
	public DefaultDecoyTree inputGraph=null;
	
	public DefaultDecoyTreeConstructor(ParametersDecoySpectrumConstruction p){
		this.p=p;
	}
	
	public String getMethodNameLong(){
		return methodNameLong;
	}
	
	public void setOriginalTree(Graph g){
		this.originalGraph=g;
		this.inputGraph=new DefaultDecoyTree(g);
	}
	
	public DecoySpectrum getDecoySpectrum(){
		if(originalGraph.getRoot()==null)return new DefaultDecoySpectrum(inputGraph.getParent(), inputGraph.getOriginalPeaks());
		return new DefaultDecoyTree(inputGraph);
	}
	
}
