package de.unijena.bioinf.decoy.model.decoytreeconstructors;

import de.unijena.bioinf.babelms.dot.Graph;
import de.unijena.bioinf.decoy.decoytrees.DecoySpectrum;

public interface DecoySpectrumConstructor {
	
	static public final String methodNameLong=null;
	
	public DecoySpectrum getDecoySpectrum();
	public void setOriginalTree(Graph g);
	public String getMethodNameLong();

}
