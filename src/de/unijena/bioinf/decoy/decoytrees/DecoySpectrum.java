package de.unijena.bioinf.decoy.decoytrees;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringWriter;
import java.io.Writer;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;

import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.babelms.dot.Graph;
import de.unijena.bioinf.babelms.dot.Vertex;
import de.unijena.bioinf.decoy.model.DecoyTreeEdge;
import de.unijena.bioinf.decoy.model.DecoyTreeVertex;
import de.unijena.bioinf.decoy.model.MassBank;
import de.unijena.bioinf.decoy.model.VertexComparator;

public interface DecoySpectrum {

	public void writeAsDot(File file) throws IOException;
	
	public boolean isTree();
	
	public void writeAsMS(File file) throws IOException;
	
	public MolecularFormula getParent();
	
	public List<double[]> getPeaks();

	public void writeAsMassbank(File outputFolderMassbank,
			List<double[]> peaks, MolecularFormula parent, double massParentIon, String acc,
			String string, String instrumentName, String inchi,
			Boolean positive, boolean contains) throws IOException ;
	
	public boolean equalsGraph(Graph g);

	

}
