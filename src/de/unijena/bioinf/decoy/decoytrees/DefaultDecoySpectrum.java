package de.unijena.bioinf.decoy.decoytrees;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.babelms.dot.Graph;
import de.unijena.bioinf.babelms.dot.Vertex;
import de.unijena.bioinf.decoy.model.DecoyTreeEdge;
import de.unijena.bioinf.decoy.model.DecoyTreeVertex;
import de.unijena.bioinf.decoy.model.MassBank;
import de.unijena.bioinf.decoy.model.VertexComparator;
import de.unijena.bioinf.deocy.Utils;

public class DefaultDecoySpectrum implements DecoySpectrum {
	MolecularFormula parent;
	List<double[]> peaks;
	
	public DefaultDecoySpectrum(){}

	public DefaultDecoySpectrum(MolecularFormula parent, List<double[]> peaks){
		this.parent=parent;
		this.peaks=peaks;
	}
	
	@Override
	public MolecularFormula getParent() {
		return parent;
	}

	@Override
	public List<double[]> getPeaks() {
		return peaks;
	}
	
	public String writeToString() {
		final StringWriter strw = new StringWriter();
		try {
			writeAsDot(strw);
		} catch (IOException e) {
			e.printStackTrace();
		}
		return strw.toString();
	}

	public void writeAsDot(Writer writer) throws IOException {	
	}
	
	public void writeAsDot(File file) throws IOException {
	}
	
	public void writeAsMassbank(File outputFolderMassbank, List<double[]> peaks, MolecularFormula parent, double massParentIon, String acc, String recordTitle, String instrumentName, String inchi, Boolean positive, boolean containsPPDiff) throws IOException {
		MassBank mbReal=new MassBank(peaks, parent, massParentIon, acc, recordTitle, instrumentName, inchi, positive, containsPPDiff);
		mbReal.writeSpectrumAsMassbank(outputFolderMassbank);
	}
	
	public void writeAsMS(File file) throws IOException {
		double hAdduct=1.00728;
		final BufferedWriter bw = new BufferedWriter(new FileWriter(file));
		bw.write(">compound "+file.getName().substring(0,file.getName().lastIndexOf("."))+"\n");
		bw.write(">formula "+getParent()+"\n");
		bw.write(">parentmass "+(getParent().getMass()+hAdduct)+"\n");
		bw.write(">ionization [M+H]+\n\n");

		List<double[]> peaks=getPeaks();		

		double tic=0;
		for (double[] p : peaks) {
			tic+=p[1];
		}

		bw.write(">collision 0\n");
		bw.write(">tic "+tic+"\n");
		for (double[] p : peaks) {
			bw.write((p[0]+hAdduct)+" "+p[1]+"\n");
		}
		bw.flush();
		bw.close();
	}
	
	public boolean isTree(){
		return false;
	}

	@Override
	public boolean equalsGraph(Graph g) {
		for(Vertex v:g.getVertices()){
			boolean found=false;
			DecoyTreeVertex dtv=new DecoyTreeVertex(v);
			Double mz=dtv.mz;
			Double intensity=dtv.intensity;			
			for(double[] d:this.peaks){
				if(mz.equals(d[0])&&intensity.equals(d[1])){
					found=true;
					break;
				}
			}
			if(!found)return false;
		}
		for(double[] d:this.peaks){
			boolean found=false;
			for(Vertex v:g.getVertices()){
				DecoyTreeVertex dtv=new DecoyTreeVertex(v);
				Double mz=dtv.mz;
				Double intensity=dtv.intensity;
				if(mz.equals(d[0])&&intensity.equals(d[1])){
					found=true;
					break;
				}
			}	
			if(!found)return false;
		}
		
		return true;
	}
	

}
