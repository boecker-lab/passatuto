package de.unijena.bioinf.decoy.decoytrees;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringWriter;
import java.io.Writer;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;

import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.babelms.dot.Edge;
import de.unijena.bioinf.babelms.dot.Graph;
import de.unijena.bioinf.babelms.dot.Vertex;
import de.unijena.bioinf.decoy.model.DecoyTreeEdge;
import de.unijena.bioinf.decoy.model.DecoyTreeVertex;
import de.unijena.bioinf.decoy.model.MassBank;
import de.unijena.bioinf.decoy.model.VertexComparator;
import de.unijena.bioinf.deocy.Utils;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction.POSTPROCESS;

public class DefaultDecoyTree extends DefaultDecoySpectrum implements Cloneable,DecoySpectrum{

	private final ArrayList<DecoyTreeVertex> vertices;
	private final ArrayList<DecoyTreeEdge> edges;
	public boolean isPositive=true;

	public DefaultDecoyTree() {
		this.vertices = new ArrayList<DecoyTreeVertex>();
		this.edges = new ArrayList<DecoyTreeEdge>();
	}

	public DefaultDecoyTree(Graph g) {
		this.vertices = new ArrayList<DecoyTreeVertex>(g.getVertices().size());
		this.edges = new ArrayList<DecoyTreeEdge>(g.getEdges().size());
		for (Vertex v : g.getVertices()) vertices.add(new DecoyTreeVertex(v));
		for (Edge e : g.getEdges()) edges.add(new DecoyTreeEdge(e));
	}
	
	public DefaultDecoyTree(DefaultDecoyTree g) {
		this.vertices = new ArrayList<DecoyTreeVertex>(g.getVertices().size());
		this.edges = new ArrayList<DecoyTreeEdge>(g.getEdges().size());
		for (Vertex v : g.getVertices()) vertices.add(new DecoyTreeVertex(v));
		for (Edge e : g.getEdges()) edges.add(new DecoyTreeEdge(e));
	}
	
	public void writeAsDot(Writer writer) throws IOException {
		final BufferedWriter bw = new BufferedWriter(writer);
		bw.write("strict digraph {\n");
		for (DecoyTreeVertex v : vertices) {
			bw.write(v.toString());
			bw.write('\n');
		}
		for (DecoyTreeEdge e : edges) {
			bw.write(e.toString());
			bw.write('\n');
		}
		bw.write("}");
		bw.flush();
	}
	
	public void writeAsDot(File f) throws IOException {
		if(!f.getParentFile().exists())f.getParentFile().mkdirs();
		final BufferedWriter bw = new BufferedWriter(new FileWriter(f));
		writeAsDot(bw);
	}
	
	

	public List<DecoyTreeEdge> getIncommingEdgesFor(DecoyTreeVertex vertex) {
		return getIncommingEdgesFor(vertex.getName());
	}
	
	public MolecularFormula getParent(){
		return getRoot().mf;
	}
	
	public boolean isTree(){
		for (DecoyTreeVertex v : vertices) {
			if(v.mf==null) return false;		
		}
		return true;
	}
	
	public DecoyTreeVertex getRoot() {
		final HashSet<String> vertices = new HashSet<String>();
		for (DecoyTreeEdge e : edges) {
			vertices.add(e.getHead());
		}
		for (DecoyTreeEdge e : edges) {
			vertices.remove(e.getTail());
		}
		if (vertices.isEmpty()){
			DecoyTreeVertex result=null;
			for (DecoyTreeVertex v : this.vertices) {
				if(v.mf!=null&&(result==null||Math.abs(result.mf.getMass()-result.mz)>Math.abs(v.mf.getMass()-v.mz)))result=v;
			}
			return result;
		}
		final String s = vertices.iterator().next();
		return getVertex(s);

	}

	public DecoyTreeVertex getVertex(String name) {
		for (DecoyTreeVertex u : vertices)
			if (u.getName().equals(name)) return u;
		return null;
	}

	public List<DecoyTreeEdge> getIncommingEdgesFor(String vertex) {
		final ArrayList<DecoyTreeEdge> neighbours = new ArrayList<DecoyTreeEdge>();
		for (DecoyTreeEdge e : edges) {
			if (e.getTail().equals(vertex)) neighbours.add(e);
		}
		return neighbours;
	}

	public List<DecoyTreeEdge> getOutgoingEdgesFor(Vertex vertex) {
		return getOutgoingEdgesFor(vertex.getName());
	}

	public List<DecoyTreeEdge> getOutgoingEdgesFor(String vertex) {
		final ArrayList<DecoyTreeEdge> neighbours = new ArrayList<DecoyTreeEdge>();
		for (DecoyTreeEdge e : edges) {
			if (e.getHead().equals(vertex)) neighbours.add(e);
		}
		return neighbours;
	}

	public DecoyTreeEdge getEdgeFor(DecoyTreeVertex u, DecoyTreeVertex v) {
		return getEdgeFor(u.getName(), v.getName());
	}

	public DecoyTreeEdge getEdgeFor(String u, String v) {
		for (DecoyTreeEdge e : edges) {
			if (e.getHead().equals(u) && e.getTail().equals(v)) return e;
		}
		return null;
	}

	public ArrayList<DecoyTreeVertex> getVertices() {
		return vertices;
	}

	public ArrayList<DecoyTreeEdge> getEdges() {
		return edges;
	}
	
	public List<double[]> getPeaks(){
		return getPeaks(Double.POSITIVE_INFINITY);
	}
	
	public List<double[]> getOriginalPeaks(){
		return getOriginalPeaks(Double.POSITIVE_INFINITY);
	}
	
	public List<double[]> getPeaks(double maxMass){
		List<double[]> peaks=new ArrayList<double[]>();
		for (DecoyTreeVertex v : getVertices()) {
			if(v.mz<maxMass)peaks.add(new double[]{v.mf.getMass(),v.intensity});
		}
		return peaks;
	}
	
	public List<double[]> getOriginalPeaks(double maxMass){
		List<double[]> peaks=new ArrayList<double[]>();
		for (DecoyTreeVertex v : getVertices()) {
			if(v.mz<maxMass)peaks.add(new double[]{v.mz,v.intensity});
		}
		return peaks;
	}
	
	public List<double[]> getPeaksIncludingMax(double maxMass){
		List<double[]> peaks=new ArrayList<double[]>();
		for (DecoyTreeVertex v : getVertices()) {		
			if(v.mz<=maxMass)peaks.add(new double[]{v.mz,v.intensity});
		}
		return peaks;
	}
	
	public static List<double[]> getPeaksIncludingMax(Graph g, double maxMass){
		List<double[]> peaks=new ArrayList<double[]>();
		for (Vertex v : g.getVertices()) {
			MolecularFormula mf=DecoyTreeVertex.getMolecularFormulaFromLabel(v);
			Double mz=DecoyTreeVertex.getExactMassFromLabel(v);
			mz=(mz==null)?mf.getMass():mz;
			double intensity=Utils.getIntensityFromVertex(v);
			
			if(mz<=maxMass)peaks.add(new double[]{mz,intensity});
		}
		return peaks;
	}
	
	@Override
	public boolean equalsGraph(Graph g) {
		for(Vertex v:g.getVertices()){
			boolean found=false;
			DecoyTreeVertex dtv=new DecoyTreeVertex(v);
			Double mz=dtv.mz;
			Double intensity=dtv.intensity;			
			for(DecoyTreeVertex dtv2:this.vertices){
				if(mz.equals(dtv2.getExactMassFromLabel())&&intensity.equals(dtv2.intensity)){
					found=true;
					break;
				}
			}
			if(!found)return false;
		}
		for(DecoyTreeVertex dtv2:this.vertices){
			boolean found=false;
			for(Vertex v:g.getVertices()){
				DecoyTreeVertex dtv=new DecoyTreeVertex(v);
				Double mz=dtv.mz;
				Double intensity=dtv.intensity;
				if(mz.equals(dtv2.getExactMassFromLabel())&&intensity.equals(dtv2.intensity)){
					found=true;
					break;
				}
			}	
			if(!found)return false;
		}
		return true;
	}
}
