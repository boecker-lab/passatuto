package de.unijena.bioinf.decoy.model;

import java.util.Iterator;
import java.util.Map;

import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.babelms.dot.Vertex;
import de.unijena.bioinf.deocy.Utils;

public class DecoyTreeVertex extends Vertex{

	public MolecularFormula mf;
	public double mz;
	public double intensity;
	
//	public DecoyTreeVertex(String name) {
//		super(name);
//	}
	
	public DecoyTreeVertex(Vertex v) {
		super(v);
		mf=this.getMolecularFormulaFromLabel();
		Double mz=this.getExactMassFromLabel();
		this.mz=(mz==null)?mf.getMass():mz;
		intensity=Utils.getIntensityFromVertex(v);
	}
	
	public String toString() {
		final StringBuilder b = new StringBuilder();
		b.append(name);
		b.append(" [");
		b.append("label=\""+mf.formatByKerstin()+"\\n"+mf.getMass()+" Da, "+intensity+" %\"");
		b.append("];");
		return b.toString();
	}
	
	public MolecularFormula getMolecularFormulaFromLabel(){
		if(!getProperties().get("label").contains("\\n"))return null;
		return MolecularFormula.parse(getProperties().get("label").replaceAll("\\\\n.*",""));
	}
	
	public static MolecularFormula getMolecularFormulaFromLabel(Vertex v){
		if(!v.getProperties().get("label").contains("\\n"))return null;
		return MolecularFormula.parse(v.getProperties().get("label").replaceAll("\\\\n.*",""));
	}
	
	public Double getExactMassFromLabel(){
		String l=getProperties().get("label");
		if(!l.contains("exact mass")) return null;
		int start=l.indexOf("exact mass")+11;
		int end=l.indexOf("Da");
		return Double.parseDouble(l.substring(start,end).trim());
	}
	
	public static Double getExactMassFromLabel(Vertex v){
		String l=v.getProperties().get("label");
		if(!l.contains("exact mass")) return null;
		int start=l.indexOf("exact mass")+11;
		int end=l.indexOf("Da");
		return Double.parseDouble(l.substring(start,end).trim());
	}

}
