package de.unijena.bioinf.decoy.model;

import java.util.Iterator;
import java.util.Map;

import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.babelms.dot.Edge;
import de.unijena.bioinf.deocy.Utils;

public class DecoyTreeEdge extends Edge{
	
	public MolecularFormula mf;

	public DecoyTreeEdge(Edge e) {
		super(e);
		mf=this.getMolecularFormulaFromLabel();
	}
	
	public void redirect(){
		String tail=getTail();
		setTail(getHead());
		setHead(tail);
	}
	
    public void setHead(String u) {
        this.u=u;
    }

    public void setTail(String v) {
        this.v=v;
    }
    
	public String toString() {
		final StringBuilder b = new StringBuilder();
        b.append(u).append(" -> ").append(v);
        b.append(" [");
        b.append("label=\""+mf.formatByKerstin()+"\"");
        b.append("];");
        return b.toString();
	}
	
	public MolecularFormula getMolecularFormulaFromLabel(){
		return MolecularFormula.parse(getProperties().get("label").replaceAll("\\\\n.*",""));
	}

}
