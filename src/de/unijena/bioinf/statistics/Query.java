package de.unijena.bioinf.statistics;

import de.unijena.bioinf.decoy.model.MassBank;

public class Query{
	String queryID;
	MassBank massbank;
	
//	public Query(String[] line){
//		queryID=line[1];		
//	}
	
	public Query(String queryID){
		this.queryID=queryID;		
	}
	
	public Query(Query q){
		this.queryID=q.queryID;
		this.massbank=q.massbank;
	}
	
	public String getDB(){
		return queryID.replaceAll("\\d","");
	}
	
	public String getNumber(){
		return queryID.replaceAll("\\D","");
	}
	
	public String getDataset(){
		return queryID.substring(0,1);
	}
	
}