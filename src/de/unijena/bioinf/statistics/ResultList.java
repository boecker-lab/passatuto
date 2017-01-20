package de.unijena.bioinf.statistics;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import de.unijena.bioinf.decoy.model.MassBank;
import de.unijena.bioinf.deocy.Utils;

public class ResultList{
	Query query;
	List<Result> results;

	public ResultList(){		
		this.results=new ArrayList<Result>();
	}

	public ResultList(ResultList rl){
		this.query=rl.query;
		this.results=new ArrayList<Result>();
		this.results.addAll(rl.results);
	}

	public static List<ResultList> getMergedResults(File[] files, Map<String, MassBank> dbFiles) throws Exception{

		List<ResultList>[] tmp= new ArrayList[files.length];

		int i=0;
		for(File f:files){
			tmp[i++]=getResults(f, dbFiles);	
		}

		List<ResultList> results=ResultList.mergeResults(tmp);

		return results;
	}

	public static List<ResultList> getResults(File f, Map<String, MassBank> dbFiles) throws Exception{
		if (f.getName().endsWith(".txt")){
			return getResultsTXT(f, dbFiles);
		}else{
			return getResultsCSV(f, dbFiles);
		}
	}

	public static List<ResultList> getResultsCSV(File f, Map<String, MassBank> dbFiles) throws Exception{
		List<ResultList> results=new ArrayList<ResultList>();

		BufferedReader br=new BufferedReader(new FileReader(f));		
		br.readLine();br.readLine();br.readLine();
		String line=br.readLine();
		List<String> resultIDs=Arrays.asList(line.split(",",Integer.MAX_VALUE));
		
		List<MassBank> resultIDsMassBank=new ArrayList<MassBank>();
		for(String s:resultIDs)resultIDsMassBank.add(dbFiles.get(s));

		ResultList rl=null;
		while((line=br.readLine())!=null){
			String l[]=line.split(",");
			rl=new ResultList();
			results.add(rl);
			Query q=new Query(l[0]);				
			q.massbank=dbFiles.get(q.queryID);
			rl.query=q;
			for(int i=1;i<l.length;i++){
				String[] res=l[i].split(" ");
				double score=Double.parseDouble(res[0]);
				if(score>0){
					Result r=new Result(resultIDs.get(i),Integer.parseInt(res[1]), score);
					r.massbank=resultIDsMassBank.get(i);
					rl.results.add(r);
					r.isTrueMatch=rl.query.massbank.inchi.equals(r.massbank.inchi);
				}
			}
			Collections.sort(rl.results);
		}
		br.close();
		return results;
	}

	public static List<ResultList> getResultsTXT(File f, Map<String, MassBank> dbFiles) throws Exception{
		List<ResultList> results=new ArrayList<ResultList>();

		BufferedReader br=new BufferedReader(new FileReader(f));
		String line;
		ResultList rl=null;
		while((line=br.readLine())!=null){
			if(line.startsWith("query")){
				rl=new ResultList();
				results.add(rl);
				Query q=new Query(line.split("\t")[1]);				
				q.massbank=dbFiles.get(q.queryID);
				rl.query=q;
			}else if(line.startsWith("result")){
				Result r=new Result(line.split("\t"));
				r.massbank=dbFiles.get(r.resultID);
				rl.results.add(r);
				r.isTrueMatch=rl.query.massbank.inchi.equals(r.massbank.inchi);
			}
		}
		br.close();


		return results;
	}

	public static List<ResultList> mergeResults(List<ResultList>[] resultLists) throws Exception{
		List<ResultList> result=new ArrayList<ResultList>();
		for(List<ResultList> resultList:resultLists){
			for(ResultList rl_tmp_1:resultList){
				ResultList found=null;
				for(ResultList rl_tmp_2:result){
					if(rl_tmp_1.query.queryID.equals(rl_tmp_2.query.queryID)){
						found=rl_tmp_2;
						break;
					}
				}
				if(found!=null){
					found.results.addAll(rl_tmp_1.results);				
				}else{
					result.add(new ResultList(rl_tmp_1));
				}
			}
		}

		for(ResultList resultList:result){
			Collections.sort(resultList.results);
		}
		return result;
	}

	public void filterResultsByMass(boolean onlySame, double ppm, double ae){
		Iterator<Result> it=results.iterator();
		double mass=query.massbank.mf.getMass();
		double error=Utils.getAbsoluteErrorForMass(mass, ppm, ae);
		while(it.hasNext()){
			Result r=it.next();
			boolean isNotSame=Math.abs(r.massbank.mf.getMass()-mass)>error;
			if((onlySame&&isNotSame)||(!onlySame&&!isNotSame))it.remove();
		}	
	}

	public void filterResultsByInChi(boolean retainSame){
		String inchiQuery=query.massbank.inchi;
		Iterator<Result> it=results.iterator();	
		while(it.hasNext()){
			Result r=it.next();
			String inchiResult=r.massbank.inchi;
			boolean isNotSame=inchiQuery==null||inchiResult==null||inchiQuery.isEmpty()||inchiResult.isEmpty()||!inchiQuery.equals(inchiResult);
			if((retainSame&&isNotSame)||(!retainSame&&!isNotSame))it.remove();
		}	
	}
}