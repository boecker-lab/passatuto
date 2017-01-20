package de.unijena.bioinf.statistics;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import de.unijena.bioinf.decoy.model.MassBank;

public class MainOld {

	static String searchMethod="MassBank";	
//	static File dbFolder=new File("U:\\MSBlast\\DB");
	static File statisticsFolder=new File("U:\\MSBlast\\statisticsOld\\"+searchMethod);
	
	static String addOutputFile="Pooled";
	
	static String addOriginal="";
	static String sepOriginal=".";
	static String fileEndingOriginal="csv";
	
	static String addDecoy="";
	static String sepDecoy="";
	static String fileEndingDecoy="";
	
	static int numberOfDatabases=1;

	public static void main(String args[]) throws Exception{
		File[] inputFilesOriginal=new File[]{
				new File("U:\\MSBlast\\searchResults\\"+searchMethod+"\\pos\\MPOriginal-APOriginal"+addOriginal+sepOriginal+fileEndingOriginal)
		};
		File[] inputFilesDecoy=new File[]{
//				new File("U:\\MSBlast\\searchResults\\"+searchMethod+"\\pos\\MPOriginal-APConditionalPeaks"+addDecoy+sepDecoy+fileEndingDecoy),
//				new File("U:\\MSBlast\\searchResults\\"+searchMethod+"\\pos\\MPOriginal-APRandomPeaks"+addDecoy+sepDecoy+fileEndingDecoy),
				new File("U:\\MSBlast\\searchResults\\"+searchMethod+"\\pos\\MPOriginal-APReroot"+addDecoy+sepDecoy+fileEndingDecoy)		
		};


		File[] comparisonFiles=new File[inputFilesOriginal.length+inputFilesDecoy.length];
		int k=0;
		for(File f:inputFilesOriginal){
			if(f.isDirectory()){
				comparisonFiles[k++]=f.listFiles()[0];
			}else{
				comparisonFiles[k++]=f;
			}
		}
		for(File f:inputFilesDecoy){
			if(f.isDirectory()){
				comparisonFiles[k++]=f.listFiles()[0];
			}else{
				comparisonFiles[k++]=f;
			}
		}

		Map<String, MassBank> dbFilesOriginal=getAllDBFiles(comparisonFiles);
		List<ResultList>[] resultsOriginal=new List[inputFilesOriginal.length];
		for(int i=0;i<inputFilesOriginal.length;i++){
			File f=inputFilesOriginal[i];
			resultsOriginal[i]=ResultList.getResults(f, dbFilesOriginal);
		}
		List<List<ResultList>>[] resultsDecoy=new List[inputFilesDecoy.length];
		for(int i=0;i<inputFilesDecoy.length;i++){			
			File folder=inputFilesDecoy[i];
			resultsDecoy[i]=new ArrayList<List<ResultList>>();			
			if(folder.isDirectory()){
				List<File> files=Arrays.asList(folder.listFiles());
				Collections.sort(files);
				files=files.subList(0, Math.min(numberOfDatabases,files.size()));
				for(File f:files){
//					resultsDecoy[i].add(ResultList.getResults(f, dbFilesOriginal));
					if(resultsDecoy[i].size()==0){
						resultsDecoy[i].add(ResultList.getResults(f, dbFilesOriginal));
					}else{
						resultsDecoy[i].get(0).addAll(ResultList.getResults(f, dbFilesOriginal));
					}
				}
			}else{
				resultsDecoy[i].add(ResultList.getResults(folder, dbFilesOriginal));
			}
		}
		List<List<ResultList>>[] resultsAll=new List[resultsOriginal.length+resultsDecoy.length];
		int i=0;
		for(List<ResultList> r:resultsOriginal){
			List<List<ResultList>> tmp=new ArrayList<List<ResultList>>();
			tmp.add(r);
			resultsAll[i++]=tmp;
		}
		for(List<List<ResultList>> r:resultsDecoy)resultsAll[i++]=r;


		if(!statisticsFolder.exists())statisticsFolder.mkdirs();
		Plot.writeScoreDistributionOfTopRank(statisticsFolder, addOutputFile, resultsAll, 0.05);		
		Plot.writeEstimatedQValueVSCalculatedQValue(statisticsFolder, addOutputFile, resultsOriginal, resultsDecoy);

		Plot.writeRankVSScore(statisticsFolder, addOutputFile, resultsAll, 10);
		Plot.writeRankVSPercentage(statisticsFolder, addOutputFile, resultsOriginal, resultsDecoy, 10);

		Plot.writeRankVSPercentage(statisticsFolder, addOutputFile+"OnlyDifferentInchi", getFilteredResultListsByInchi(resultsOriginal, false), getFilteredResultListsByInchi(resultsDecoy, false), 10);
		
		Plot.writeEstimatedQValueVSCalculatedQValue(statisticsFolder, addOutputFile+"OnlySameMass", getFilteredResultListsByMass(resultsOriginal, true), getFilteredResultListsByMass(resultsDecoy, true));

	}

	public static List[] getFilteredResultListsByInchi(List[] results, boolean retainSame){
		List[] result=new ArrayList[results.length];
		for(int i=0;i<results.length;i++){			
			for(Object rl:results[i]){				
				if(rl.getClass().equals(ResultList.class)){
					if(result[i]==null)result[i]=new ArrayList<ResultList>();
					ResultList newRL=new ResultList((ResultList)rl);
					newRL.filterResultsByInChi(retainSame);
					if(!newRL.results.isEmpty())result[i].add(newRL);
				}else{
					if(result[i]==null)result[i]=new ArrayList<List<ResultList>>();
					List<ResultList> tmp=new ArrayList<ResultList>();
					for(ResultList rl2:(List<ResultList>) rl){
						ResultList newRL=new ResultList(rl2);
						newRL.filterResultsByInChi(retainSame);
						if(!newRL.results.isEmpty())tmp.add(newRL);
					}
					result[i].add(tmp);
				}
			}
		}
		return result;
	}
	
	public static List[] getFilteredResultListsByMass(List[] results, boolean retainSame){
		List[] result=new ArrayList[results.length];
		for(int i=0;i<results.length;i++){			
			for(Object rl:results[i]){				
				if(rl.getClass().equals(ResultList.class)){
					if(result[i]==null)result[i]=new ArrayList<ResultList>();
					ResultList newRL=new ResultList((ResultList)rl);
					newRL.filterResultsByMass(retainSame,20,5);
					if(!newRL.results.isEmpty())result[i].add(newRL);
				}else{
					if(result[i]==null)result[i]=new ArrayList<List<ResultList>>();
					List<ResultList> tmp=new ArrayList<ResultList>();
					for(ResultList rl2:(List<ResultList>) rl){
						ResultList newRL=new ResultList(rl2);
						newRL.filterResultsByMass(retainSame,20,5);
						if(!newRL.results.isEmpty())tmp.add(newRL);
					}
					result[i].add(tmp);
				}
			}
		}
		return result;
	}

	private static Map<String, MassBank> getAllDBFiles(File[] inputFiles)throws Exception {
		Map<String, MassBank> dbFiles=new HashMap<String,MassBank>();
		Set<File> dbFolders=new HashSet<File>();
		dbFolders.addAll(getFolders(inputFiles,"Query"));
		dbFolders.addAll(getFolders(inputFiles,"DB"));
		getAllDBFiles(dbFolders, dbFiles);
		return dbFiles;
	}

	public static Set<File> getFolders(File[] inputFiles,String queryOrDB) throws Exception{
		Set<File> result=new HashSet<File>();
		for(File f:inputFiles)result.addAll(getFolders(f,queryOrDB));
		return result;
	}

	public static Set<File> getFolders(File inputFile, String queryOrDB) throws Exception{		
		BufferedReader br=new BufferedReader(new FileReader(inputFile));
		String line;
		while((line=br.readLine())!=null){
			if(line.startsWith(queryOrDB)){
				Set<File> result=new HashSet<File>();
				String[] queryFiles=line.substring(queryOrDB.length()+2).split(";");
				for(String f:queryFiles)result.add(new File(f));
				return result;
			}
		}
		return null;
	}

	public static void getAllDBFiles(Set<File> dbFolders, Map<String, MassBank> results) throws Exception{
		for(File dbFolder:dbFolders)getAllDBFiles(dbFolder, results);
	}

	public static void getAllDBFiles(File dbFolder, Map<String, MassBank> results) throws Exception{	
		for(File f:dbFolder.listFiles()){
			if(f.isFile()&&f.getName().endsWith(".txt")){
				results.put(f.getName().replaceAll(".txt",""), new MassBank(f));
				System.out.println(f.getAbsolutePath());
			}else if (f.isDirectory()){
				getAllDBFiles(f,results);
			}
		}
	}
}
