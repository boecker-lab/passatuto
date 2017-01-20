package de.unijena.bioinf.statistics;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Map.Entry;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.statistics.DefaultBoxAndWhiskerCategoryDataset;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import de.unijena.bioinf.decoy.model.MassBank;

public class SimilarityMatrix {

	public static enum DatabaseType{Decoy, Original};
	public static enum MatchType{TruePositiveMatch, FalsePositiveMatch, DecoyMatch};

	public List<File> similarityFiles=new ArrayList<File>();
	public Map<String,DatabaseType> typesQuery=new TreeMap<String, DatabaseType>();
	public Map<String,DatabaseType> typesDB=new TreeMap<String, DatabaseType>();
	public List<String> methodsQuery=new ArrayList<String>();
	public List<String> methodsDB=new ArrayList<String>();
	public List<MassBank> massbankQuery=new ArrayList<MassBank>();
	public List<MassBank> massbankDB=new ArrayList<MassBank>();
	double[][] similarity;
	public int[][] commonPeaks;
	public boolean[][] exclusionListDB;
	
	static final String sep=System.getProperty("file.separator");

	public SimilarityMatrix(File similarityFile, TreeMap<String,MassBank> massbank) throws Exception{
		this.similarityFiles.add(similarityFile);		
		BufferedReader br=new BufferedReader(new FileReader(similarityFile));
		String methodsQueryTmp=getMethodFromFilename(br.readLine().replaceAll(";", "").split(" ")[1]);
		String methodsDBTmp=getMethodFromFilename(br.readLine().replaceAll(";", "").split(" ")[1]);
		this.typesQuery.put(methodsQueryTmp,methodsQueryTmp.contains("Original")?DatabaseType.Original:DatabaseType.Decoy);
		this.typesDB.put(methodsDBTmp, methodsDBTmp.contains("Original")?DatabaseType.Original:DatabaseType.Decoy);
		br.readLine();

		if(similarityFile.getName().endsWith(".csv")){
			String[] header=br.readLine().substring(1).split(",");
			for(String s:header){
				methodsDB.add(methodsDBTmp);
				massbankDB.add(massbank.get(getIDDatasetChargeOrigin(s)));
			}
			String line;
			while((line=br.readLine())!=null){
				String tmp=line.substring(0,line.indexOf(","));
				methodsQuery.add(methodsQueryTmp);
				massbankQuery.add(massbank.get(getIDDatasetChargeOrigin(tmp)));
			}
			br.close();
			similarity=new double[massbankQuery.size()][massbankDB.size()];
			commonPeaks=new int[massbankQuery.size()][massbankDB.size()];
			br=new BufferedReader(new FileReader(similarityFile));
			br.readLine();br.readLine();br.readLine();br.readLine();
			int r=0;
			while((line=br.readLine())!=null){
				String[] tmp=line.substring(line.indexOf(",")+1).split(",");
				int c=0;
				for(String s:tmp){
					if(s.contains(" ")){
						similarity[r][c]=Double.parseDouble(s.split(" ")[0]);
						commonPeaks[r][c]=Integer.parseInt(s.split(" ")[1]);
					}else{
						similarity[r][c]=Double.parseDouble(s);
						commonPeaks[r][c]=0;
					}
					c++;
				}
				r++;
			}
		}else{
			String line;
			while((line=br.readLine())!=null){
				if(line.startsWith("query")){
					MassBank mb=massbank.get(getIDDatasetChargeOrigin(line.split("\\t")[1]));
					if(!massbankQuery.contains(mb)){
						massbankQuery.add(mb);
						methodsQuery.add(methodsQueryTmp);
					}
				}
				if(line.startsWith("result")){
					MassBank mb=massbank.get(getIDDatasetChargeOrigin(line.split("\\t")[1]));
					if(!massbankDB.contains(mb)){
						massbankDB.add(mb);
						methodsDB.add(methodsDBTmp);
					}
				}
			}
			String querySub="";
			if(massbankQuery.size()==0){
				querySub=methodsQueryTmp.startsWith("MetlinPositive")?methodsQueryTmp.replaceAll("MetlinPositiv", "MP"):methodsQueryTmp.replaceAll("AgilentPositiv", "AP");
			}else{
				querySub=massbankQuery.get(0).massbankID.replaceAll("\\d", "");
			}
			querySub=getIDDatasetChargeOrigin(querySub);

			for(MassBank mb:massbank.values())
				if(mb.massbankID.startsWith(querySub)&&!massbankQuery.contains(mb)){
					massbankQuery.add(mb);
					methodsQuery.add(methodsQueryTmp);
				}	

			String dbSub="";
			if(massbankDB.size()==0){
				dbSub=methodsDBTmp.startsWith("MetlinPositive")?methodsDBTmp.replaceAll("MetlinPositive", "MP"):methodsDBTmp.replaceAll("AgilentPositive", "AP");
			}else{
				dbSub=massbankDB.get(0).massbankID.replaceAll("\\d", "");
			}
			dbSub=getIDDatasetChargeOrigin(dbSub);

			for(MassBank mb:massbank.values())
				if(mb.massbankID.startsWith(dbSub)&&!massbankDB.contains(mb)){
					massbankDB.add(mb);
					methodsDB.add(methodsDBTmp);
				}

			br.close();
			similarity=new double[massbankQuery.size()][massbankDB.size()];
			commonPeaks=new int[massbankQuery.size()][massbankDB.size()];
			br=new BufferedReader(new FileReader(similarityFile));
			int indexQuery=Integer.MAX_VALUE;
			while((line=br.readLine())!=null){
				if(line.startsWith("query")){
					MassBank mb=massbank.get(getIDDatasetChargeOrigin(line.split("\\t")[1]));
					indexQuery=massbankQuery.indexOf(mb);
				}
				if(line.startsWith("result")){
					String l[]=line.split("\\t");
					MassBank mb=massbank.get(getIDDatasetChargeOrigin(l[1]));
					int indexDB=massbankDB.indexOf(mb);
					similarity[indexQuery][indexDB]=Double.parseDouble(l[4]);
					commonPeaks[indexQuery][indexDB]=Integer.parseInt(l[3]);
				}
			}
		}
		br.close();


		clearExclusionList();
	}

	public SimilarityMatrix(){}

	public SimilarityMatrix(List<SimilarityMatrix> matrices) throws Exception{		
		for(SimilarityMatrix m:matrices){
			similarityFiles.addAll(m.similarityFiles);
			typesQuery.putAll(m.typesQuery);
			typesDB.putAll(m.typesDB);			
			massbankDB.addAll(m.massbankDB);
			methodsDB.addAll(m.methodsDB);
			if(massbankQuery.isEmpty()){
				massbankQuery.addAll(m.massbankQuery);
				methodsQuery.addAll(m.methodsQuery);
			}else{
				if(!massbankQuery.containsAll(m.massbankQuery)||
						!m.massbankQuery.containsAll(massbankQuery)||
						!methodsQuery.containsAll(m.methodsQuery)||
						!m.methodsQuery.containsAll(methodsQuery)){
					throw new Exception("No merge possible as the query was not the same!");
				}
			}
		}

		similarity=new double[massbankQuery.size()][massbankDB.size()];
		commonPeaks=new int[massbankQuery.size()][massbankDB.size()];

		clearExclusionList();

		for(MassBank mb:massbankQuery){
			int c=0;
			for(SimilarityMatrix m:matrices){
				int r=m.massbankQuery.indexOf(mb);
				for(int i=0;i<m.similarity[r].length;i++){
					similarity[r][c]=m.similarity[r][i];
					commonPeaks[r][c]=m.commonPeaks[r][i];
					exclusionListDB[r][c]=m.exclusionListDB[r][i];
					c++;
				}

			}
		}

	}

	public void clearExclusionList(){
		exclusionListDB=new boolean[similarity.length][similarity[0].length];
		for(int i=0;i<exclusionListDB.length;i++){
			for(int j=0;j<exclusionListDB[i].length;j++){
				exclusionListDB[i][j]=Boolean.FALSE;
			}
		}
	}

	public static void clearExclusionLists(List<SimilarityMatrix> matrices){
		for(SimilarityMatrix m:matrices){
			m.clearExclusionList();			
		}
	}

	public void excludeSameInChi(){
		for(int i=0;i<massbankQuery.size();i++){
			for(int j=0;j<massbankDB.size();j++){
				if(massbankQuery.get(i).hasEqualInChiKey(massbankDB.get(j))){
					exclusionListDB[i][j]=Boolean.TRUE;
				}
			}
		}
	}

	public static void excludeSameInChi(List<SimilarityMatrix> matrices){
		for(SimilarityMatrix m:matrices){
			m.excludeSameInChi();			
		}
	}

	public void excludeZeroEntries(){
		for(int i=0;i<massbankQuery.size();i++){
			for(int j=0;j<massbankDB.size();j++){
				if(similarity[i][j]==0){
					exclusionListDB[i][j]=Boolean.TRUE;
				}
			}
		}
	}

	public void excludeDifferentMass(double ppm, double ae){
		for(int i=0;i<massbankQuery.size();i++){
			for(int j=0;j<massbankDB.size();j++){		
				if(!massbankQuery.get(i).hasEqualMass(massbankDB.get(j), ppm, ae)){
					exclusionListDB[i][j]=Boolean.TRUE;
				}
			}
		}
	}
	
	public void excludeMassBanksFromDB(Map<MassBank, List<MassBank>> mb){
		for(int i=0;i<massbankQuery.size();i++){
			List<MassBank> currMB=mb.get(massbankQuery.get(i));
			if(currMB==null)continue;
			for(int j=0;j<massbankDB.size();j++){		
				if(currMB.contains(massbankDB.get(j))){
					exclusionListDB[i][j]=Boolean.TRUE;
				}
			}
		}
	}

	public void excludeDifferentFormulas(){
		for(int i=0;i<massbankQuery.size();i++){
			for(int j=0;j<massbankDB.size();j++){
				if(!massbankQuery.get(i).mf.equals(massbankDB.get(j).mf)){
					exclusionListDB[i][j]=Boolean.TRUE;
				}
			}
		}
	}

	public static void excludeDifferentMasses(List<SimilarityMatrix> matrices, double ppm, double ae){
		for(SimilarityMatrix m:matrices){
			m.excludeDifferentMass(ppm, ae);			
		}
	}

	public static double getAbsoluteErrorForMass(double precursor, double ppm, double ae) {
		return Math.max(ppm*1e-6*precursor, ae*1e-3);
	}

	public static void excludeZeroEntries(List<SimilarityMatrix> matrices){
		for(SimilarityMatrix m:matrices){
			m.excludeZeroEntries();
		}
	}

	public boolean isExcluded(MassBank query, MassBank db){
		return exclusionListDB[massbankQuery.indexOf(query)][massbankDB.indexOf(db)];
	}

	public double getSimilarityValue(int i,int j){
		if(exclusionListDB[i][j])return Double.NaN;
		return similarity[i][j];
	}

	static public String getIDDatasetChargeOrigin(String methodID){
		return methodID.substring(0,7)+methodID.replaceAll("\\D", "");
	}

	static public String getIDDatasetCharge(String methodID){
		return methodID.substring(0,2)+methodID.replaceAll("\\D", "");
	}

	public int[][] getIndizesOfSortedScores(){
		int[][] result=new int[methodsQuery.size()][methodsDB.size()];
		for(int i=0;i<result.length;i++){
			for(int k=0;k<result[i].length;k++){
				result[i][k]=k;
			}
		}
		for(int i=0;i<result.length;i++){
			double[] similarityTmp=new double[methodsDB.size()];
			for(int k=0;k<similarityTmp.length;k++){
				similarityTmp[k]=getSimilarityValue(i,k);
			}
			Utils.quicksort(similarityTmp,result[i]);
		}
		return result;
	}

	public XYSeries getEstimatedQValueVSCalculatedQValue(boolean FDR, List<MassBank> compoundsOfInterest, String number, double minCalcValue, double maxCalcValue) throws IOException{		
		int[][] sortedIndizes=getIndizesOfSortedScores();
		List<Double> bestScores=new ArrayList<Double>();
		List<MatchType> matches=new ArrayList<MatchType>();
		for(int i=0;i<sortedIndizes.length;i++){
			if(!compoundsOfInterest.contains(massbankQuery.get(i)))continue;
			boolean foundOriginal=false;
			boolean foundDecoy=false;			
			for(int j=0;j<sortedIndizes[i].length;j++){
				int currCol=sortedIndizes[i][j];
				if(Double.isNaN(getSimilarityValue(i,currCol)))continue;
				String method=methodsDB.get(currCol);
				if(!foundOriginal&&typesDB.get(method)==DatabaseType.Original){
					bestScores.add(getSimilarityValue(i,currCol));
					if(massbankQuery.get(i).hasEqualInChiKey(massbankDB.get(currCol))){
						matches.add(MatchType.TruePositiveMatch);
					}else{
						matches.add(MatchType.FalsePositiveMatch);
					}
					foundOriginal=true;
				}
				if(!foundDecoy&&typesDB.get(method)==DatabaseType.Decoy){
					bestScores.add(getSimilarityValue(i,currCol));					
					matches.add(MatchType.DecoyMatch);
					foundDecoy=true;
				}
				if(foundOriginal&&foundDecoy)break;
			}
		}
		int[] indizes=new int[bestScores.size()];
		double[] scores=new double[bestScores.size()]; 
		for(int i=0;i<indizes.length;i++)indizes[i]=i;
		for(int i=0;i<scores.length;i++)scores[i]=bestScores.get(i);
		Utils.quicksort(scores,indizes);

		List<Double> FDRByDecoy=new ArrayList<Double>();
		List<Double> FDRCalculated=new ArrayList<Double>();
		int countTrueMatches=0;
		int countFalseMatches=0;
		int countDecoyMatches=0;
		for(int i:indizes){
			if(matches.get(i)==MatchType.DecoyMatch){
				countDecoyMatches++;
			}else{
				if(matches.get(i)==MatchType.TruePositiveMatch){
					countTrueMatches++;
				}else{
					countFalseMatches++;
				}
				//				FDRByDecoy.add(2.0*countDecoyMatches/(countDecoyMatches+countTrueMatches+countFalseMatches));
				FDRByDecoy.add(1.0*countDecoyMatches/(countTrueMatches+countFalseMatches));
				FDRCalculated.add(1.0*countFalseMatches/(countFalseMatches+countTrueMatches));
			}
		}
		FDRByDecoy.add(2.0*countDecoyMatches/(countDecoyMatches+countTrueMatches+countFalseMatches));
		FDRCalculated.add(1.0*countFalseMatches/(countFalseMatches+countTrueMatches));

		if(!FDR){
			double min=Double.POSITIVE_INFINITY;
			for(int i=FDRByDecoy.size()-1;i>=0;i--){
				min=Math.min(FDRByDecoy.get(i),min);
				FDRByDecoy.set(i, min);
			}

			min=Double.POSITIVE_INFINITY;
			for(int i=FDRCalculated.size()-1;i>=0;i--){
				min=Math.min(FDRCalculated.get(i),min);
				FDRCalculated.set(i, min);
			}
		}

		for(int i=FDRCalculated.size()-1;i>=0;i--){
			if(FDRCalculated.get(i)<minCalcValue||FDRCalculated.get(i)>maxCalcValue){
				FDRCalculated.remove(i);
				FDRByDecoy.remove(i);
			}
		}

		XYSeries series=new XYSeries((FDR?"FDR":"q-Values")+" " +number+" "+ getMethodsQueryAndDBString());
		for(int i=0;i<FDRByDecoy.size();i++){
			series.add(FDRCalculated.get(i),FDRByDecoy.get(i));
		}

		return series;
	}

	public static XYSeries getBisectingLine(double min, double max){
		XYSeries seriesWH=new XYSeries("bisecting line");
		for(double i=min;i<=max;i+=(max-min)/10)seriesWH.add(i,i);
		return seriesWH;
	}

	public static void writeEstimatedQValueVSCalculatedQValue(File outputFile, List<SimilarityMatrix> matrices, String add, double percentage, int repeat, double minCalcValue, double maxCalcValue) throws IOException{
		writeEstimatedQValueVSCalculatedQValue(outputFile, matrices, add, false, percentage, repeat, minCalcValue, maxCalcValue);
	}
	public static void writeEstimatedFDRVSCalculatedFDR(File outputFile, List<SimilarityMatrix> matrices, String add, double percentage, int repeat, double minCalcValue, double maxCalcValue) throws IOException{
		writeEstimatedQValueVSCalculatedQValue(outputFile, matrices, add, true, percentage, repeat, minCalcValue, maxCalcValue);
	}

	public static void writeEstimatedQValueVSCalculatedQValue(File outputFile, List<SimilarityMatrix> matrices, String add, boolean FDR, double percentage, int repeat, double minCalcValue, double maxCalcValue) throws IOException{
		Random r=new Random();	
		List<List<MassBank>> compoundsOfInterest=new ArrayList<List<MassBank>>();
		for(int i=0;i<repeat;i++){
			List<MassBank> allMassBanks=new ArrayList<MassBank>(matrices.get(0).massbankQuery);
			int size=allMassBanks.size();
			List<MassBank> finalMassBanks=new ArrayList<MassBank>();
			while(finalMassBanks.size()<percentage*size){
				int n=r.nextInt(allMassBanks.size());
				finalMassBanks.add(allMassBanks.get(n));
				allMassBanks.remove(n);
			}
			compoundsOfInterest.add(allMassBanks);
		}

		for(SimilarityMatrix m:matrices){
			m.writeEstimatedQValueVSCalculatedQValue(outputFile, add,  FDR, compoundsOfInterest, minCalcValue, maxCalcValue);
		}
	}

	public void writeEstimatedQValueVSCalculatedQValue(File outputFile, String add, boolean FDR, List<List<MassBank>> compoundsOfInterest, double minCalcValue, double maxCalcValue) throws IOException{
		XYSeriesCollection dataset = new XYSeriesCollection();
		List<XYSeries> series=new ArrayList<XYSeries>();
		for(int i=0;i<compoundsOfInterest.size();i++){
			//			System.out.println(i);
			List<MassBank> finalMassBanks=compoundsOfInterest.get(i);
			//			series.add(this.getEstimatedQValueVSCalculatedQValue(FDR, finalMassBanks, Integer.toString(i), minCalcValue, maxCalcValue));
			series.add(this.getEstimatedQValueVSCalculatedQValue(FDR, finalMassBanks, "", minCalcValue ,maxCalcValue));
		}

		String desc=FDR?"FDRsComparisonMultiple_":"qValuesComparisonMultiple_";
		String desc2=FDR?"FDRs":"qValues";

		BufferedWriter bw=new BufferedWriter(new FileWriter(new File(outputFile.getAbsolutePath()+sep+add+desc+getMethodsQueryAndDBString()+".txt")));

		XYSeries xyseries=new XYSeries(series.get(0).getKey());

		double sum=0;
		int count=0;
		for(XYSeries s:series){
			double sum2=0;
			int count2=0;
			for(Object iTmp:s.getItems()){
				XYDataItem i=(XYDataItem) iTmp;
				xyseries.add(i.getXValue(), i.getYValue());
				double curr=Math.pow(i.getXValue()-i.getYValue(),2);
				sum2+=curr;
				count2++;
			}
			bw.write(s.getKey()+"\t"+count2+"\t"+sum2+"\t"+1.0*sum2/count2);bw.newLine();
			if(count2!=0){
				sum+=1.0*sum2/count2;
				count++;
			}
		}
		bw.write("average sum\t"+count+"\t"+sum+"\t"+1.0*sum/count);bw.newLine();
		bw.close();

		dataset.addSeries(xyseries);

		XYSeries seriesWH=getBisectingLine(minCalcValue, maxCalcValue);
		dataset.addSeries(seriesWH);

		final JFreeChart chart = ChartFactory.createScatterPlot(desc2, desc2 +" calculated", desc2+" by decoy database", dataset);
		ChartUtilities.saveChartAsJPEG(new File(outputFile.getAbsolutePath()+sep+add+desc+getMethodsQueryAndDBString()+".jpg"), chart, 1000, 1000);
	}

	public static void writeEstimatedQValueVSCalculatedQValue(File outputFile, List<SimilarityMatrix> matrices, String add, boolean FDR) throws IOException{
		XYSeriesCollection dataset = new XYSeriesCollection();
		for(SimilarityMatrix m:matrices){
			dataset.addSeries(m.getEstimatedQValueVSCalculatedQValue(FDR, m.massbankQuery,"",Double.NEGATIVE_INFINITY,Double.POSITIVE_INFINITY));				
		}

		XYSeries seriesWH=getBisectingLine(0, 1);
		dataset.addSeries(seriesWH);

		String desc=FDR?"FDRsComparison_":"qValuesComparison_";
		String desc2=FDR?"FDRs":"qValues";

		final JFreeChart chart = ChartFactory.createScatterPlot(desc2, desc2+" calculated", desc2+" by decoy database", dataset);
		ChartUtilities.saveChartAsJPEG(new File(outputFile.getAbsolutePath()+sep+add+desc+getMethodsQueryStringFromList(matrices)+".jpg"), chart, 1000, 1000);
	}

	public void writeEstimatedQValueVSCalculatedQValue(File outputFile, String add) throws IOException{
		SimilarityMatrix.writeEstimatedQValueVSCalculatedQValue(outputFile, Arrays.asList(new SimilarityMatrix[]{this}), add, false);
	}

	public void writeEstimatedFDRVSCalculatedFDR(File outputFile, String add) throws IOException{
		SimilarityMatrix.writeEstimatedQValueVSCalculatedQValue(outputFile, Arrays.asList(new SimilarityMatrix[]{this}), add, true);	
	}

	public static void writeEstimatedQValueVSCalculatedQValue(File outputFile, List<SimilarityMatrix> matrices, String add) throws IOException{
		SimilarityMatrix.writeEstimatedQValueVSCalculatedQValue(outputFile, matrices, add, false);
	}

	public static void writeEstimatedFDRVSCalculatedFDR(File outputFile, List<SimilarityMatrix> matrices, String add) throws IOException{
		SimilarityMatrix.writeEstimatedQValueVSCalculatedQValue(outputFile, matrices, add, true);
	}

	public String getMethodsQueryAndDBString(){
		return getMethodsQueryString()+"---"+getMethodsDBString();
	}

	public String getMethodsQueryString(){
		String result="";
		for(String s:typesQuery.keySet())result+="_"+s;
		return result.substring(1).replaceAll("[a-z]","");
	}

	public String getMethodsDBString(){
		String result="";
		for(String s:typesDB.keySet())result+="_"+s;
		return result.substring(1).replaceAll("[a-z]","");
	}

	private String getMethodFromFilename(String s){
		String[] tmp=s.split(s.contains("/")?"/":"\\\\");
		return tmp[tmp.length-2];
	}


	//	public void writeScoresOfSameCompounds(File outputFile, String add) throws Exception{
	//		HashMap<String, Map<MassBank,Integer>> massbank2indexDB=new HashMap<String, Map<MassBank,Integer>>();
	//		List<String> massbankMethods=new ArrayList<String>();
	//		List<MassBank> massbankEntries=new ArrayList<MassBank>();
	//		for(String t:methodsDB)if(!massbank2indexDB.containsKey(t))massbank2indexDB.put(t, new HashMap<MassBank,Integer>());
	//		for(int i=0;i<massbankDB.size();i++){
	//			MassBank mb=massbankDB.get(i);
	//			String method=methodsDB.get(i);
	//			massbank2indexDB.get(method).put(mb,i);
	//			if(!massbankMethods.contains(method))massbankMethods.add(method);
	//			if(!massbankEntries.contains(mb))massbankEntries.add(mb);
	//		}
	//
	//		XYSeriesCollection dataset = new XYSeriesCollection();
	//		Map<String, XYSeries> seriesMap=new HashMap<String, XYSeries>();
	//		int[][] sortedIndizes=getIndizesOfSortedScores();
	//		for(int i=0;i<sortedIndizes.length;i++){
	//			for(int m=0;m<1;m++){
	//				int ind=sortedIndizes[i][m];
	//				MassBank mb=massbankDB.get(ind);
	//				double valueOrig=getSimilarityValue(i,ind);
	//				if(Double.isNaN(valueOrig))continue;
	//				for(int j=0;j<massbankMethods.size()-1;j++){
	//					for(int k=j+1;k<massbankMethods.size();k++){				
	//						String methodsA=(methodsDB.get(j).compareTo(massbankMethods.get(k))<0)?massbankMethods.get(j):massbankMethods.get(k);
	//						String methodsB=(methodsDB.get(j).compareTo(massbankMethods.get(k))<0)?massbankMethods.get(k):massbankMethods.get(j);
	//						String key=methodsQuery.get(i)+" "+methodsA+" "+methodsB;				
	//						int indexA=massbank2indexDB.get(methodsA).get(mb);
	//						int indexB=massbank2indexDB.get(methodsB).get(mb);
	//						if(!seriesMap.containsKey(key)){
	//							String description="DB: "+methodsQuery.get(i)+", "+methodsA+" vs. "+methodsB;
	//							XYSeries series=new XYSeries(description);
	//							seriesMap.put(key, series);
	//							dataset.addSeries(series);
	//						}
	//						seriesMap.get(key).add(getSimilarityValue(i,indexA),getSimilarityValue(i,indexB));
	//					}
	//				}
	//			}
	//		}
	//		XYSeries seriesWH=getBisectingLine(0,1);
	//		dataset.addSeries(seriesWH);
	//
	//		final JFreeChart chart = ChartFactory.createScatterPlot("Score 1 vs. Score 2", "Score 1", "Score2", dataset);
	//		ChartUtilities.saveChartAsJPEG(new File(outputFile.getAbsolutePath()+sep+add+"ScoresOfSameCompounds_"+getMethodsQueryAndDBString()+".jpg"), chart, 1000, 1000);
	//	}

	//	public void writeScoreDistributionOfTopRank(File outputFile, double step, String add) throws Exception{
	//		Map<Double, Integer> topScoringMatchesAll=new TreeMap<Double,Integer>();
	//		Map<Double, Integer> topScoringMatchesTrue=new TreeMap<Double,Integer>();
	//		Map<Double, Integer> topScoringMatchesFalse=new TreeMap<Double,Integer>();
	//		int numberTrueMatches=0;
	//		int numberFalseMatches=0;
	//
	//		int[][] sortedIndizes=getIndizesOfSortedScores();
	//
	//		for(int i=0;i<sortedIndizes.length;i++){
	//			int ind=sortedIndizes[i][0];
	//			double valueOrig=getSimilarityValue(i,ind);
	//			if(Double.isNaN(valueOrig))continue;			
	//			double value=Math.round(valueOrig/step)*step;
	//			if(!topScoringMatchesAll.containsKey(value))topScoringMatchesAll.put(value, 0);
	//			topScoringMatchesAll.put(value,topScoringMatchesAll.get(value)+1);
	//			if(massbankQuery.get(i).hasEqualInChiKey(massbankDB.get(ind))){
	//				if(!topScoringMatchesTrue.containsKey(value))topScoringMatchesTrue.put(value, 0);
	//				topScoringMatchesTrue.put(value,topScoringMatchesTrue.get(value)+1);
	//				numberTrueMatches++;
	//			}else{
	//				if(!topScoringMatchesFalse.containsKey(value))topScoringMatchesFalse.put(value, 0);
	//				topScoringMatchesFalse.put(value,topScoringMatchesFalse.get(value)+1);
	//				numberFalseMatches++;
	//			}
	//		}
	//
	//		XYSeries seriesTrue = new XYSeries("True Hits ("+numberTrueMatches+")");
	//		for(Entry<Double,Integer> m:topScoringMatchesTrue.entrySet()){
	//			seriesTrue.add(m.getKey(),new Double(1.0*m.getValue()/numberTrueMatches));
	//		}
	//		XYSeries seriesFalse = new XYSeries("False Hits ("+numberFalseMatches+")");
	//		for(Entry<Double,Integer> m:topScoringMatchesFalse.entrySet()){
	//			seriesFalse.add(m.getKey(),new Double(1.0*m.getValue()/numberFalseMatches));
	//		}
	//		XYSeries seriesAll= new XYSeries("All Hits ("+(numberTrueMatches+numberFalseMatches)+")");
	//		for(Entry<Double,Integer> m:topScoringMatchesAll.entrySet()){
	//			seriesAll.add(m.getKey(),new Double(1.0*m.getValue()/(numberTrueMatches+numberFalseMatches)));
	//		}
	//
	//		XYSeriesCollection dataset = new XYSeriesCollection();
	//		dataset.addSeries(seriesTrue);
	//		dataset.addSeries(seriesFalse);
	//		dataset.addSeries(seriesAll);
	//
	//		final JFreeChart chart =ChartFactory.createXYLineChart("Score Distribution", "score", "percentage", dataset);
	//		ChartUtilities.saveChartAsJPEG(new File(outputFile.getPath()+sep+add+"ScoreDistributionTopRank"+"_"+getMethodsQueryAndDBString()+".jpg"), chart, 1000, 400);
	//	}

	public double[] writeROCCurveOfTopRank(File outputFile, double step, String add) throws Exception{
		double AUC=0;
		double F=0;
		int countF=0;
		TreeMap<Double, List<MatchType>> topScoringMatches=new TreeMap<Double,List<MatchType>>(Collections.reverseOrder());
		int[][] sortedIndizes=getIndizesOfSortedScores();
		int allPositives=0;
		int allNegatives=0;

		for(int i=0;i<sortedIndizes.length;i++){
			int ind=sortedIndizes[i][0];
			double value=getSimilarityValue(i,ind);
			if(Double.isNaN(value))continue;
			if(!topScoringMatches.containsKey(value))topScoringMatches.put(value, new ArrayList<MatchType>());			
			if(massbankQuery.get(i).hasEqualInChiKey(massbankDB.get(ind))){
				topScoringMatches.get(value).add(MatchType.TruePositiveMatch);
				allPositives++;
			}else{
				topScoringMatches.get(value).add(MatchType.FalsePositiveMatch);
				allNegatives++;
			}
		}

		List<double[]> values=new ArrayList<double[]>();
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries series=new XYSeries("ROCCurve " + getMethodsQueryAndDBString());
		double currentNumberTruePositives=0;
		double currentNumberFalsePositives=0;
		for(Entry<Double, List<MatchType>> matches:topScoringMatches.entrySet()){
			for(MatchType m:matches.getValue()){
				if(m.equals(MatchType.TruePositiveMatch))currentNumberTruePositives++;
				if(m.equals(MatchType.FalsePositiveMatch))currentNumberFalsePositives++;
			}

			double sensitivity=1.0*currentNumberTruePositives/allPositives;
			double EinsMinusSpecificity=1.0*currentNumberFalsePositives/allNegatives;
			values.add(new double[]{EinsMinusSpecificity,sensitivity});
			series.add(EinsMinusSpecificity,sensitivity);
			double precision=1.0*currentNumberTruePositives/(currentNumberTruePositives+currentNumberFalsePositives);
			F+=matches.getValue().size()*2*precision*sensitivity/(precision+sensitivity);
			countF+=matches.getValue().size();
		}
		dataset.addSeries(series);
		F/=countF;

		double lastX=0;		
		for(int i=1;i<values.size();i++){
			if(values.get(i-1)[0]!=values.get(i)[0]){
				AUC+=values.get(i-1)[1]*(values.get(i)[0]-lastX);
				lastX=values.get(i)[0];
			}
		}

		final JFreeChart chart = ChartFactory.createScatterPlot("ROCCurve","False Positive Rate", "True Positive Rate", dataset);
		ChartUtilities.saveChartAsJPEG(new File(outputFile.getAbsolutePath()+sep+add+"ROCCurveOfTopRank_"+getMethodsQueryAndDBString()+".jpg"), chart, 1000, 1000);

		return new double[]{AUC,F};
	}

	public void writeRankVSPercentage(File outputFile, int length, String add) throws Exception{

		int[][] sortedIndizes=getIndizesOfSortedScores();
		Map<String,int[]> numbers=new TreeMap<String,int[]>();
		numbers.put("Decoy", new int[Math.min(length,sortedIndizes[0].length)]);
		numbers.put("Original", new int[Math.min(length,sortedIndizes[0].length)]);

		for(int i=0;i<sortedIndizes.length;i++){
			int l=0;
			for(int j=0;j<sortedIndizes[i].length;j++){
				int currCol=sortedIndizes[i][j];
				if(Double.isNaN(getSimilarityValue(i,currCol)))continue;
				String method=methodsDB.get(currCol);
				if(typesDB.get(method)==DatabaseType.Original){
					numbers.get("Original")[l]=numbers.get("Original")[l]+1;
				}else if(typesDB.get(method)==DatabaseType.Decoy){
					numbers.get("Decoy")[l]=numbers.get("Decoy")[l]+1;					
				}
				l++;
				if(l>=numbers.get("Original").length)break;
			}			
		}
		final DefaultCategoryDataset dataset = new DefaultCategoryDataset();
		for(Entry<String, int[]> n:numbers.entrySet()){
			for(int i=0;i<n.getValue().length;i++){
				dataset.addValue(n.getValue()[i], n.getKey(), Integer.toString(i));
			}
		}

		final JFreeChart chart = ChartFactory.createBarChart("BarChart","Ranks","Percentage",dataset);
		ChartUtilities.saveChartAsJPEG(new File(outputFile.getPath()+sep+add+"RankVSpercentage_"+getMethodsQueryAndDBString()+".jpg"), chart, Math.min(2000,length*100), 1000);			
	}

	//	public void writeNumberHitsVSPercentage(File outputFile, String add) throws Exception{
	//
	//		TreeMap<Integer,Integer> numbers=new TreeMap<Integer,Integer>();
	//
	//		for(int i=0;i<massbankQuery.size();i++){
	//			int countHits=0;
	//			for(int j=0;j<massbankDB.size();j++){
	//				if(Double.isNaN(getSimilarityValue(i,j)))continue;
	//				countHits++;
	//			}
	//			if(!numbers.containsKey(countHits))numbers.put(countHits, 0);
	//			numbers.put(countHits, numbers.get(countHits)+1);
	//		}
	//
	//		final DefaultCategoryDataset dataset = new DefaultCategoryDataset();
	//		for(Entry<Integer, Integer> n:numbers.entrySet()){			
	//			dataset.addValue(n.getValue(), "numbers", Integer.toString(n.getKey()));
	//		}
	//
	//		final JFreeChart chart = ChartFactory.createBarChart("BarChart","Number of Hits","Percentage",dataset);
	//		ChartUtilities.saveChartAsJPEG(new File(outputFile.getPath()+sep+add+"NumberHitsVSPercentage_"+getMethodsQueryAndDBString()+".jpg"), chart, 2000, 1000);			
	//	}

	public static void writeRankVSScore(File outputFile, List<SimilarityMatrix> matrices, int length, String add) throws Exception{		
		final DefaultBoxAndWhiskerCategoryDataset dataset = new DefaultBoxAndWhiskerCategoryDataset();
		for(SimilarityMatrix m:matrices){
			int[][] sortedIndizes=m.getIndizesOfSortedScores();
			List<Double>[] scoreList=new List[Math.min(length,sortedIndizes[0].length)];
			for(int i=0;i<scoreList.length;i++)scoreList[i]=new ArrayList<Double>();
			for(int i=0;i<sortedIndizes.length;i++){
				int l=0;
				for(int j=0;j<sortedIndizes[i].length;j++){
					double v=m.getSimilarityValue(i,sortedIndizes[i][j]);
					if(!Double.isNaN(v)){
						scoreList[l].add(v);
						l++;						
					}
					if(l>=scoreList.length)break;
				}			
			}
			for(int i=0;i<scoreList.length;i++)
				dataset.add(scoreList[i], m.getMethodsQueryAndDBString(), i);

		}

		JFreeChart chart=ChartFactory.createBoxAndWhiskerChart("Whiskerplot", "Ranks", "Scores", dataset, true);
		ChartUtilities.saveChartAsJPEG(new File(outputFile.getPath()+sep+add+"RankVSScore_"+getMethodsQueryStringFromList(matrices)+".jpg"), chart, Math.min(2000,length*100), 1000);
	}

	public static String getMethodsQueryAndDBStringFromList(List<SimilarityMatrix> matrices){
		String name="";
		for(SimilarityMatrix m:matrices){
			name+="___"+m.getMethodsQueryAndDBString();
		}
		return name.substring(3);
	}

	public static String getMethodsDBStringFromList(List<SimilarityMatrix> matrices){
		String name="";
		for(SimilarityMatrix m:matrices){
			name+="___"+m.getMethodsDBString();
		}
		return name.substring(3);
	}

	public static String getMethodsQueryStringFromList(List<SimilarityMatrix> matrices){
		String name="";
		List<String> methods=new ArrayList<String>();
		for(SimilarityMatrix m:matrices){
			String me = m.getMethodsQueryString();
			if(!methods.contains(me)){
				methods.add(me);
				name+="___"+me;
			}			
		}
		return name.substring(3);
	}

	public int getNumberOfPossibleHits(){
		int number=0;
		for(int i=0;i<similarity.length;i++){
			boolean found=false;
			for(int j=0;j<similarity[i].length;j++){
				double currValue=getSimilarityValue(i,j);
				if(Double.isNaN(currValue))continue;
				if(massbankQuery.get(i).hasEqualInChiKey(massbankDB.get(j))){		
					found=true;
					break;
				}
			}
			if(found)number++;
		}
		return number;
	}

	public int getNumberOfNonemptySearchLists(){
		int number=0;
		for(int i=0;i<similarity.length;i++){
			boolean found=false;
			for(int j=0;j<similarity[i].length;j++){
				double currValue=getSimilarityValue(i,j);
				if(Double.isNaN(currValue))continue;
				found=true;
				break;
			}
			if(found)number++;
		}
		return number;
	}


	public int getNumberOfTrueHits(){
		int number=0;
		int[][] sortedIndizes=getIndizesOfSortedScores();
		for(int i=0;i<sortedIndizes.length;i++){
			for(int j=0;j<sortedIndizes[i].length;j++){
				int currCol=sortedIndizes[i][j];
				double currValue=getSimilarityValue(i,currCol);
				if(Double.isNaN(currValue))continue;
				if(massbankQuery.get(i).hasEqualInChiKey(massbankDB.get(currCol))){
					number++;
				}
				break;
			}
		}
		return number;
	}

	public List<MassBank> getFirstSortedMassbankEntries(int i, int[] sortedIndizes){
		List<MassBank> result=new ArrayList<MassBank>();

		double maxValue=Double.NaN;
		for(int j=0;j<sortedIndizes.length;j++){
			int currCol=sortedIndizes[j];
			double currValue=getSimilarityValue(i,currCol);
			if(Double.isNaN(currValue))continue;
			if(Double.isNaN(maxValue)||maxValue==currValue){
				result.add(massbankDB.get(currCol));
				maxValue=currValue;
			}else{
				break;
			}
		}
		return result;
	}

	public List<List<MassBank>> getSortedMassbankEntries(int i, int[] sortedIndizes){
		List<List<MassBank>> result=new ArrayList<List<MassBank>>();
		List<MassBank> currResult=new ArrayList<MassBank>();

		double currMaxValue=Double.NaN;
		for(int j=0;j<sortedIndizes.length;j++){
			int currCol=sortedIndizes[j];
			double currValue=getSimilarityValue(i,currCol);
			if(Double.isNaN(currValue))continue;
			if(Double.isNaN(currMaxValue)||Double.compare(currMaxValue,currValue)!=0){
				currResult=new ArrayList<MassBank>();
				result.add(currResult);
			}
			currResult.add(massbankDB.get(currCol));
			currMaxValue=currValue;		
		}
		return result;
	}

	public List<List<DatabaseType>> getSortedDatabaseTypes(int i, int[] sortedIndizes){
		List<List<DatabaseType>> result=new ArrayList<List<DatabaseType>>();
		List<DatabaseType> currResult=new ArrayList<DatabaseType>();

		double currMaxValue=Double.NaN;
		for(int j=0;j<sortedIndizes.length;j++){
			int currCol=sortedIndizes[j];
			double currValue=getSimilarityValue(i,currCol);
			if(Double.isNaN(currValue))continue;
			if(Double.isNaN(currMaxValue)||Double.compare(currMaxValue,currValue)!=0){
				currResult=new ArrayList<DatabaseType>();
				result.add(currResult);
			}
			currResult.add(typesDB.get(methodsDB.get(currCol)));
			currMaxValue=currValue;		
		}
		return result;
	}

	public List<MassBank> getHitsInDB(MassBank mbQuery){
		List<MassBank> result=new ArrayList<MassBank>();
		for(MassBank mbDB:massbankDB){
			if(mbQuery.hasEqualInChiKey(mbDB)&&!isExcluded(mbQuery, mbDB))
				result.add(mbDB);
		}
		return result;
	}

	public static Map<DatabaseType, Double>[] getRankDistribution(File outputFile, SimilarityMatrix mOriginal, List<SimilarityMatrix> mDecoys, int length) throws Exception{
		Map<DatabaseType, Double>[] result=new TreeMap[length];
		for(int i=0;i<result.length;i++)result[i]=new TreeMap<DatabaseType, Double>();
		for(SimilarityMatrix matrix:mDecoys){
			SimilarityMatrix merged=new SimilarityMatrix(Arrays.asList(new SimilarityMatrix[]{mOriginal, matrix}));					
			int[][] sortedIndizes=merged.getIndizesOfSortedScores();
			for(int i=0;i<merged.massbankQuery.size();i++){
				List<List<DatabaseType>> sortedDatabaseType = merged.getSortedDatabaseTypes(i, sortedIndizes[i]);
				int j=0;
				int index=0;
				while(j<length&&index<sortedDatabaseType.size()){
					Map<DatabaseType, Integer> count=new HashMap<DatabaseType,Integer>();
					List<DatabaseType> currentTypes=sortedDatabaseType.get(index);
					for(int k=0;k<currentTypes.size();k++){
						if(!count.containsKey(currentTypes.get(k)))count.put(currentTypes.get(k),0);
						count.put(currentTypes.get(k),count.get(currentTypes.get(k))+1);
					}
					for(int k=0;k<currentTypes.size();k++){
						if(j<result.length){
							for(Entry<DatabaseType,Integer> c:count.entrySet()){
								if(!result[j].containsKey(c.getKey()))result[j].put(c.getKey(),0.0);
								result[j].put(c.getKey(),result[j].get(c.getKey())+1.0*c.getValue()/currentTypes.size());
							}
							j++;
						}
					}
					index++;				
				}
			}
		}

		final DefaultCategoryDataset dataset = new DefaultCategoryDataset();
		File txtFile=new File(outputFile.getAbsolutePath().replaceAll(".jpg", ".txt"));
		BufferedWriter bw=new BufferedWriter(new FileWriter(txtFile));
		List<DatabaseType> header=new ArrayList<DatabaseType>();
		for(DatabaseType db:DatabaseType.values()){
			bw.write("\t"+db.toString());
			header.add(db);
		}
		bw.newLine();
		for(int i=0;i<result.length;i++){
			bw.write(Integer.toString(i));
			int sum=0;
			for(Entry<DatabaseType, Double> e:result[i].entrySet()){
				sum+=e.getValue();
			}
			double[] values=new double[header.size()];
			for(Entry<DatabaseType, Double> e:result[i].entrySet()){
				dataset.addValue(e.getValue(), e.getKey(), Integer.toString(i));
				values[header.indexOf(e.getKey())]=e.getValue()/sum;
			}
			for(double d:values)bw.write("\t"+d);
			bw.newLine();
		}
		final JFreeChart chart = ChartFactory.createBarChart("Boxplot","Ranks","Percentage",dataset);
		ChartUtilities.saveChartAsJPEG(outputFile, chart, Math.min(2000,length*100), 1000);
		bw.close();
		return result;
	}

	public HitStatistic getHitStatistic(){
		HitStatistic hs=new HitStatistic();
		hs.m=this;
		int[][] sortedIndizes=getIndizesOfSortedScores();
		for(int i=0;i<sortedIndizes.length;i++){
			Hit h=new Hit();
			h.massbankQuery=massbankQuery.get(i);
			h.hitsInDB=getHitsInDB(massbankQuery.get(i));
			if(!h.hitsInDB.isEmpty()){
				h.minRankBestHit=Integer.MAX_VALUE;
				h.maxRankBestHit=Integer.MAX_VALUE;
				h.scoreBestHit=Double.NEGATIVE_INFINITY;
			}
			int currentNumberHits=0;
			List<List<MassBank>> sortedEntries=getSortedMassbankEntries(i, sortedIndizes[i]);
			if(sortedEntries.size()>0){
				for(int j=0;j<sortedEntries.size();j++){
					List<MassBank> group=sortedEntries.get(j);
					if(h.bestNonHits.isEmpty()){
						for(MassBank mb:group){
							if(!h.hitsInDB.contains(mb)){
								if(h.scoreBestNonHit==null)h.scoreBestNonHit=Double.NEGATIVE_INFINITY;
								h.bestNonHits.add(mb);
								h.scoreBestNonHit=Math.max(h.scoreBestNonHit,getSimilarityValue(i,massbankDB.indexOf(mb)));
							}
						}
					}
					int numberCorrectHits=0;
					for(MassBank mb:group){
						if(h.hitsInDB.contains(mb)){
							numberCorrectHits++;
						}
					}
					for(MassBank mbDB:h.hitsInDB){				
						if(group.contains(mbDB)){
							h.minRankBestHit=Math.min(h.minRankBestHit,currentNumberHits+1);
							h.maxRankBestHit=Math.min(h.maxRankBestHit,currentNumberHits+1+(group.size()-numberCorrectHits));
							h.scoreBestHit=Math.max(h.scoreBestHit, getSimilarityValue(i,massbankDB.indexOf(mbDB)));
						}
					}				
					currentNumberHits+=group.size();
				}
			}
			h.numberEntries=currentNumberHits;
			hs.hits.add(h);
		}

		return hs;
	}

	//	public int[] getNumberOfHits(){
	//		int numberTrueHits=0;
	//		int numberPossibleHits=0;
	//		int numberNonemptyLists=0;
	//		int[][] sortedIndizes=getIndizesOfSortedScores();
	//		for(int i=0;i<sortedIndizes.length;i++){
	//			boolean firstEntry=true;
	//			boolean foundTrueHits=false;
	//			boolean foundPossibleHits=false;
	//			boolean foundNonempty=false;
	//			for(int j=0;j<sortedIndizes[i].length;j++){
	//				int currCol=sortedIndizes[i][j];
	//				double currValue=getSimilarityValue(i,currCol);
	//				if(Double.isNaN(currValue))continue;
	//				if(massbankQuery.get(i).hasEqualInChiKey(massbankDB.get(currCol))){					
	//					if(firstEntry)foundTrueHits=true;
	//					foundPossibleHits=true;
	//				}
	//				double nextValue=j+1<sortedIndizes[i].length?getSimilarityValue(i,sortedIndizes[i][j+1]):Double.NaN;
	//				if(!Double.isNaN(nextValue)&&currValue!=nextValue)firstEntry=false;
	//				foundNonempty=true;
	//				if(foundTrueHits&&foundPossibleHits&&foundNonempty)break;
	//			}
	//			if(foundTrueHits)numberTrueHits++;
	//			if(foundPossibleHits)numberPossibleHits++;
	//			if(foundNonempty)numberNonemptyLists++;
	//		}
	//		return new int[]{numberTrueHits, numberPossibleHits, numberNonemptyLists};
	//	}
}	
