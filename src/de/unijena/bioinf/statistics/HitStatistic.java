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
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import de.unijena.bioinf.decoy.model.MassBank;
import de.unijena.bioinf.statistics.SimilarityMatrix.DatabaseType;
import de.unijena.bioinf.statistics.SimilarityMatrix.MatchType;

public class HitStatistic {
	enum HitDescription{BestWorst, Best, NotFirstMatch, NotInDB, NotInQuery}
	enum MatchType{TruePositiveMatch, FalsePositiveMatch, DecoyMatch};
	enum HitMerge{Add, Merge};
	
	SimilarityMatrix m=null;
	
	Double averageRankMin=null;
	Double averageRankMedium=null;
	Double averageRankMax=null;
	Integer numberTopHitsMin=null;
	Integer numberTopHitsMax=null;
	Double numberTopHitsMedium=null;
	Integer numberPossibleHits=null;
	Integer numberSearches=null;
	Integer numberEntriesWithAtLeastOneHit=null;
	
	Double AUCMax=null;
	Double AUCAverage=null;
	Double AUCAverageSmall=null;
	Double AUCMin=null;
	Double F1ScoreMin=null;
	Double F1ScoreMax=null;
	
	List<Hit> hits=new ArrayList<Hit>();
	
	static final String sep=System.getProperty("file.separator");
	
	public HitStatistic(){}
	
	public HitStatistic(File inputFile, TreeMap<String, MassBank> massbankFilesOriginal) throws Exception{
		BufferedReader br=new BufferedReader(new FileReader(inputFile));
		String line=br.readLine();
		List<String> header=Arrays.asList(line.split("\t",Integer.MAX_VALUE));
		while((line=br.readLine())!=null){
			Hit h=new Hit();
			hits.add(h);
			String l[]=line.split("\t",Integer.MAX_VALUE);
			h.massbankQuery=massbankFilesOriginal.get(SimilarityMatrix.getIDDatasetChargeOrigin(l[header.indexOf("massbankID query")]));
			h.numberEntries=Integer.parseInt(l[header.indexOf("number entries")]);
			if(!l[header.indexOf("massbankIDs hits")].isEmpty()){
				for(String s:l[header.indexOf("massbankIDs hits")].split(";"))
					h.hitsInDB.add(massbankFilesOriginal.get(SimilarityMatrix.getIDDatasetChargeOrigin(s)));
			}
			if(!h.hitsInDB.isEmpty()){
				h.scoreBestHit=Double.parseDouble(l[header.indexOf("score best hit")]);
				h.minRankBestHit=Integer.parseInt(l[header.indexOf("min rank best hit")]);
				h.maxRankBestHit=Integer.parseInt(l[header.indexOf("max rank best hit")]);
			}
			if(!l[header.indexOf("massbankIDs non hits")].isEmpty()){
				for(String s:l[header.indexOf("massbankIDs non hits")].split(";"))
					h.bestNonHits.add(massbankFilesOriginal.get(SimilarityMatrix.getIDDatasetChargeOrigin(s)));
			}
			if(!h.bestNonHits.isEmpty()){
				h.scoreBestNonHit=Double.parseDouble(l[header.indexOf("score best non hit")]);
			}
		}
		br.close();
	}
	
	
	public static String getHitDescriptionString(HitDescription hd){
		switch (hd){
			case BestWorst: return "best case worst case";
			case Best: return "best case";
			case NotFirstMatch: return "not first match";
			case NotInDB: return "not in DB";
			case NotInQuery: return "not in query";
			default: return "";
		}	
	}
	
	public static int getHitDescriptionRank(HitDescription hd){
		switch (hd){
			case BestWorst: return 3;
			case Best: return 2;
			case NotFirstMatch: return 1;
			case NotInDB: return 0;
			case NotInQuery: return 0;
			default: return -1;
		}	
	}
	
	public static Integer compareHitDescription(HitDescription hd1, HitDescription hd2){
		int rank1=getHitDescriptionRank(hd1);
		int rank2=getHitDescriptionRank(hd2);
		if(rank1==0||rank2==0)return null;
		return Integer.compare(rank1,rank2);
	}
	
	public double getAverageRankMin(){
		if(averageRankMin!=null)return averageRankMin;
		double sum=0;
		int k=0;
		for(Hit h:hits){
			if(h.isEmpty()||h.noCorrectHit())continue;
			sum+=h.minRankBestHit;
			k++;
		}
		averageRankMin=sum/k;
		return averageRankMin;
	}
	
	public double getAverageRankMax(){
		if(averageRankMax!=null)return averageRankMax;
		double sum=0;
		int k=0;
		for(Hit h:hits){
			if(h.isEmpty()||h.noCorrectHit())continue;
			sum+=h.maxRankBestHit;
			k++;
		}
		averageRankMax=sum/k;
		return averageRankMax;
	}
	
	public double getAverageRankMedium(){
		if(averageRankMedium!=null)return averageRankMedium;
		double sum=0;
		int k=0;
		for(Hit h:hits){
			if(h.isEmpty()||h.noCorrectHit())continue;
			sum+=h.getMediumRankBestHit();
			k++;
		}
		averageRankMedium=sum/k;
		return averageRankMedium;
	}
	
	public double getNumberTopHitsMax(){
		if(numberTopHitsMax!=null)return numberTopHitsMax;
		int k=0;
		for(Hit h:hits){
			if(h.isEmpty()||h.noCorrectHit())continue;
			if(h.correctHitMinRank())k++;
		}
		numberTopHitsMax=k;
		return numberTopHitsMax;
	}
	
	public double getNumberTopHitsMin(){
		if(numberTopHitsMin!=null)return numberTopHitsMin;
		int k=0;
		for(Hit h:hits){
			if(h.isEmpty()||h.noCorrectHit())continue;
			if(h.correctHitMaxRank())k++;
		}
		numberTopHitsMin=k;
		return numberTopHitsMin;
	}
	
	public double getNumberTopHitsMedium(){
		if(numberTopHitsMedium!=null)return numberTopHitsMedium;
		double k=0;
		for(Hit h:hits){
			if(h.isEmpty()||h.noCorrectHit()||!h.correctHitMinRank())continue;	
			k+=h.correctHitMediumRank();
		}
		numberTopHitsMedium=k;
		return numberTopHitsMedium;
	}
	
	public int getNumberPossibleHits(){
		if(numberPossibleHits!=null)return numberPossibleHits;
		int k=0;
		for(Hit h:hits){
			if(h.isEmpty())continue;
			if(!h.noCorrectHit())k++;
		}
		numberPossibleHits=k;
		return numberPossibleHits;
	}
	
	public int getNumberEntriesWithAtLeastOneHit(){
		if(numberEntriesWithAtLeastOneHit!=null)return numberEntriesWithAtLeastOneHit;
		int k=0;
		for(Hit h:hits){
			if(h.bestNonHits.isEmpty()&&h.noCorrectHit())continue;
			k++;
		}
		numberEntriesWithAtLeastOneHit=k;
		return numberEntriesWithAtLeastOneHit;
	}
	
	public int getNumberSearches(){
		if(numberSearches!=null)return numberSearches;
		numberSearches=hits.size();
		return numberSearches;
	}
	
	public String toString(){
		String result="";
		result+="ranks:\t"+getAverageRankMin()+"\t"+getAverageRankMedium()+"\t"+getAverageRankMax()+"\n";
		result+="top hits:\t"+getNumberTopHitsMin()+"\t"+getNumberTopHitsMax()+"\n";
		result+="possible hits:\t"+getNumberPossibleHits()+"\n";
		result+="searches:\t"+getNumberSearches()+"\n";
		return result; 
	}
	
	public static String getHeaderString(){
		return "number searches\tpossible hits\tnumber top hits min\tnumber top hits medium\tnumber top hits max\taverage rank min\taverage rank medium\taverage rank max\tAUC min\tAUC average\tAUC max\tAUC average 0.1\tF1 min\tF1 max";
	}
	
	public String toEntryString(){
		String result="";
		result+=getNumberSearches()+"\t";
		result+=getNumberPossibleHits()+"\t";
		result+=getNumberTopHitsMin()+"\t"+getNumberTopHitsMedium()+"\t"+getNumberTopHitsMax()+"\t";
		result+=getAverageRankMin()+"\t"+getAverageRankMedium()+"\t"+getAverageRankMax()+"\t";		
		result+=getAUCMin()+"\t"+getAUCAverage()+"\t"+getAUCMax()+"\t";
		result+=getAUCAverageSmall()+"\t";
		result+=getF1ScoreMin()+"\t"+getF1ScoreMax();
		return result; 
	}
	
	public double getAUCMax(){
		if(AUCMax!=null)return AUCMax;
		List<Hit> hits=new ArrayList<Hit>();
		for(Hit h:this.hits)if(!h.isEmpty())hits.add(h);
		Collections.sort(hits,new HitComparator(HitComparator.SortAtEqualScore.BESTCASE));
		AUCMax=getAUC(hits);
		return AUCMax;
	}
	
	public double getAUCAverage(){
		if(AUCAverage!=null)return AUCAverage;
		AUCAverage=getAUCAverageAll(Double.POSITIVE_INFINITY);
		return AUCAverage;
	}
	
	public double getAUCAverageSmall(){
		if(AUCAverageSmall!=null)return AUCAverageSmall;
		AUCAverageSmall=getAUCAverageAll(0.1);
		return AUCAverageSmall;
	}
	
	public double getAUCAverageAll(double fpr){
		List<Hit> hits=new ArrayList<Hit>();
		for(Hit h:this.hits)if(!h.isEmpty())hits.add(h);
		Collections.sort(hits,new HitComparator(HitComparator.SortAtEqualScore.BESTCASE));
		return getAUCAverage(hits,fpr);
	}
	
	public List<double[]> getROCMax(){		
		List<Hit> hits=new ArrayList<Hit>();
		for(Hit h:this.hits)if(!h.isEmpty())hits.add(h);
		Collections.sort(hits,new HitComparator(HitComparator.SortAtEqualScore.BESTCASE));		
		return getROC(hits);
	}
	
	public List<double[]> getROCAverage(){		
		List<Hit> hits=new ArrayList<Hit>();
		for(Hit h:this.hits)if(!h.isEmpty())hits.add(h);
		Collections.sort(hits,new HitComparator(HitComparator.SortAtEqualScore.BESTCASE));		
		return getROCAverage(hits);
	}
	
	public double getAUCMin(){
		if(AUCMin!=null)return AUCMin;
		List<Hit> hits=new ArrayList<Hit>();
		for(Hit h:this.hits)if(!h.isEmpty())hits.add(h);
		Collections.sort(hits,new HitComparator(HitComparator.SortAtEqualScore.WORSTCASE));
		AUCMin=getAUC(hits);
		return AUCMin;
	}
	
	public double getF1ScoreMin(){
		if(F1ScoreMin!=null)return F1ScoreMin;
		List<Hit> hits=new ArrayList<Hit>();
		for(Hit h:this.hits)if(!h.isEmpty())hits.add(h);
		Collections.sort(hits,new HitComparator(HitComparator.SortAtEqualScore.WORSTCASE));
		F1ScoreMin=getF1Score(hits);		
		return F1ScoreMin;
	}
	
	
	public double getF1ScoreMax(){
		if(F1ScoreMax!=null)return F1ScoreMax;
		List<Hit> hits=new ArrayList<Hit>();
		for(Hit h:this.hits)if(!h.isEmpty())hits.add(h);
		Collections.sort(hits,new HitComparator(HitComparator.SortAtEqualScore.BESTCASE));
		F1ScoreMax=getF1Score(hits);
		return F1ScoreMax;
	}
	
	public double getF1Score(List<Hit> hits){
		int allPositives=0;
		for(Hit h:hits){
			if(h.correctHitMinRank())allPositives++;
		}
		int currentNumberTruePositives=0;
		int currentNumberFalsePositives=0;
		double F=0;
		int countF=0;		
		for(Hit h:hits){
			if(h.correctHitMinRank())currentNumberTruePositives++;
			else currentNumberFalsePositives++;
			double sensitivity=1.0*currentNumberTruePositives/allPositives;
			double precision=1.0*currentNumberTruePositives/(currentNumberTruePositives+currentNumberFalsePositives);
			if(precision==0&&sensitivity==0){
				F+=0;
			}else{
				F+=2*precision*sensitivity/(precision+sensitivity);
			}
			countF++;
		}
		return F/countF;
	}

	private double getAUC(List<Hit> hits) {
		List<Integer> numbers=new ArrayList<Integer>();
		numbers.add(0);
		for(Hit h:hits){
			if(h.correctHitMinRank())numbers.add(numbers.get(numbers.size()-1)+1);
			else numbers.add(numbers.get(numbers.size()-1));
		}
		double result=0;
		int correctHits=0;
		int inCorrectHits=0;
		for(int i=1;i<numbers.size();i++){
			if(numbers.get(i-1)==numbers.get(i)){
				result+=numbers.get(i-1);
				inCorrectHits++;
			}else{
				correctHits++;
			}
		}
		result/=(inCorrectHits*correctHits);
		return result;
	}
	
	private double getAUCAverage(List<Hit> hits, double fpr) {
		List<double[]> numbers=getROCAverage(hits);
		double result=0;
		double maxX=numbers.get(numbers.size()-1)[0];
		double maxY=numbers.get(numbers.size()-1)[1];
		for(int i=1;i<numbers.size();i++){
			if(numbers.get(i)[0]/maxX<fpr){
				double deltaX=numbers.get(i)[0]-numbers.get(i-1)[0];
				double deltaY=numbers.get(i)[1]-numbers.get(i-1)[1];
				result+=deltaX*deltaY/2+numbers.get(i-1)[1]*deltaX;				
			}else{
				double deltaX=fpr*maxX-numbers.get(i-1)[0];
				double deltaY=numbers.get(i)[1]-numbers.get(i-1)[1];
				result+=deltaX*deltaY/2+numbers.get(i-1)[1]*deltaX;
				break;
			}
		}
		result/=(maxX*maxY);
		return result;
	}
	
	private List<double[]> getROCAverage(List<Hit> hits) {
		List<List<Hit>> hitsTmp=new ArrayList<List<Hit>>();
		List<Hit> hTmp=new ArrayList<Hit>();
		for(int i=0;i<hits.size();i++){
			if(i==0||hits.get(i).getBestScore()!=hits.get(i-1).getBestScore()){
				hTmp=new ArrayList<Hit>();
				hitsTmp.add(hTmp);
			}
			hTmp.add(hits.get(i));
		}
		
		List<double[]> numbers=new ArrayList<double[]>();
		numbers.add(new double[]{0.0,0.0});
		
		double numberCorrectHits=0;
		double numberWrongHits=0;
		for(List<Hit> hitList:hitsTmp){			
			for(Hit h:hitList){
				double correctHits=0;
				if(h.correctHitMinRank())correctHits=h.correctHitMediumRank();
				numberCorrectHits+=correctHits;
				numberWrongHits+=1-correctHits;
			}
			numbers.add(new double[]{numberWrongHits, numberCorrectHits});
		}
		double maxX=numbers.get(numbers.size()-1)[0];
		double maxY=numbers.get(numbers.size()-1)[1];
		for(int i=0;i<numbers.size();i++){
			numbers.get(i)[0]=numbers.get(i)[0]/maxX;
			numbers.get(i)[1]=numbers.get(i)[1]/maxY;
		}
		
		return numbers;
	}
	
	private List<double[]> getROC(List<Hit> hits) {
		List<double[]> result=new ArrayList<double[]>();
		int an=0;
		int ap=0;
		for(Hit h:hits){
			if(h.correctHitMinRank())ap++;
			else an++;
		}
		int fp=0;		
		int rp=0;
		result.add(new double[]{0.0,0.0});
		for(Hit h:hits){
			if(h.correctHitMinRank())rp++;
			else fp++;
			result.add(new double[]{1.0*fp/(an),1.0*rp/(ap)});
		}		
		return result;
	}
	
	public void writeHitStatisticToFile(File outputFile) throws Exception{
		Map<MassBank, HitDescription> result=getHitDescription();
		Map<MassBank, Hit> massbank2Hit=new HashMap<MassBank, Hit>();
		for(Hit h:hits){
			massbank2Hit.put(h.massbankQuery, h);			
		}
		
		BufferedWriter bw=new BufferedWriter(new FileWriter(outputFile));
		bw.write("massbankID query\tnumber entries\tmassbankIDs hits\tscore best hit\tmin rank best hit\tmax rank best hit\tmassbankIDs non hits\tscore best non hit\t hit description");
		bw.newLine();
		for(Entry<MassBank, HitDescription> e:result.entrySet()){			
			Hit h=massbank2Hit.get(e.getKey());
			String mbHits="";
			for(int i=0;i<h.hitsInDB.size();i++){
				if(i!=0)mbHits+=";";
				mbHits+=h.hitsInDB.get(i).massbankID;
			}
			
			String mbNonHits="";
			for(int i=0;i<h.bestNonHits.size();i++){
				if(i!=0)mbNonHits+=";";
				mbNonHits+=h.bestNonHits.get(i).massbankID;
			}
			
			bw.write(
			h.massbankQuery.massbankID+"\t"+
			h.numberEntries+"\t"+
			mbHits+"\t"+
			h.scoreBestHit+"\t"+
			h.minRankBestHit+"\t"+
			h.maxRankBestHit+"\t"+
			mbNonHits+"\t"+
			h.scoreBestNonHit+"\t"+
			getHitDescriptionString(e.getValue()));
			bw.newLine();
		}
		bw.close();
	}
	
	public Map<MassBank, HitDescription> getHitDescription() throws Exception{
		Map<MassBank, HitDescription> result=new TreeMap<MassBank, HitDescription>();
		for(Hit h:hits){
			HitDescription description=null;
			if(h.noCorrectHit())description=HitDescription.NotInDB;
			else{
				if(h.correctHitMinRank())description=HitDescription.Best;
				else description=HitDescription.NotFirstMatch;
				if(h.correctHitMaxRank()) description=HitDescription.BestWorst;
			}
			result.put(h.massbankQuery,description);
		}
		return result;
	}
	
	public static void writeROCCurves(File outputFile, TreeMap<String, HitStatistic> hits) throws Exception{

		XYSeriesCollection dataset = new XYSeriesCollection();
		
		for(Entry<String, HitStatistic> e:hits.entrySet()){
			File txtFile=new File(outputFile.getAbsolutePath().replaceAll(".jpg","_"+e.getKey().split("-")[2]+".txt"));
			BufferedWriter bw=new BufferedWriter(new FileWriter(txtFile));
			XYSeries series= new XYSeries(e.getKey());
			bw.write("false positive rate\tsensitivity");
			bw.newLine();
			List<double[]> roc=e.getValue().getROCAverage();
			for(double[] r:roc){
				bw.write(r[0]+"\t"+r[1]);
				bw.newLine();
				series.add(r[0],r[1]);
			}
			dataset.addSeries(series);
			bw.close();
		}
		
		
		final JFreeChart chart =ChartFactory.createXYLineChart("ROCCurve",  "false positive rate", "sensitivity", dataset);
		ChartUtilities.saveChartAsJPEG(outputFile, chart, 1000, 400);
	}
	
	public void writeScoreDistributionOfTopRank(File outputFile, double step, String add) throws Exception{
		File jpg=new File(outputFile.getPath()+sep+add+"ScoreDistributionTopRank"+"_"+m.getMethodsQueryAndDBString()+".jpg");
		File txt=new File(outputFile.getPath()+sep+add+"ScoreDistributionTopRank"+"_"+m.getMethodsQueryAndDBString()+".txt");
		
		Map<Double, Integer> topScoringMatchesAll=new TreeMap<Double,Integer>();
		Map<Double, Integer> topScoringMatchesTrue=new TreeMap<Double,Integer>();
		Map<Double, Integer> topScoringMatchesFalse=new TreeMap<Double,Integer>();
		int numberTrueMatches=0;
		int numberFalseMatches=0;
		
		double maxValue=Double.NEGATIVE_INFINITY;
		for(Hit h:hits){
			if(h.isEmpty())continue;
			Double valueOrig=h.getBestScore();		
			maxValue=Math.max(Math.round(valueOrig/step)*step,maxValue);
		}

		for(Hit h:hits){
			if(h.isEmpty())continue;
			Double valueOrig=h.getBestScore();
			valueOrig/=maxValue;
			double value=Math.round(valueOrig/step)*step;
			if(!topScoringMatchesAll.containsKey(value))topScoringMatchesAll.put(value, 0);
			topScoringMatchesAll.put(value,topScoringMatchesAll.get(value)+1);
			if(h.correctHitMinRank()){
				if(!topScoringMatchesTrue.containsKey(value))topScoringMatchesTrue.put(value, 0);
				topScoringMatchesTrue.put(value,topScoringMatchesTrue.get(value)+1);
				numberTrueMatches++;
			}else{
				if(!topScoringMatchesFalse.containsKey(value))topScoringMatchesFalse.put(value, 0);
				topScoringMatchesFalse.put(value,topScoringMatchesFalse.get(value)+1);
				numberFalseMatches++;
			}
		}

		BufferedWriter bw=new BufferedWriter(new FileWriter(txt));
		Set<Double> v=new TreeSet<Double>();
		bw.write("\t"+m.getMethodsQueryAndDBString()+"-true\t"+m.getMethodsQueryAndDBString()+"-wrong");
		bw.newLine();
		v.addAll(topScoringMatchesTrue.keySet());
		v.addAll(topScoringMatchesFalse.keySet());
		for(Double d:v){
			bw.write(d+"\t");
			Integer t=topScoringMatchesTrue.get(d);
			Integer f=topScoringMatchesFalse.get(d);
			bw.write((t!=null?1.0*t/numberTrueMatches:"")+"\t");
			bw.write(""+(f!=null?1.0*f/numberFalseMatches:""));
			bw.newLine();
		}
		bw.close();
		
		
		XYSeries seriesTrue = new XYSeries("True Hits ("+numberTrueMatches+")");
		for(Entry<Double,Integer> m:topScoringMatchesTrue.entrySet()){
			seriesTrue.add(m.getKey(),new Double(1.0*m.getValue()/numberTrueMatches));
		}
		XYSeries seriesFalse = new XYSeries("False Hits ("+numberFalseMatches+")");
		for(Entry<Double,Integer> m:topScoringMatchesFalse.entrySet()){
			seriesFalse.add(m.getKey(),new Double(1.0*m.getValue()/numberFalseMatches));
		}
		XYSeries seriesAll= new XYSeries("All Hits ("+(numberTrueMatches+numberFalseMatches)+")");
		for(Entry<Double,Integer> m:topScoringMatchesAll.entrySet()){
			seriesAll.add(m.getKey(),new Double(1.0*m.getValue()/(numberTrueMatches+numberFalseMatches)));
		}

		XYSeriesCollection dataset = new XYSeriesCollection();
		dataset.addSeries(seriesTrue);
		dataset.addSeries(seriesFalse);
		dataset.addSeries(seriesAll);

		final JFreeChart chart =ChartFactory.createXYLineChart("Score Distribution", "score", "percentage", dataset);
		ChartUtilities.saveChartAsJPEG(jpg, chart, 1000, 400);
	}
	
	public static List<Hit> mergeHits(HitStatistic hsOriginal, HitStatistic hsDecoy, HitMerge hm){
		

		for(Hit h:hsOriginal.hits)h.db=Hit.DataBase.Original;
		for(Hit h:hsDecoy.hits)h.db=Hit.DataBase.Decoy;
		
		List<Hit> tmp=new  ArrayList<Hit>();
		for(Hit h:hsOriginal.hits)if(!h.isEmpty()){
			tmp.add(h);
		}
		System.out.println();
		for(Hit h:hsDecoy.hits)if(!h.isEmpty()){
			tmp.add(h);
		}
		
		List<Hit> mergedHits=new  ArrayList<Hit>();
		if(hm.equals(HitMerge.Add)){
			mergedHits=tmp;
		}else if(hm.equals(HitMerge.Merge)){
			Map<MassBank, Hit> mb=new HashMap<MassBank, Hit>();
			for(Hit h:tmp){
				if(!mb.containsKey(h.massbankQuery))mb.put(h.massbankQuery,h);
				else{
					Hit h2=mb.get(h.massbankQuery);
					if(h.getBestScore()>h2.getBestScore()||(h.getBestScore()==h2.getBestScore()&&h.db.equals(Hit.DataBase.Original))){
						mb.put(h.massbankQuery,h);
					}
				}
			}
			mergedHits.addAll(mb.values());
		}
		return mergedHits;
	}
	
	public static void writeEstimatedQValueVSCalculatedQValueToFile(File outputFile,File outputFileMean, File outputFileHitlist, HitStatistic hsOriginal, HitStatistic hsDecoy, HitMerge hm,  boolean FDR, List<MassBank> compoundsOfInterest) throws IOException{
		
		for(Hit h:hsOriginal.hits)h.db=Hit.DataBase.Original;
		for(Hit h:hsDecoy.hits)h.db=Hit.DataBase.Decoy;
		
		List<Hit> mergedHits=mergeHits(hsOriginal, hsDecoy, hm);
		Collections.sort(mergedHits,new HitComparator(HitComparator.SortAtEqualScore.BESTCASE));
		
		List<Double> bestScores=new ArrayList<Double>();
		List<MatchType> matches=new ArrayList<MatchType>();
		
		BufferedWriter bwHitlist=new BufferedWriter(new FileWriter(outputFileHitlist));
		if(compoundsOfInterest!=null){
			for(int i=0;i<compoundsOfInterest.size();i++){				
				if(i!=0)bwHitlist.write(";");
				bwHitlist.write(compoundsOfInterest.get(i).massbankID);					
			}
			bwHitlist.newLine();
		}
		
		bwHitlist.write("massbank\tmatch type\tscore");
		bwHitlist.newLine();
		
		List<Double> FDRByDecoy=new ArrayList<Double>();
		List<Double> FDRCalculated=new ArrayList<Double>();
		
		int allTrueMatches=0;
		int allFalseMatches=0;
		for(Hit h:mergedHits){
			if(compoundsOfInterest!=null&&!compoundsOfInterest.contains(h.massbankQuery))continue;
			if(h.isEmpty())continue;
			bestScores.add(h.getBestScore());
			MatchType mt;
			if(h.db.equals(Hit.DataBase.Decoy))mt=MatchType.DecoyMatch;
			else if(h.correctHitMinRank()){
				mt=MatchType.TruePositiveMatch;
				allTrueMatches++;
			}else{
				mt=MatchType.FalsePositiveMatch;
				allFalseMatches++;
			}
			matches.add(mt);
			bwHitlist.write(h.massbankQuery.massbankID+"\t"+mt+"\t"+h.getBestScore());
			bwHitlist.newLine();
		}
		bwHitlist.close();

		int countTrueMatches=0;
		int countFalseMatches=0;
		int countDecoyMatches=0;
		for(int i=0;i<matches.size();i++){
			if(matches.get(i)==MatchType.DecoyMatch){
				countDecoyMatches++;
			}else{
				if(matches.get(i)==MatchType.TruePositiveMatch){
					countTrueMatches++;
				}else{
					countFalseMatches++;
				}
				double fdrDecoy=hm.equals(HitMerge.Merge)?
								2.0*countDecoyMatches/(countDecoyMatches+countTrueMatches+countFalseMatches):
								1.0*countDecoyMatches/(countTrueMatches+countFalseMatches)*allFalseMatches/(allFalseMatches+allTrueMatches);
				FDRByDecoy.add(fdrDecoy);
				FDRCalculated.add(1.0*countFalseMatches/(countFalseMatches+countTrueMatches));
			}
		}
		double fdrDecoy=
			hm.equals(HitMerge.Merge)?
				2.0*countDecoyMatches/(countDecoyMatches+countTrueMatches+countFalseMatches):
				1.0*countDecoyMatches/(countTrueMatches+countFalseMatches)*allFalseMatches/(allFalseMatches+allTrueMatches);
		FDRByDecoy.add(fdrDecoy);			
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

		BufferedWriter bw=new BufferedWriter(new FileWriter(outputFile));
		BufferedWriter bwMean=new BufferedWriter(new FileWriter(outputFileMean));
		if(compoundsOfInterest!=null){
			for(int i=0;i<compoundsOfInterest.size();i++){				
				if(i!=0){
					bw.write(";");
					bwMean.write(";");
				}
				bw.write(compoundsOfInterest.get(i).massbankID);
				bwMean.write(compoundsOfInterest.get(i).massbankID);
			}
			bw.newLine();
			bwMean.newLine();
		}
		
		bw.write("calculated qValue\testimated qValue");
		bwMean.write("calculated qValue\testimated qValue");
		bw.newLine();
		bwMean.newLine();
		for(int i=0;i<FDRCalculated.size();i++){
			bw.write(FDRCalculated.get(i)+"\t"+FDRByDecoy.get(i));
			bw.newLine();
		}
		bw.close();
		
		Map<Double, List<Double>> means=new TreeMap<Double, List<Double>>();
		for(int i=0;i<FDRCalculated.size();i++){
			if(!means.containsKey(FDRCalculated.get(i)))means.put(FDRCalculated.get(i), new ArrayList<Double>());
			means.get(FDRCalculated.get(i)).add(FDRByDecoy.get(i));			
		}
		
		for(Entry<Double, List<Double>> e: means.entrySet()){
			double m=0;
			for(double d:e.getValue())m+=d;
			m/=e.getValue().size();
			bwMean.write(e.getKey()+"\t"+m);
			bwMean.newLine();
		}
		bwMean.close();		
	}
	
	public static void writeQValuesDecoyHits(File decoyFolder, File outputFile, HitStatistic hsOriginal, List<MassBank> compoundsOfInterest, boolean FDR) throws Exception{
		if(!decoyFolder.exists())decoyFolder.mkdir();
		List<Hit> currentHits=new ArrayList<Hit>();
		for(Hit h:hsOriginal.hits)
			if((compoundsOfInterest==null||compoundsOfInterest.contains(h.massbankQuery))&&!h.isEmpty()){
				h.db=Hit.DataBase.Original;
				currentHits.add(h);
			}
		
		Collections.sort(currentHits,new HitComparator(HitComparator.SortAtEqualScore.BESTCASE));
		List<MatchType> matches=new ArrayList<MatchType>();
		List<Double> FDRCalculated=new ArrayList<Double>();
		int allTrueMatches=0;
		int allFalseMatches=0;
		for(Hit h:currentHits){
			if(compoundsOfInterest!=null&&!compoundsOfInterest.contains(h.massbankQuery))continue;
			if(h.isEmpty())continue;
			MatchType mt;
			if(h.db.equals(Hit.DataBase.Decoy))mt=MatchType.DecoyMatch;
			else if(h.correctHitMinRank()){
				mt=MatchType.TruePositiveMatch;
				allTrueMatches++;
			}else{
				mt=MatchType.FalsePositiveMatch;
				allFalseMatches++;
			}
			matches.add(mt);
		}

		int countTrueMatches=0;
		int countFalseMatches=0;
		for(int i=0;i<matches.size();i++){
			if(matches.get(i)==MatchType.TruePositiveMatch){
				countTrueMatches++;
			}else{
				countFalseMatches++;
			}				
			FDRCalculated.add(1.0*countFalseMatches/(countFalseMatches+countTrueMatches));
		}		
		FDRCalculated.add(1.0*countFalseMatches/(countFalseMatches+countTrueMatches));

		if(!FDR){		
			double min=Double.POSITIVE_INFINITY;
			for(int i=FDRCalculated.size()-1;i>=0;i--){
				min=Math.min(FDRCalculated.get(i),min);
				FDRCalculated.set(i, min);
			}
		}
		
		int n=currentHits.size();
		File randomDecoyFile=new File(decoyFolder.getAbsolutePath()+sep+"RandomDecoyValues_"+n+".txt");
		if(!randomDecoyFile.exists()){
			createRandomDecoyQValues(n,randomDecoyFile);
		}
		List<double[]> avg=getAverageAndStandardDeviation(randomDecoyFile);
		BufferedWriter bw=new BufferedWriter(new FileWriter(outputFile));
		bw.write("calculated qValue\testimated qValue");
		bw.newLine();
		for(int i=0;i<Math.max(avg.size(), FDRCalculated.size());i++){
			double v1=i>=FDRCalculated.size()?1:FDRCalculated.get(i);
			double v2=i>=avg.size()?1:avg.get(i)[0];
			bw.write(v1+"\t"+calculateEstimatedQVAlue(v2,allFalseMatches,allTrueMatches));
			bw.newLine();
		}
		bw.close();
	}
	
	private static double calculateEstimatedQVAlue(double value, int allFalseMatches, int allTrueMatches){
		return value*allFalseMatches/(allFalseMatches+allTrueMatches);
	}
	
	public static List<double[]> getAverageAndStandardDeviation(File randomDecoyFile) throws IOException{
		List<double[]> result=new ArrayList<double[]>();
		BufferedReader br=new BufferedReader(new FileReader(randomDecoyFile));
		List<String> header=Arrays.asList(br.readLine().split("\t"));
		List<Double> qValues=new ArrayList<Double>();
		for(int i=1;i<header.size();i++){
			qValues.add(Double.parseDouble(header.get(i)));
		}		
		int n=0;
		while(br.readLine()!=null){
			n++;
		}
		br.close();
		
		double[][] probs=new double[n][qValues.size()];
		br=new BufferedReader(new FileReader(randomDecoyFile));
		br.readLine();
		String line;
		n=0;
		while((line=br.readLine())!=null){
			String[] l=line.split("\t");
			for(int i=1;i<l.length;i++){
				probs[n][i-1]=l[i].equals("null")?0.0:Double.parseDouble(l[i]);
			}
			n++;
		}
		
		for(int i=0;i<probs.length;i++){
			double avg=0;
			double sum=0;
			for(int j=0;j<probs[i].length;j++){
				avg+=qValues.get(j)*probs[i][j];
				sum+=probs[i][j];
			}
			avg/=sum;
			
			double s=0;
			sum=0;
			for(int j=0;j<probs[i].length;j++){
				s+=Math.pow(qValues.get(j),2)*probs[i][j];
				sum+=probs[i][j];
			}
			s=Math.pow(s/sum-Math.pow(avg,2),0.5);
			result.add(new double[]{avg, s});
		}
		br.close();
		
//		BufferedWriter bw=new BufferedWriter(new FileWriter(outputFile));
//		bw.write("\testimated qValue\tstandard deviation");
//		bw.newLine();
//		for(int i=0;i<result.size();i++){
//			bw.write(i+"\t"+result.get(i)[0]+"\t"+result.get(i)[1]);
//			bw.newLine();
//		}
//		bw.close();
		return result;
	}
	
	private static void outputMatrix(double[][] m, File outputFile) throws IOException{
		BufferedWriter bw=new BufferedWriter(new FileWriter(outputFile));
		for(int j=0;j<m[0].length;j++){
			bw.write("\t"+j);
		}
		bw.newLine();
		for(int i=0;i<m.length;i++){
			bw.write(Integer.toString(i));
			for(int j=0;j<m[i].length;j++){
				bw.write("\t"+m[i][j]);			
			}
			bw.newLine();
		}
		
		
		bw.close();
		
	}

	static void createRandomDecoyQValues(int n, File randomDecoyFile) throws IOException {
		double[][] P=new double[n+1][n+1];
		P[0][0]=1;
		for(int i=0;i<P.length;i++){
			for(int k=0;k<P[i].length;k++){
				if(i>0)P[i][k]+=P[i-1][k]*(n-i+1)/(2*n-i-k+1);
				if(k>0)P[i][k]+=P[i][k-1]*(n-k+1)/(2*n-i-k+1);
			}
		}
		
		double[][] D=new double[n][n+1];
		for(int i=0;i<D.length;i++){
			for(int k=0;k<D[i].length;k++){
//				D[i][k]=P[i][k]*P[n-i-1][n-k];
				D[i][k]=P[i][k]*(n-i)/(2*n-i-k);
			}
		}
		
		double[][] FDR=new double[n][n+1];
		for(int i=0;i<FDR.length;i++){
			for(int k=0;k<FDR[i].length;k++){
				FDR[i][k]=1.0*k/(i+1);
			}
		}
		
		Map<Double, Double>[] valuesANDprobs=new Map[D.length];
		for(int i=FDR.length-1;i>=0;i--){
			long time=System.currentTimeMillis();
			Map<Double, Double> currentMap=new HashMap<Double, Double>();
			valuesANDprobs[i]=currentMap;
			Map<Double, Double> previousMap=new HashMap<Double, Double>();
			if(i==FDR.length-1){
				previousMap.put(Double.POSITIVE_INFINITY,1.0);
			}else{
				previousMap=valuesANDprobs[i+1];
			}
			for(Entry<Double,Double> entry:previousMap.entrySet()){
				double d1=entry.getKey();
				double prob1=entry.getValue();
				for(int j=0;j<D[i].length;j++){
					double d2=FDR[i][j];
					double prob2=D[i][j];
					double min=Math.round(Math.min(d1,d2)*1000)/1000.0;
					double currentProb=currentMap.containsKey(min)?currentMap.get(min):0.0;
					currentMap.put(min, currentProb+prob1*prob2);
				}
			}
			System.out.print(i+" "+currentMap.size()+" "+(System.currentTimeMillis()-time)/1000.0);
			System.out.println();
		}
		
		BufferedWriter bw=new BufferedWriter(new FileWriter(randomDecoyFile));
		Set<Double> allKeys=new TreeSet<Double>();
		for(Map<Double,Double> v:valuesANDprobs){
			allKeys.addAll(v.keySet());
		}
		for(double d:allKeys)bw.write("\t"+d);
		bw.newLine();
		for(int i=0;i<valuesANDprobs.length;i++){
			Map<Double,Double> v=valuesANDprobs[i];
			bw.write(Integer.toString(i));
			for(double d:allKeys)
				bw.write("\t"+v.get(d));
			bw.newLine();
		}
		
		
		bw.close();
	}
	
}

class Hit{
	enum DataBase{Original, Decoy}
	
	DataBase db;
	
	MassBank massbankQuery;
	int numberEntries;
	
	List<MassBank> hitsInDB=new ArrayList<MassBank>();
	Double scoreBestHit;
	Integer minRankBestHit;
	Integer maxRankBestHit;		
	
	List<MassBank> bestNonHits=new ArrayList<MassBank>();
	Double scoreBestNonHit;
	
	
//	List<MassBank> bestHits=null;
		
	public boolean isEmpty(){
		return numberEntries==0;
	}
	
	public boolean noCorrectHit(){
		return hitsInDB.isEmpty();
	}
	
	public Double getBestScore(){
		if(isEmpty())return null;
		return Math.max(scoreBestHit==null?Double.MIN_VALUE:scoreBestHit, scoreBestNonHit==null?Double.MIN_VALUE:scoreBestNonHit);
	}
	
	public Boolean correctHitMinRank(){
		if(noCorrectHit())return false;
		return minRankBestHit==1;
	}
	
	public Double correctHitMediumRank(){
		if(noCorrectHit()||!correctHitMinRank())return null;
		return 1.0/maxRankBestHit;
	}
	
	public Double getMediumRankBestHit(){
		if(noCorrectHit())return null;
		return minRankBestHit+(maxRankBestHit-minRankBestHit)/2.0;
	}
	
	public Boolean correctHitMaxRank(){
		if(noCorrectHit())return false;
		return maxRankBestHit==1;
	}
}

class HitComparator implements Comparator<Hit>{
	
	enum SortAtEqualScore{BESTCASE, WORSTCASE};
	
	SortAtEqualScore s=null;
	public HitComparator(SortAtEqualScore s){
		this.s=s;
	}
	
	@Override
	public int compare(Hit o1, Hit o2) {
		int c=o2.getBestScore().compareTo(o1.getBestScore());
		if(c!=0)return c;
		if(o1.db!=null&&o2.db!=null&&o1.db!=o2.db){
			if(o2.db==Hit.DataBase.Decoy)return -1;
			return 1;
		}
		c=Boolean.compare(o2.correctHitMinRank(), o1.correctHitMinRank());
		if(s.equals(SortAtEqualScore.BESTCASE)){
			return c;
		}
		return -c;
		
	}
	
}
