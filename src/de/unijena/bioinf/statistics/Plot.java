package de.unijena.bioinf.statistics;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.statistics.DefaultBoxAndWhiskerCategoryDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import de.unijena.bioinf.decoy.model.MassBank;

public class Plot {

	public static void writeRankVSScore(File outputFile, String add, List<List<ResultList>>[] results, int length) throws Exception{

		Set<String>[] legendQuery=new Set[results.length];
		Set<String>[] legendDB=new Set[results.length];
		
		String name="RankVSscore";
		boolean isMean=false;
		
		for(int i=0;i<results.length;i++){
			if(results[i].size()>1)isMean=true;
			if(legendQuery[i]==null)legendQuery[i]=new HashSet<String>();
			if(legendDB[i]==null)legendDB[i]=new HashSet<String>();
			List<ResultList> resultLists=results[i].get(0);
			for(ResultList resultList:resultLists){
				legendQuery[i].add(resultList.query.getDB());
				for(Result r:resultList.results)legendDB[i].add(r.getDB());
			}
		}
		if(isMean)name+="Mean";

		String[] legend=new String[results.length];
		for(int i=0;i<legendQuery.length;i++){
			String tmp="";
			for(String s:legendQuery[i])tmp+=s+"_";
			tmp=tmp.substring(0,tmp.length()-1);
			tmp+="-";
			for(String s:legendDB[i])tmp+=s+"_";
			tmp=tmp.substring(0,tmp.length()-1);
			legend[i]=tmp;
		}

		final DefaultBoxAndWhiskerCategoryDataset dataset = new DefaultBoxAndWhiskerCategoryDataset();
		for (int i = 0; i < length; i++) {
			for (int j = 0; j < results.length; j++) {
				List<Double> scoreList=new ArrayList<Double>();
				for(List<ResultList> rls:results[j]){
					for(ResultList rl:rls){
						if(i<rl.results.size())scoreList.add(rl.results.get(i).score);
					}
				}
				dataset.add(scoreList, legend[j], i);
			}
		}

		String fileName="";
		for(String s:legend)fileName+=(s+"---");
		fileName=fileName.substring(0,fileName.length()-3);

		JFreeChart chart=ChartFactory.createBoxAndWhiskerChart("Whiskerplot", "Ranks", "Scores", dataset, true);
		ChartUtilities.saveChartAsJPEG(new File(outputFile.getPath()+"\\"+name+"_"+fileName+add+".jpg"), chart, Math.min(2000,length*100), 1000);
	}

	public static void writeScoreDistributionOfTopRank(File outputFile, String add, List<List<ResultList>>[] resultLists, double step) throws Exception{
		for(List<List<ResultList>> resultList:resultLists){
			Map<Double, Integer> topScoringMatchesAll=new TreeMap<Double,Integer>();
			Map<Double, Integer> topScoringMatchesTrue=new TreeMap<Double,Integer>();
			Map<Double, Integer> topScoringMatchesFalse=new TreeMap<Double,Integer>();
			int numberTrueMatches=0;
			int numberFalseMatches=0;
			String fileName=null;

			String name="ScoreDistributionTopRank";
			if(resultList.size()>1)name+="Mean";
			
			for(List<ResultList> result:resultList){
				if(fileName==null){
					for(ResultList rl:result){
						if(!rl.results.isEmpty())fileName=rl.query.getDB()+"-"+rl.results.get(0).getDB();
					}
				}

				for(ResultList rl:result){
					Collections.sort(rl.results);
					if(!rl.results.isEmpty()&&rl.results.get(0).massbank.inchi!=null&&rl.query.massbank.inchi!=null){
						double value=Math.round(rl.results.get(0).score/step)*step;
						if(!topScoringMatchesAll.containsKey(value))topScoringMatchesAll.put(value, 0);
						topScoringMatchesAll.put(value,topScoringMatchesAll.get(value)+1);
						if(rl.results.get(0).massbank.inchi.equals(rl.query.massbank.inchi)){
							if(!topScoringMatchesTrue.containsKey(value))topScoringMatchesTrue.put(value, 0);
							topScoringMatchesTrue.put(value,topScoringMatchesTrue.get(value)+1);
							numberTrueMatches++;
						}else{
							if(!topScoringMatchesFalse.containsKey(value))topScoringMatchesFalse.put(value, 0);
							topScoringMatchesFalse.put(value,topScoringMatchesFalse.get(value)+1);
							numberFalseMatches++;
						}
					}
				}
			}

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
			ChartUtilities.saveChartAsJPEG(new File(outputFile.getPath()+"\\"+name+"_"+fileName+add+".jpg"), chart, 1000, 400);
		}

	}

	public static void writeRankVSPercentage(File outputFile, String add, List<ResultList>[] resultsOriginal, List<List<ResultList>>[] resultsDecoy, int length) throws Exception{

		for(List<ResultList> o:resultsOriginal){
			for(List<List<ResultList>> ds:resultsDecoy){
				String name="RankVSpercentage";
				if(ds.size()>1)name+="Mean";
				Map<String,int[]> numbers=new TreeMap<String,int[]>();
				for(List<ResultList> d:ds){
					List<ResultList> merged=ResultList.mergeResults(new List[]{o,d});
					for(int i=0;i<merged.size();i++){
						ResultList rl=merged.get(i);
						for(int j=0;j<Math.min(length,rl.results.size());j++){
							String key=rl.query.getDB()+"-"+rl.results.get(j).getDB();
							if(!numbers.containsKey(key))numbers.put(key, new int[length]);
							numbers.get(key)[j]=numbers.get(key)[j]+1;
						}
					}
				}

				String fileName="";
				final DefaultCategoryDataset dataset = new DefaultCategoryDataset();
				for(Entry<String, int[]> n:numbers.entrySet()){
					double all=0;
					for(int i=0;i<n.getValue().length;i++){
						dataset.addValue(n.getValue()[i], n.getKey(), Integer.toString(i));
					}
					fileName+=n.getKey()+"---";
				}

				fileName=fileName.substring(0,fileName.length()-3);

				final JFreeChart chart = ChartFactory.createBarChart("Boxplot","Ranks","Percentage",dataset);
				ChartUtilities.saveChartAsJPEG(new File(outputFile.getPath()+"\\"+name+"_"+fileName+add+".jpg"), chart, Math.min(2000,length*100), 1000);
			}
		}
	}
	
//	public static void writeRankVSPercentageFilteredByInchi(File outputFile, String add, List<ResultList>[] resultsOriginal, List<List<ResultList>>[] resultsDecoy, int length) throws Exception{
//
//		for(List<ResultList> o:resultsOriginal){
//			for(List<List<ResultList>> ds:resultsDecoy){
//				String name="RankVSpercentageFiltered";
//				if(ds.size()>1)name+="Mean";
//				Map<String,int[]> numbers=new TreeMap<String,int[]>();
//				for(List<ResultList> d:ds){
//					List<ResultList> merged=ResultList.mergeResults(new List[]{o,d});
//					for(int i=0;i<merged.size();i++){
//						ResultList rl=new ResultList(merged.get(i));
//						rl.filterResultsByInChi();
//						merged.set(i,rl);
//						for(int j=0;j<Math.min(length,rl.results.size());j++){
//							String key=rl.query.getDB()+"-"+rl.results.get(j).getDB();
//							if(!numbers.containsKey(key))numbers.put(key, new int[length]);
//							numbers.get(key)[j]=numbers.get(key)[j]+1;
//						}
//					}
//				}
//
//				String fileName="";
//				final DefaultCategoryDataset dataset = new DefaultCategoryDataset();
//				for(Entry<String, int[]> n:numbers.entrySet()){
//					double all=0;
//					for(int i=0;i<n.getValue().length;i++){
//						dataset.addValue(n.getValue()[i], n.getKey(), Integer.toString(i));
//					}
//					fileName+=n.getKey()+"---";
//				}
//
//				fileName=fileName.substring(0,fileName.length()-3);
//
//				final JFreeChart chart = ChartFactory.createBarChart("Boxplot","Ranks","Percentage",dataset);
//				ChartUtilities.saveChartAsJPEG(new File(outputFile.getPath()+"\\"+name+"_"+fileName+add+".jpg"), chart, Math.min(2000,length*100), 1000);
//			}
//		}
//	}

	public static void writeEstimatedQValueVSCalculatedQValue(File outputFile, String add, List<ResultList>[] resultsOriginal, List<List<ResultList>>[] resultsDecoy) throws Exception{

		for(List<ResultList> o:resultsOriginal){
			for(List<List<ResultList>> ds:resultsDecoy){
				String name="qValuesComparison";
				if(ds.size()>1)name+="Mean";
				XYSeries series = null;
				String fileName=null;
				for(List<ResultList> d:ds){
					List<Result> topScores=new ArrayList<Result>();
					String fileNameLeft="";
					for(ResultList rl:o){
						if(!rl.results.isEmpty()){
							topScores.add(rl.results.get(0));
							fileNameLeft=rl.query.getDB()+"_"+rl.results.get(0).getDB();
						}
					}
					String fileNameRight="";
					for(ResultList rl:d){
						if(!rl.results.isEmpty()){
							topScores.add(rl.results.get(0));
							fileNameRight=rl.query.getDB()+"_"+rl.results.get(0).getDB();
						}
					}
					Collections.sort(topScores);
					if(fileName==null)fileName=fileNameLeft+"---"+fileNameRight;

					List<Double> FDRByDecoy=new ArrayList<Double>();
					List<Double> FDRCalculated=new ArrayList<Double>();
					int countTrueMatches=0;
					int countFalseMatches=0;
					int countDecoyMatches=0;
					for(Result r:topScores){
						if(r.massbank.isDecoy()){
							countDecoyMatches++;
						}else{
							if(r.isTrueMatch){
								countTrueMatches++;
							}else{
								countFalseMatches++;
							}
							FDRByDecoy.add(2.0*countDecoyMatches/(countDecoyMatches+countTrueMatches+countFalseMatches));
							FDRCalculated.add(1.0*countFalseMatches/(countFalseMatches+countTrueMatches));
						}			
					}

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

					if(series==null)series=new XYSeries("q-Values " + fileName);
					for(int i=0;i<FDRByDecoy.size();i++){
						series.add(FDRByDecoy.get(i), FDRCalculated.get(i));
					}

				}

				XYSeriesCollection dataset = new XYSeriesCollection(series);
				
				XYSeries seriesWH=new XYSeries("bisecting line");
				for(double i=0;i<=1;i+=0.05)seriesWH.add(i,i);
				dataset.addSeries(seriesWH);
				
				final JFreeChart chart = ChartFactory.createScatterPlot("q-Values", "q-values by decoy database", "calculated q-values", dataset);
				ChartUtilities.saveChartAsJPEG(new File(outputFile.getAbsolutePath()+"\\"+name+"_"+fileName+add+".jpg"), chart, 1000, 1000);				
			}
		}
	}
}
