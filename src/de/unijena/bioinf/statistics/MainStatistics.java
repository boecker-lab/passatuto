package de.unijena.bioinf.statistics;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Shape;
import java.awt.Stroke;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import javax.swing.filechooser.FileFilter;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.DeviationRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.statistics.DefaultBoxAndWhiskerCategoryDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import de.unijena.bioinf.decoy.model.MassBank;
import de.unijena.bioinf.spectralcomparison.ParametersSpectralComparison.COMPARISONMETHOD;
import de.unijena.bioinf.spectralcomparison.ParametersSpectralComparison.COMPARISONPOSTPROCESS;
import de.unijena.bioinf.statistics.HitStatistic.HitDescription;
import de.unijena.bioinf.statistics.SimilarityMatrix.DatabaseType;

//generating statistics

public class MainStatistics {

	static String[] querys=new String[]{"APFilesOriginal-MPFiles","APFilesOriginal-MPTrees","APTreesOriginal-MPTrees","MPFilesOriginal-APFiles","MPFilesOriginal-APTrees","MPTreesOriginal-APTrees"};
	static String[] searchMethodsStatisticsOriginal=new String[]{"CosineDistance"/*,"DBSearch"*/,"MassBank","OberacherWithIntensities","OberacherWithIntensitiesAndMasses","TreeAlignment"};
	static String[] searchMethods=new String[]{"CosineDistance"/*,"DBSearch"*/,"MassBank"/*,"OberacherWithIntensities","OberacherWithIntensitiesAndMasses","TreeAlignment"*/};	
	static String[] decoyMethods=new String[]{"Original","RandomPeaks","MixedSpectrum"/*,"ConditionalPeaks"*/,"Reroot",/*"RandomTree",*/"ConditionalFast"/*,"RandomDecoy"*/};

	static double ppm=10;
	static double ae=2;

	//	static String msBlastFolder="U:\\MSBlast\\";
	static String msBlastFolder="D:\\metabo_tandem_ms\\MSBlast\\";
	static String addFolder="";

	static String ser;
	static File[] originalDBFolders;

	static boolean getHitLists=false;
	static boolean statisticsOriginal=false;
	static boolean getRankStatistic=false;
	static boolean checkEquality=false;
	static boolean getQValues=false;
	static boolean getQValuesAverage=false;
	static boolean getQValueSummary=false;
	static boolean getQValueSummaryPerQValue=false;
	static boolean getQValueSummaryPerQValueEstimated=false;
	static boolean getQValueSummaryPerPosition=false;

	static final String sep=System.getProperty("file.separator");

	public static void init(String[] args){

		int i=0;
		if(args!=null){
			while(i<args.length){
				if(args[i].equals("-base"))msBlastFolder=args[++i];			
				if(args[i].equals("-addFolder"))addFolder=args[++i];
				if(args[i].equals("-ppm"))ppm=Double.parseDouble(args[++i]);
				if(args[i].equals("-ae"))ae=Double.parseDouble(args[++i]);
				if(args[i].equals("-queries"))querys=args[++i].split(",");
				if(args[i].equals("-searchMethods"))searchMethods=args[++i].split(",");
				if(args[i].equals("-searchMethodsStatisticsOriginal"))searchMethodsStatisticsOriginal=args[++i].split(",");
				if(args[i].equals("-decoyMethods"))decoyMethods=args[++i].split(",");

				if(args[i].equals("-getHitLists"))getHitLists=true;
				if(args[i].equals("-statisticsOriginal"))statisticsOriginal=true;
				if(args[i].equals("-getRankStatistic"))getRankStatistic=true;
				if(args[i].equals("-checkEquality"))checkEquality=true;
				if(args[i].equals("-getQValues"))getQValues=true;
				if(args[i].equals("-getQValuesAverage"))getQValuesAverage=true;
				if(args[i].equals("-getQValueSummary"))getQValueSummary=true;
				if(args[i].equals("-getQValueSummaryPerQValue"))getQValueSummaryPerQValue=true;
				if(args[i].equals("-getQValueSummaryPerQValueEstimated"))getQValueSummaryPerQValueEstimated=true;
				if(args[i].equals("-getQValueSummaryPerPosition"))getQValueSummaryPerPosition=true;
				i++;
			}
		} 

		ser=msBlastFolder+"Data"+sep+"originalDB.ser";
		originalDBFolders=new File[]{
				new File(msBlastFolder+"Data"+sep+"DataBases"+sep+"Agilent"+sep+"pos"+sep+"AgilentPositiveTreesOriginal"+sep+"massbank"),
				new File(msBlastFolder+"Data"+sep+"DataBases"+sep+"Agilent"+sep+"pos"+sep+"AgilentPositiveFilesOriginal"+sep+"massbank"),
				new File(msBlastFolder+"Data"+sep+"DataBases"+sep+"Massbank"+sep+"pos"+sep+"MassbankPositiveTreesOriginal"+sep+"massbank"),
				new File(msBlastFolder+"Data"+sep+"DataBases"+sep+"Massbank"+sep+"pos"+sep+"MassbankPositiveFilesOriginal"+sep+"massbank"),
				new File(msBlastFolder+"Data"+sep+"DataBases"+sep+"MassbankOrbi"+sep+"pos"+sep+"MassbankOrbiPositiveTreesOriginal"+sep+"massbank"),
				new File(msBlastFolder+"Data"+sep+"DataBases"+sep+"MassbankOrbi"+sep+"pos"+sep+"MassbankOrbiPositiveFilesOriginal"+sep+"massbank"),
				new File(msBlastFolder+"Data"+sep+"DataBases"+sep+"MassbankQTof"+sep+"pos"+sep+"MassbankQTofPositiveTreesOriginal"+sep+"massbank"),
				new File(msBlastFolder+"Data"+sep+"DataBases"+sep+"MassbankQTof"+sep+"pos"+sep+"MassbankQTofPositiveFilesOriginal"+sep+"massbank"),
				new File(msBlastFolder+"Data"+sep+"DataBases"+sep+"Gnps"+sep+"pos"+sep+"GnpsPositiveTreesOriginal"+sep+"massbank"),
				new File(msBlastFolder+"Data"+sep+"DataBases"+sep+"Gnps"+sep+"pos"+sep+"GnpsPositiveFilesOriginal"+sep+"massbank")	
		};
	}

	public static void main(String args[]) throws Exception{
		init(args);

		File searchResultsFile=new File(msBlastFolder+"searchResults"+addFolder);
		if(!searchResultsFile.exists())searchResultsFile.mkdirs();
		File statisticsFile=new File(msBlastFolder+"statistics"+addFolder+sep+"pos");
		if(!statisticsFile.exists())statisticsFile.mkdirs();
		File graphicsFile=new File(msBlastFolder+sep+"graphics_Statistics"+addFolder+sep+"pos");
		if(!graphicsFile.exists())graphicsFile.mkdirs();

		if(getHitLists)getHitLists(searchResultsFile,statisticsFile);
		if(statisticsOriginal){
			statisticsOriginal(statisticsFile, new File(statisticsFile.getAbsolutePath()+sep+"statOriginal.txt"), new File (statisticsFile.getAbsolutePath()+sep+"changesOriginal.txt"), graphicsFile);
			statisticsHighFDR(statisticsFile, new File(statisticsFile.getAbsolutePath()+sep+"statOriginalDMEHighFDR.txt"), new File (statisticsFile.getAbsolutePath()+sep+"changesOriginalDMEHighFDR.txt"), graphicsFile);
		}
		if(getRankStatistic)getRankStatistic(searchResultsFile,graphicsFile);

		if(checkEquality)checkEquality(searchResultsFile,graphicsFile);

		for(int i:new int[]{1,5}){
			String add="";
			if(i==1)add="_AllCompounds";
			if(i==5)add="_AllSelections";

			File qValuesFileAdd=new File(msBlastFolder+"qValues"+add+addFolder+sep+"pos");
			if(!qValuesFileAdd.exists())qValuesFileAdd.mkdirs();
			File graphicsFileAdd=new File(msBlastFolder+"graphics"+add+addFolder+sep+"pos");
			if(!graphicsFileAdd.exists())graphicsFileAdd.mkdirs();

			if(getQValues)getQValues(i, statisticsFile, qValuesFileAdd);
			if(getQValuesAverage)getQValuesAverage(qValuesFileAdd);

			if(getQValueSummary)getQValueSummary(new double[]{/*0,0,*/0,0},new double[]{/*0.01,0.05,*/0.1,1},qValuesFileAdd,graphicsFileAdd);
			if(getQValueSummaryPerQValue)getQValueSummaryPerQValue(new double[]{/*0,0,*/0,0},new double[]{/*0.01,0.05,*/0.1,1}, qValuesFileAdd, graphicsFileAdd, graphicsFileAdd, true);
			if(getQValueSummaryPerQValueEstimated)getQValueSummaryPerQValue(new double[]{/*0,0,*/0,0},new double[]{/*0.01,0.05,*/0.1,1}, qValuesFileAdd, graphicsFileAdd, graphicsFileAdd, false);
			if(getQValueSummaryPerPosition)getQValueSummaryPerPosition(new double[]{/*0,0,*/0,0},new double[]{/*0.1,0.3,*/0.5,1},qValuesFileAdd, graphicsFileAdd, graphicsFileAdd);
		}

		System.out.println("FINISHED");
	}

	private static void getQValuesAverage(File outputFileQValues) throws Exception{
		Locale.setDefault(Locale.US);
		DecimalFormat df=new DecimalFormat("0.00");		
		for(String masses:new String[]{"HitStatistic_DifferentMassesExcluded_","HitStatistic_DMEHighFDR_"}){
			//				BufferedWriter bw=new BufferedWriter(new FileWriter(new File(outputFileQValues.getAbsolutePath()+sep+masses+"qValueSummary_"+df.format(minQValue)+"-"+df.format(maxQValue)+".txt")));

			List<String> dataset=new ArrayList<String>();
			List<String> lines=new ArrayList<String>();
			List<String> rows=new ArrayList<String>();

			//			searchMethods=new String[]{"CosineDistance"};
			//			querys=new String[]{"APFilesOriginal-MPFiles"};
			//			decoyMethods=new String[]{"RandomPeaks"};


			for(String query:querys)dataset.add(query);
			for(String searchMethod:searchMethods)lines.add(searchMethod);
			for(String s:decoyMethods)if(!s.equals("Original"))rows.add(s);


			for (int i=0;i<querys.length;i++){
				String queryLong=querys[i];

				for(String searchMethod:searchMethods){

					for(String s:decoyMethods){
						if(s.equals("Original"))continue;

						File qValueFolderSelection=new File(outputFileQValues.getAbsolutePath()+sep+searchMethod+sep+masses+queryLong+s);
						if(!qValueFolderSelection.exists())continue;	
						System.out.println(qValueFolderSelection.getAbsolutePath());
						getQValuesAverage(Arrays.asList(qValueFolderSelection.listFiles(new FileExtensionFilter("qValueMean"))),"Mean");
						getQValuesAverage(Arrays.asList(qValueFolderSelection.listFiles(new FileExtensionFilter("qValueMergeMean"))),"MergeMean");
					}														
				}
			}

		}
	}

	private static void getQValueSummary(double[] minQValues, double[] maxQValues, File outputFileQValues, File outputFileGraphics) throws Exception{
		Locale.setDefault(Locale.US);
		DecimalFormat df=new DecimalFormat("0.00");
		for(int q=0;q<minQValues.length;q++){
			double minQValue=minQValues[q];
			double maxQValue=maxQValues[q];
			for(boolean normalizedEstimated:new boolean[]{true, false}){
				for(String extension:new String[]{"qValue","qValueMean","qValueMerge","qValueMergeMean","qValueMeanAverage","qValueMergeMeanAverage"}){
					if(extension.contains("Merge")&&normalizedEstimated)continue;
					for(String masses:new String[]{"HitStatistic_DifferentMassesExcluded_","HitStatistic_DMEHighFDR_"}){
						//				BufferedWriter bw=new BufferedWriter(new FileWriter(new File(outputFileQValues.getAbsolutePath()+sep+masses+"qValueSummary_"+df.format(minQValue)+"-"+df.format(maxQValue)+".txt")));

						List<String> dataset=new ArrayList<String>();
						List<String> lines=new ArrayList<String>();
						List<String> rows=new ArrayList<String>();

						//			searchMethods=new String[]{"CosineDistance"};
						//			querys=new String[]{"APFilesOriginal-MPFiles"};
						//			decoyMethods=new String[]{"RandomPeaks"};


						for(String query:querys)dataset.add(query);
						for(String searchMethod:searchMethods)lines.add(searchMethod);
						for(String s:decoyMethods)if(!s.equals("Original"))rows.add(s);


						for (int i=0;i<querys.length;i++){
							String queryLong=querys[i];

							for(String searchMethod:searchMethods){

								XYSeriesCollection ds = new XYSeriesCollection();
								XYSeries seriesWH=SimilarityMatrix.getBisectingLine(minQValue, maxQValue);
								ds.addSeries(seriesWH);

								for(String s:decoyMethods){
									if(s.equals("Original"))continue;

									File qValueFolderSelection=new File(outputFileQValues.getAbsolutePath()+sep+searchMethod+sep+masses+queryLong+s);
									if(!qValueFolderSelection.exists())continue;	
									System.out.println(qValueFolderSelection.getName());
									double[] v=getAverageDistance(Arrays.asList(qValueFolderSelection.listFiles(new FileExtensionFilter(extension))), minQValue, maxQValue, normalizedEstimated);
									//							results[dataset.indexOf(queryLong)][lines.indexOf(searchMethod)][rows.indexOf(s)]=v;

									File f=new File(outputFileGraphics.getAbsolutePath()+sep+searchMethod+sep+"qValues"+(normalizedEstimated?"_norm":"")+sep+extension+sep+masses+queryLong+s+"_"+df.format(minQValue)+"-"+df.format(maxQValue)+".jpg");
									if(!f.getParentFile().exists())f.getParentFile().mkdirs();
									List<double[]> averageQValues=writeQValues(f, Arrays.asList(qValueFolderSelection.listFiles(new FileExtensionFilter(extension))), minQValue, maxQValue, normalizedEstimated);

									//									File txt =new File(outputFileGraphics.getAbsoluteFile()+sep+searchMethod+sep+"qValues"+(normalizedEstimated?"_norm":"")+sep+extension+sep+masses+"Average_"+s+"_"+queryLong+"_"+df.format(minQValue)+"-"+df.format(maxQValue)+".txt");
									//									BufferedWriter bw=new BufferedWriter(new FileWriter(txt)); 
									//									bw.write("real qValue\testimated qValue ("+s+")");
									//									bw.newLine();

									XYSeries xyseries=new XYSeries(s);
									for(double[] e:averageQValues){
										xyseries.add(e[0], e[1]);
										//										bw.write(e[0]+"\t"+e[1]);
										//										bw.newLine();
									}
									ds.addSeries(xyseries);
									//									bw.close();

									System.out.println("getQValueSummary: " + (normalizedEstimated?"_norm":"") + " "+searchMethod+" "+queryLong+s+" "+v[0]+" "+v[1]+" "+v[2]);
								}							

								if(ds.getSeries().size()>1){
									File jpg =new File(outputFileGraphics.getAbsoluteFile()+sep+searchMethod+sep+"qValues"+(normalizedEstimated?"_norm":"")+sep+extension+sep+masses+"Joint_"+queryLong+"_"+df.format(minQValue)+"-"+df.format(maxQValue)+".jpg");
									if(!jpg.getParentFile().exists())jpg.getParentFile().mkdirs();
									final JFreeChart chart =ChartFactory.createScatterPlot("averageQValues", "qValue calculated", "qValue by decoy database",ds);
									ChartUtilities.saveChartAsJPEG(jpg, chart, 1000, 1000);
								}
								System.out.println();
							}
						}
					}
				}
			}
		}

	}

	private static void getQValueSummaryPerQValue(double[] minQValues, double[] maxQValues, File inputFileQValues, File outputFileQValues, File outputFileGraphics, boolean useRealQValue) throws Exception{
		Locale.setDefault(Locale.US);

		Stroke dashedStroke = new BasicStroke(
				1.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND,
				1.0f, new float[] {10.0f, 6.0f}, 0.0f );

		DecimalFormat df=new DecimalFormat("0.00");
		for(int q=0;q<minQValues.length;q++){
			double minQValue=minQValues[q];
			double maxQValue=maxQValues[q];	
			for(boolean normalizedEstimated:new boolean[]{false, true}){
				for(String extension:new String[]{"qValue","qValueMerge","qValueAverage"}){
					if(extension.equals("qValueMerge")&&normalizedEstimated)continue;
					for(String masses:new String[]{"HitStatistic_DifferentMassesExcluded_","HitStatistic_DMEHighFDR_"}){			

						List<String> dataset=new ArrayList<String>();
						List<String> lines=new ArrayList<String>();
						List<String> rows=new ArrayList<String>();

						//			searchMethods=new String[]{"CosineDistance"};
						//			querys=new String[]{"APFilesOriginal-MPFiles"};
						//			decoyMethods=new String[]{"RandomPeaks"};


						for(String query:querys)dataset.add(query);
						for(String searchMethod:searchMethods)lines.add(searchMethod);
						for(String s:decoyMethods)if(!s.equals("Original"))rows.add(s);

						Map<String,TreeMap<Double,double[]>> results[][][]=new Map[dataset.size()][lines.size()][rows.size()];


						for (int k=0;k<dataset.size();k++){
							String queryLong=dataset.get(k);

							for(int i=0;i<lines.size();i++){
								String searchMethod=lines.get(i);

								Set<String> outputs=new HashSet<String>();
								for(int j=0;j<rows.size();j++){
									String s=rows.get(j);
									if(s.equals("Original"))continue;

									File qValueFolderSelection=new File(inputFileQValues.getAbsolutePath()+sep+searchMethod+sep+masses+queryLong+s);
									if(!qValueFolderSelection.exists())continue;						
									Map<String,TreeMap<Double,double[]>> v=getAverageDistancePerQValue(minQValue, maxQValue, useRealQValue, Arrays.asList(qValueFolderSelection.listFiles(new FileExtensionFilter(extension))), normalizedEstimated);
									results[k][i][j]=v;
									outputs.addAll(v.keySet());

									System.out.println("getQValueSummaryPerQValue: " +(normalizedEstimated?"_norm":"")+useRealQValue+" "+ searchMethod+" "+queryLong+s);
								}

								for(String o:outputs){

									File txtFile=new File(outputFileQValues.getAbsolutePath()+sep+lines.get(i)+sep+"qValuePerQValue"+(useRealQValue?"":"Estimated")+(normalizedEstimated?"_norm":"")+sep+extension+sep+masses+"qValuePerQValueSummary_"+dataset.get(k)+"_"+lines.get(i)+"_"+df.format(minQValue)+"-"+df.format(maxQValue)+o+".txt");
									if(!txtFile.getParentFile().exists())txtFile.getParentFile().mkdirs();
									BufferedWriter bw=new BufferedWriter(new FileWriter(txtFile));
									for(int j=0;j<rows.size();j++){
										if(results[k][i][j]!=null&&results[k][i][j]!=null)bw.write("\t"+lines.get(i)+"_"+rows.get(j)+" average\t"+lines.get(i)+"_"+rows.get(j)+" standard deviation");
									}
									bw.newLine();

									XYSeriesCollection col = new XYSeriesCollection();
									XYSeriesCollection colAvgDev = new XYSeriesCollection();

									TreeSet<Double> s=new TreeSet<Double>();
									XYSeries[][] series = new XYSeries[lines.size()][rows.size()];
									XYSeries[][][] seriesAvgDev = new XYSeries[lines.size()][rows.size()][3];	
									for(int j=0;j<rows.size();j++){
										if(results[k][i][j]==null||results[k][i][j].get(o)==null)continue;
										s.addAll(results[k][i][j].get(o).keySet());
										series[i][j]=new XYSeries(lines.get(i)+"_"+rows.get(j));
										if(!rows.get(j).equals("RandomDecoy"))col.addSeries(series[i][j]);

										seriesAvgDev[i][j][0]=new XYSeries(lines.get(i)+"_"+rows.get(j)+" standard deviation -");
										seriesAvgDev[i][j][1]=new XYSeries(lines.get(i)+"_"+rows.get(j)+" average");
										seriesAvgDev[i][j][2]=new XYSeries(lines.get(i)+"_"+rows.get(j)+" standard deviation +");
										if(!rows.get(j).equals("RandomDecoy"))colAvgDev.addSeries(seriesAvgDev[i][j][0]);
										if(!rows.get(j).equals("RandomDecoy"))colAvgDev.addSeries(seriesAvgDev[i][j][1]);
										if(!rows.get(j).equals("RandomDecoy"))colAvgDev.addSeries(seriesAvgDev[i][j][2]);
									}


									for(Double key:s){
										bw.write(Double.toString(key));				
										for(int j=0;j<rows.size();j++){								
											if(results[k][i][j]==null||results[k][i][j].get(o)==null)continue;
											bw.write("\t");
											Double[] d =new Double[]{null,null,null};
											Entry<Double, double[]> tmp=results[k][i][j].get(o).floorEntry(key);
											if(tmp!=null){ 
												double v[] =tmp.getValue();
												d=new Double[]{v[0],v[1]};
												bw.write(Double.toString(d[0])+"\t"+Double.toString(d[1]));
												series[i][j].add(key,new Double(d[0]));
												seriesAvgDev[i][j][0].add(key,new Double(d[0]-d[1]));
												seriesAvgDev[i][j][1].add(key,new Double(d[0]));
												seriesAvgDev[i][j][2].add(key,new Double(d[0]+d[1]));
											}else{
												bw.write("null");
											}
										}
										bw.newLine();
									}
									bw.close();

									File jpg =new File(outputFileGraphics.getAbsolutePath()+sep+lines.get(i)+sep+"qValuePerQValue"+(useRealQValue?"":"Estimated")+(normalizedEstimated?"_norm":"")+sep+extension+sep+masses+"qValuePerQValue_"+dataset.get(k)+"_"+lines.get(i)+"_"+df.format(minQValue)+"-"+df.format(maxQValue)+o+".jpg");
									if(!jpg.getParentFile().exists())jpg.getParentFile().mkdirs();
									JFreeChart chart =ChartFactory.createXYLineChart("qValuePerQValue", (useRealQValue?"real":"estimated")+" qValue", "difference", col);								
									ChartUtilities.saveChartAsJPEG(jpg, chart, 1000, 1000);

									jpg =new File(outputFileGraphics.getAbsolutePath()+sep+lines.get(i)+sep+"qValuePerQValue"+(useRealQValue?"":"Estimated")+(normalizedEstimated?"_norm":"")+sep+extension+sep+masses+"qValuePerQValue_"+dataset.get(k)+"_"+lines.get(i)+"_"+df.format(minQValue)+"-"+df.format(maxQValue)+o+"_AvgDev.jpg");
									if(!jpg.getParentFile().exists())jpg.getParentFile().mkdirs();
									chart =ChartFactory.createXYLineChart("qValuePerQValue", (useRealQValue?"real":"estimated")+" qValue", "difference", colAvgDev);

									XYPlot plot = (XYPlot) chart.getPlot(); 
									XYLineAndShapeRenderer renderer = new DeviationRenderer(true, false);
									Color[] c=new Color[]{Color.BLUE, Color.CYAN, Color.RED, Color.GREEN, Color.ORANGE, Color.PINK, Color.YELLOW, Color.MAGENTA};
									for(int j=0;j<rows.size();j++){
										renderer.setSeriesPaint(j*3, c[j%c.length] ); 
										renderer.setSeriesStroke(j*3, dashedStroke);
										renderer.setSeriesPaint(j*3+1, c[j%c.length]);
										renderer.setSeriesStroke(j*3+2, dashedStroke);
										renderer.setSeriesPaint(j*3+2, c[j%c.length]);
									}
									plot.setRenderer(renderer); 

									ChartUtilities.saveChartAsJPEG(jpg, chart, 1000, 1000);	
								}
							}
						}
					}
				}
			}
		}
	}

	private static void getQValueSummaryPerPosition(double[] minQValues, double[] maxQValues, File inputFileQValues, File outputFileQValues, File outputFileGraphics) throws Exception{
		Locale.setDefault(Locale.US);

		Stroke dashedStroke = new BasicStroke(
				1.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND,
				1.0f, new float[] {10.0f, 6.0f}, 0.0f );

		DecimalFormat df=new DecimalFormat("0.00");
		for(int q=0;q<minQValues.length;q++){
			double minQValue=minQValues[q];
			double maxQValue=maxQValues[q];	
			for(boolean normalizedEstimated:new boolean[]{false, true}){
				for(String extension:new String[]{"qValue","qValueMerge","qValueAverage"}){
					if(extension.equals("qValueMerge")&&normalizedEstimated)continue;
					for(String masses:new String[]{"HitStatistic_DifferentMassesExcluded_","HitStatistic_DMEHighFDR_"}){

						List<String> dataset=new ArrayList<String>();
						List<String> lines=new ArrayList<String>();
						List<String> rows=new ArrayList<String>();

						//			searchMethods=new String[]{"CosineDistance"};
						//			querys=new String[]{"APFilesOriginal-MPFiles"};
						//			decoyMethods=new String[]{"RandomPeaks"};


						for(String query:querys)dataset.add(query);
						for(String searchMethod:searchMethods)lines.add(searchMethod);
						for(String s:decoyMethods)if(!s.equals("Original"))rows.add(s);

						Map<String, TreeMap<Double, double[]>> results[][][]=new HashMap[dataset.size()][lines.size()][rows.size()];


						for (int k=0;k<dataset.size();k++){
							String queryLong=dataset.get(k);

							for(int i=0;i<lines.size();i++){
								String searchMethod=lines.get(i);

								Set<String> outputs=new HashSet<String>();

								for(int j=0;j<rows.size();j++){
									String s=rows.get(j);
									if(s.equals("Original"))continue;

									File qValueFolderSelection=new File(inputFileQValues.getAbsoluteFile()+sep+searchMethod+sep+masses+queryLong+s);
									if(!qValueFolderSelection.exists())continue;						
									Map<String,TreeMap<Double,double[]>> v=getAverageDistancePerPosition(minQValue, maxQValue, Arrays.asList(qValueFolderSelection.listFiles(new FileExtensionFilter(extension))), normalizedEstimated);
									results[k][i][j]=v;

									outputs.addAll(v.keySet());

									System.out.println("getQValueSummaryPerPosition: "+(normalizedEstimated?"_norm":"")+ searchMethod+" "+queryLong+s);
								}

								for(String o:outputs){

									File txtFile=new File(outputFileQValues.getAbsolutePath()+sep+lines.get(i)+sep+"qValuePerPosition"+(normalizedEstimated?"_norm":"")+sep+extension+sep+masses+"qValuePerPositionSummary_"+dataset.get(k)+"_"+lines.get(i)+"_"+df.format(minQValue)+"-"+df.format(maxQValue)+o+".txt");
									if(!txtFile.getParentFile().exists())txtFile.getParentFile().mkdirs();
									BufferedWriter bw=new BufferedWriter(new FileWriter(txtFile));
									for(int j=0;j<rows.size();j++){
										if(results[k][i][j]!=null&&results[k][i][j].get(o)!=null)bw.write("\t"+lines.get(i)+"_"+rows.get(j)+" average\t"+lines.get(i)+"_"+rows.get(j)+" standard deviation");
									}
									bw.newLine();

									XYSeriesCollection col = new XYSeriesCollection();
									XYSeriesCollection colAvgDev = new XYSeriesCollection();


									double max=0;
									XYSeries[][] series = new XYSeries[lines.size()][rows.size()];
									XYSeries[][][] seriesAvgDev = new XYSeries[lines.size()][rows.size()][3];
									for(int j=0;j<rows.size();j++){
										if(results[k][i][j]==null||results[k][i][j].get(o)==null)continue;
										max=Math.max(max,results[k][i][j].get(o).size());
										series[i][j]=new XYSeries(lines.get(i)+"_"+rows.get(j));
										if(!rows.get(j).equals("RandomDecoy"))col.addSeries(series[i][j]);

										seriesAvgDev[i][j][0]=new XYSeries(lines.get(i)+"_"+rows.get(j)+" standard deviation -");
										seriesAvgDev[i][j][1]=new XYSeries(lines.get(i)+"_"+rows.get(j)+" average");
										seriesAvgDev[i][j][2]=new XYSeries(lines.get(i)+"_"+rows.get(j)+" standard deviation +");
										if(!rows.get(j).equals("RandomDecoy"))colAvgDev.addSeries(seriesAvgDev[i][j][0]);
										if(!rows.get(j).equals("RandomDecoy"))colAvgDev.addSeries(seriesAvgDev[i][j][1]);
										if(!rows.get(j).equals("RandomDecoy"))colAvgDev.addSeries(seriesAvgDev[i][j][2]);
									}


									for(double s=0;s<max;s++){
										bw.write(Double.toString(s));					
										for(int j=0;j<rows.size();j++){								
											if(results[k][i][j]==null||results[k][i][j].get(o)==null)continue;
											bw.write("\t");
											if(s<results[k][i][j].get(o).size()){
												double[] d=results[k][i][j].get(o).get(s);	
												if(d!=null){
													bw.write(Double.toString(d[0])+"\t"+Double.toString(d[1]));
													series[i][j].add(s,new Double(d[0]));
													seriesAvgDev[i][j][0].add(s,new Double(d[0]-d[1]));
													seriesAvgDev[i][j][1].add(s,new Double(d[0]));
													seriesAvgDev[i][j][2].add(s,new Double(d[0]+d[1]));
												}
											}
										}
										bw.newLine();
									}
									bw.close();

									File jpg =new File(outputFileGraphics.getAbsolutePath()+sep+lines.get(i)+sep+"qValuePerPosition"+(normalizedEstimated?"_norm":"")+sep+extension+sep+masses+"qValuePerPosition_"+dataset.get(k)+"_"+lines.get(i)+"_"+df.format(minQValue)+"-"+df.format(maxQValue)+o+".jpg");
									if(!jpg.getParentFile().exists())jpg.getParentFile().mkdirs();
									JFreeChart chart =ChartFactory.createXYLineChart("qValuePerPosition", "qValue", "difference", col);								
									ChartUtilities.saveChartAsJPEG(jpg, chart, 1000, 1000);

									jpg =new File(outputFileGraphics.getAbsolutePath()+sep+lines.get(i)+sep+"qValuePerPosition"+(normalizedEstimated?"_norm":"")+sep+extension+sep+masses+"qValuePerPosition_"+dataset.get(k)+"_"+lines.get(i)+"_"+df.format(minQValue)+"-"+df.format(maxQValue)+o+"_AvgDev.jpg");
									if(!jpg.getParentFile().exists())jpg.getParentFile().mkdirs();
									chart =ChartFactory.createXYLineChart("qValuePerPosition", "qValue", "difference", colAvgDev);

									XYPlot plot = (XYPlot) chart.getPlot(); 
									XYLineAndShapeRenderer renderer = new DeviationRenderer(true, false);
									Color[] c=new Color[]{Color.BLUE, Color.CYAN, Color.RED, Color.GREEN, Color.ORANGE, Color.PINK, Color.YELLOW, Color.MAGENTA};
									for(int j=0;j<rows.size();j++){
										renderer.setSeriesPaint(j*3, c[j%c.length] ); 
										renderer.setSeriesStroke(j*3, dashedStroke);
										renderer.setSeriesPaint(j*3+1, c[j%c.length]);
										renderer.setSeriesStroke(j*3+2, dashedStroke);
										renderer.setSeriesPaint(j*3+2, c[j%c.length]);
									}
									plot.setRenderer(renderer); 

									ChartUtilities.saveChartAsJPEG(jpg, chart, 1000, 1000);	

								}
							}
						}
					}
				}
			}
		}
	}

	public static List<List<double[]>> getPoints(List<File> inputFiles, double minQValue, double maxQValue, boolean realQValue, boolean decoyQValue, String op, boolean percentage, boolean normalizedEstimated) throws Exception{
		List<List<double[]>> result=new ArrayList<List<double[]>>();
		int all=0;
		for(File f:inputFiles){			
			BufferedReader br=new BufferedReader(new FileReader(f));
			String line=br.readLine();
			if(!line.contains("\t"))line=br.readLine();
			int allTmp=0;
			while((line=br.readLine())!=null)allTmp++;
			all=Math.max(all, allTmp);
			br.close();
		}
		for(File f:inputFiles){
			List<double[]> tmp=new ArrayList<double[]>();
			result.add(tmp);
			BufferedReader br=new BufferedReader(new FileReader(f));
			String line=br.readLine();
			if(!line.contains("\t")||line.split("\t").length>2)line=br.readLine();
			List<String> header=Arrays.asList(line.split("\t"));
			int indexCalc=header.indexOf("calculated qValue");
			int indexEst=header.indexOf("estimated qValue");
			int n=0;
			double maxCalc=0;
			List<double[]> ds=new ArrayList<double[]>();
			while((line=br.readLine())!=null){
				String[] l=line.split("\t");
				double calc=Double.parseDouble(l[indexCalc]);
				double est=Double.parseDouble(l[indexEst]);
				ds.add(new double[]{calc, est});
				maxCalc=Math.max(maxCalc,calc);
			}
			if(normalizedEstimated){
				for(double[] d:ds)d[1]/=maxCalc;
			}
			for(double[] d:ds){
				double calc=d[0];
				double est=d[1];
				boolean add=op.equals("or")?
						((!decoyQValue||(est>=minQValue&&est<=maxQValue))||(!realQValue||(calc>=minQValue&&calc<=maxQValue))):
							((!decoyQValue||(est>=minQValue&&est<=maxQValue))&&(!realQValue||(calc>=minQValue&&calc<=maxQValue)));
						if(percentage){
							double p=1.0*n/all;
							add=(p>=minQValue&&p<=maxQValue);
						}
						if(add){
							tmp.add(d);							
						}
						n++;
			}
			br.close();
		}		
		return result;
	}

	private static List<double[]> writeQValues(File outputFile, List<File> inputFiles, double minQValue, double maxQValue, boolean normalizedEstimated) throws Exception{
		List<List<double[]>> r=getPoints(inputFiles, minQValue, maxQValue, true, false, "and", false, normalizedEstimated);
		Map<Integer,Map<Double, List<Double>>> allValues=new HashMap<Integer,Map<Double, List<Double>>>();
		XYSeriesCollection dataset = new XYSeriesCollection();
		List<double[]> result=new ArrayList<double[]>();
		for(int i=0;i<inputFiles.size();i++){
			XYSeries xyseries=new XYSeries(inputFiles.get(i).getName());
			for(int j=0;j<r.get(i).size();j++){
				double[] d=r.get(i).get(j);
				xyseries.add(d[0],d[1]);
				result.add(new double[]{d[0],d[1]});
				if(!allValues.containsKey(j))allValues.put(j,new TreeMap<Double,List<Double>>());
				if(!allValues.get(j).containsKey(d[0]))allValues.get(j).put(d[0], new ArrayList<Double>());
				allValues.get(j).get(d[0]).add(d[1]);
			}
			dataset.addSeries(xyseries);
		}


		Collections.sort(result,new Comparator<double[]>(){
			@Override
			public int compare(double[] o1, double[] o2) {
				int c=Double.compare(o1[0],o2[0]);
				if(c!=0)return c;
				return Double.compare(o1[1],o2[1]);
			}
		});

		XYSeries seriesWH=SimilarityMatrix.getBisectingLine(minQValue, maxQValue);
		dataset.addSeries(seriesWH);

		final JFreeChart chart = ChartFactory.createScatterPlot("qValue", "qValue calculated", "qValue by decoy database", dataset, PlotOrientation.VERTICAL, false, false, false);
		ChartUtilities.saveChartAsJPEG(outputFile, chart, 1000, 1000);

		return result;
	}

	private static List<List<double[]>> getQValuesAverage(List<File> inputFiles, String add) throws Exception{
		Map<String,List<File>> sortedFiles=new HashMap<String, List<File>>();
		for(File f:inputFiles){
			String name=f.getName();
			name=name.replaceAll("Selection\\d*_","");
			String id=f.getParent()+sep+name.substring(0,name.lastIndexOf('_'))+".qValue"+add+"Average";
			if(!sortedFiles.containsKey(id))sortedFiles.put(id,new ArrayList<File>());
			sortedFiles.get(id).add(f);
		}


		List<List<double[]>> result=new ArrayList<List<double[]>>();
		for(Entry<String,List<File>> f:sortedFiles.entrySet()){			
			List<List<double[]>> r=getPoints(f.getValue(), Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, true, true, "and", false, false);

			Set<Double> keys=new TreeSet<Double>();
			for(List<double[]> dl:r){
				for(double d[]:dl){
					keys.add(d[0]);
				}
			}

			List<double[]> res=new ArrayList<double[]>();
			List<List<Double>> p=new ArrayList<List<Double>>();
			for(double k:keys){
				List<Double> points=new ArrayList<Double>();
				for(List<double[]> dl:r){
					for(int i=1;i<dl.size();i++){
						double k_i=dl.get(i-1)[0];
						double v_i=dl.get(i-1)[1];						
						double k_j=dl.get(i)[0];
						double v_j=dl.get(i)[1];

						if(k_i==k_j)continue;

						if(k_j>=k){
							double m=(v_i-v_j)/(k_i-k_j);
							double n=v_i-m*k_i;
							points.add(m*k+n);
							break;
						}
						if(i==dl.size()-1){
							points.add(v_j);
						}
					}
				}

				//				for(double d:points)System.out.println(d);
				//				System.out.println();

				double[] as=getAverageAndStandardDeviation(points);
				p.add(points);
				res.add(new double[]{k,as[0],as[1]});
			}
			BufferedReader br=new BufferedReader(new FileReader(f.getValue().get(0)));			
			BufferedWriter bw=new BufferedWriter(new FileWriter(f.getKey()));
			bw.write(br.readLine());
			br.close();
			bw.newLine();
			bw.write("calculated qValue\testimated qValue\tstandard deviation\tall values");
			bw.newLine();			
			for(int i=0;i<res.size();i++){
				double[] entry=res.get(i);
				List<Double> entryPoints=p.get(i);
				bw.write(entry[0]+"\t"+entry[1]+"\t"+entry[2]+"\t");
				for(int j=0;j<entryPoints.size();j++){
					double d=entryPoints.get(j);
					bw.write(d+"");
					if(j!=entryPoints.size()-1)bw.write(",");
				}
				bw.newLine();
			}
			bw.close();
			result.add(res);
		}		

		return result;
	}

	private static double[] getAverageDistance(List<File> inputFiles, double minQValue, double maxQValue, boolean normalizedEstimated) throws Exception{
		List<Double> avgs=new ArrayList<Double>();
		List<List<double[]>> r=getPoints(inputFiles, minQValue, maxQValue, true, false, "and", false, normalizedEstimated);
		int points=0;
		for(List<double[]> l:r){
			double sumTmp=0;
			int numberTmp=0;
			for(double[] d:l){
				double calc=d[0];
				double est=d[1];	
				double abw=Math.pow(calc-est,2);				
				sumTmp+=abw;
				numberTmp++;
			}
			double avg=sumTmp/numberTmp;
			points+=numberTmp;
			avgs.add(avg);
		}

		double sum=0;
		for(double a:avgs){
			sum+=a;
		}
		double avg=sum/avgs.size();

		double s=0;
		for(double a:avgs){							
			s+=Math.pow(a-avg,2);
		}
		s=Math.pow(s/(avgs.size()-1), 0.5);
		return new double[]{1.0*points/avgs.size(), avg, s};
	}


	private static Map<String, TreeMap<Double, double[]>> getAverageDistancePerQValue(double minQValue, double maxQValue, boolean useRealQValue, List<File> inputFiles, boolean normalizedEstimated) throws Exception{
		List<List<double[]>> r=null;
		if(useRealQValue)r=getPoints(inputFiles, minQValue, maxQValue, true, false, "and", false, normalizedEstimated);
		if(!useRealQValue){
			r=getPoints(inputFiles, minQValue, maxQValue, false, true, "and", false, normalizedEstimated);
			for(List<double[]> dl:r){
				for(double[] d:dl){
					double tmp=d[0];
					d[0]=d[1];
					d[1]=tmp;
				}
			}
		}
		Map<String, TreeMap<Double, double[]>> result = getResultDataStructure(r);

		for(Double qValueReal:result.get("_square").keySet()){
			Map<String, Map<String, List<Double>>> points = getPointsDataStructure(inputFiles);

			for(int i=0;i<inputFiles.size();i++){
				List<double[]> dl=r.get(i);
				String[] n=inputFiles.get(i).getName().substring(0,inputFiles.get(i).getName().lastIndexOf(".")).split("_");
				double[] point=new double[]{0,0};
				Iterator<double[]> it=dl.iterator();
				while(it.hasNext()){
					double d[]=it.next();
					if(d[0]>qValueReal)
						break;
					point=d;
					it.remove();
				}
				dl.add(0,point);		
				double square=Math.pow(point[1]-point[0],2);
				double abs=Math.abs(point[1]-point[0]);
				double minus=point[1]-point[0];
				points.get("minus").get("").add(minus);
				points.get("square").get("").add(square);
				points.get("abs").get("").add(abs);
				points.get("absPerSelection").get(n[1]).add(abs);
				points.get("squarePerSelection").get(n[1]).add(square);
				points.get("minusPerSelection").get(n[1]).add(minus);
				if(n.length==5){
					points.get("absPerDB").get(n[4]).add(abs);
					points.get("squarePerDB").get(n[4]).add(square);
					points.get("minusPerDB").get(n[4]).add(minus);
				}else{
					points.get("absPerDB").get("").add(abs);
					points.get("squarePerDB").get("").add(square);
					points.get("minusPerDB").get("").add(minus);
				}
			}
			result.get("_square").put(qValueReal, getAverageAndStandarddeviation(points.get("square")));
			result.get("_squarePerDB").put(qValueReal, getAverageAndStandarddeviation(points.get("squarePerDB")));
			result.get("_squarePerSelection").put(qValueReal, getAverageAndStandarddeviation(points.get("squarePerSelection")));
			result.get("_abs").put(qValueReal, getAverageAndStandarddeviation(points.get("abs")));
			result.get("_absPerDB").put(qValueReal, getAverageAndStandarddeviation(points.get("absPerDB")));
			result.get("_absPerSelection").put(qValueReal, getAverageAndStandarddeviation(points.get("absPerSelection")));
			result.get("_minus").put(qValueReal, getAverageAndStandarddeviation(points.get("minus")));
			result.get("_minusPerDB").put(qValueReal, getAverageAndStandarddeviation(points.get("minusPerDB")));
			result.get("_minusPerSelection").put(qValueReal, getAverageAndStandarddeviation(points.get("minusPerSelection")));
		}
		return result;
	}

	private static Map<String, Map<String, List<Double>>> getPointsDataStructure(List<File> inputFiles) {
		Map<String, Map<String,List<Double>>> points=new HashMap<String, Map<String,List<Double>>>();
		for(String s:new String[]{"abs","absPerDB","absPerSelection","square","squarePerDB","squarePerSelection","minus","minusPerDB","minusPerSelection"}){
			points.put(s,new HashMap<String,List<Double>>());
		}
		points.get("abs").put("",new ArrayList<Double>());
		points.get("square").put("",new ArrayList<Double>());
		points.get("minus").put("",new ArrayList<Double>());
		for(File f:inputFiles){
			String[] n=f.getName().substring(0,f.getName().lastIndexOf(".")).split("_");
			points.get("absPerSelection").put(n[1],new ArrayList<Double>());
			points.get("squarePerSelection").put(n[1],new ArrayList<Double>());
			points.get("minusPerSelection").put(n[1],new ArrayList<Double>());
			if(n.length==5){
				points.get("absPerDB").put(n[4],new ArrayList<Double>());
				points.get("squarePerDB").put(n[4],new ArrayList<Double>());
				points.get("minusPerDB").put(n[4],new ArrayList<Double>());
			}else{
				points.get("absPerDB").put("",new ArrayList<Double>());
				points.get("squarePerDB").put("",new ArrayList<Double>());
				points.get("minusPerDB").put("",new ArrayList<Double>());
			}
		}
		return points;
	}

	private static Map<String, TreeMap<Double, double[]>> getResultDataStructure(
			List<List<double[]>> r) {
		Map<String,TreeMap<Double, double[]>> result=new HashMap<String, TreeMap<Double, double[]>>();
		for(String s:new String[]{"_square","_squarePerDB","_squarePerSelection","_abs","_absPerDB","_absPerSelection","_minus","_minusPerDB","_minusPerSelection"}){
			TreeMap<Double, double[]> resultTmp=new TreeMap<Double, double[]>();
			result.put(s, resultTmp);
			for(List<double[]> dl:r){
				for(double[] d:dl){
					resultTmp.put(d[0],new double[2]);
				}
			}
		}
		return result;
	}

	public static double[] getAverageAndStandarddeviation(List<Double> allPoints){
		Map<String, List<Double>> tmp=new HashMap<String, List<Double>>();
		tmp.put("", allPoints);
		return getAverageAndStandarddeviation(tmp);
	}

	public static double[] getAverageAndStandarddeviation(Map<String,List<Double>> allPoints){
		double allAvg=0;
		double allS=0;
		for(List<Double> points:allPoints.values()){
			double sum=0;
			for(double a:points){							
				sum+=a;
			}
			double avg=sum/points.size();
			allAvg+=avg;

			double s=0;
			for(double a:points){							
				s+=Math.pow(a-avg,2);
			}
			s=Math.pow(s/(points.size()), 0.5);
			allS+=s;
		}
		return new double[]{allAvg/allPoints.size(), allS/allPoints.size()};
	}

	public static double[] getAverageAndStandardDeviation(List<Double> points){
		double[] result=new double[2];
		double sum=0;
		for(double a:points){							
			sum+=a;
		}
		double avg=sum/points.size();
		result[0]=avg;

		double s=0;
		for(double a:points){							
			s+=Math.pow(a-avg,2);
		}
		s=Math.pow(s/(points.size()), 0.5);
		result[1]=s;
		return result;
	}

	private static Map<String, TreeMap<Double, double[]>> getAverageDistancePerPosition(double minQValue, double maxQValue, List<File> inputFiles, boolean normalizedEstimated) throws Exception{
		List<List<double[]>> r=getPoints(inputFiles, minQValue, maxQValue, true, true, "and", true, normalizedEstimated);
		int max=0;
		for(List<double[]> dl:r){
			max=Math.max(max,dl.size());		
		}

		Map<String, TreeMap<Double, double[]>> result = getResultDataStructure(r);

		for(int i=0;i<max;i++){

			Map<String, Map<String, List<Double>>> points = getPointsDataStructure(inputFiles);

			for(int k=0;k<inputFiles.size();k++){
				List<double[]> dl=r.get(k);
				String[] n=inputFiles.get(k).getName().substring(0,inputFiles.get(k).getName().lastIndexOf(".")).split("_");
				if(i<dl.size()){
					double[] d=dl.get(i);					
					double square=Math.pow(d[1]-d[0],2);
					double abs=Math.abs(d[1]-d[0]);
					double minus=d[1]-d[0];
					points.get("minus").get("").add(minus);
					points.get("square").get("").add(square);
					points.get("abs").get("").add(abs);
					points.get("absPerSelection").get(n[1]).add(abs);
					points.get("squarePerSelection").get(n[1]).add(square);
					points.get("minusPerSelection").get(n[1]).add(minus);
					if(n.length==5){
						points.get("absPerDB").get(n[4]).add(abs);
						points.get("squarePerDB").get(n[4]).add(square);
						points.get("minusPerDB").get(n[4]).add(minus);
					}else{
						points.get("absPerDB").get("").add(abs);
						points.get("squarePerDB").get("").add(square);
						points.get("minusPerDB").get("").add(minus);
					}
				}
			}

			result.get("_square").put(1.0*i, getAverageAndStandarddeviation(points.get("square")));
			result.get("_squarePerDB").put(1.0*i, getAverageAndStandarddeviation(points.get("squarePerDB")));
			result.get("_squarePerSelection").put(1.0*i, getAverageAndStandarddeviation(points.get("squarePerSelection")));
			result.get("_abs").put(1.0*i, getAverageAndStandarddeviation(points.get("abs")));
			result.get("_absPerDB").put(1.0*i, getAverageAndStandarddeviation(points.get("absPerDB")));
			result.get("_absPerSelection").put(1.0*i, getAverageAndStandarddeviation(points.get("absPerSelection")));
			result.get("_minus").put(1.0*i, getAverageAndStandarddeviation(points.get("minus")));
			result.get("_minusPerDB").put(1.0*i, getAverageAndStandarddeviation(points.get("minusPerDB")));
			result.get("_minusPerSelection").put(1.0*i, getAverageAndStandarddeviation(points.get("minusPerSelection")));
		}
		return result;
	}

	private static void getQValues(int numberSelections, File inputFile, File outputFile) throws Exception{
		TreeMap<String, MassBank> massbankFilesOriginal=Utils.getAllDBFiles(originalDBFolders, ser, false, false);

		int repetitions=1;
		if(numberSelections>1)repetitions=2;

		Map<String,List<List<MassBank>>> selectionFinal=new HashMap<String, List<List<MassBank>>>();
		for(int rep=0;rep<repetitions;rep++){

			Map<String,Set<MassBank>> selectionTmp=new HashMap<String, Set<MassBank>>();
			for(MassBank mb:massbankFilesOriginal.values()){
				String v=mb.massbankID.replaceAll("\\d", "");
				if(!selectionTmp.containsKey(v))selectionTmp.put(v, new HashSet<MassBank>());
				selectionTmp.get(v).add(mb);
			}

			Random r=new Random();
			Map<String,List<List<MassBank>>> selectionTmp2=new HashMap<String, List<List<MassBank>>>();
			for(Entry<String, Set<MassBank>> e:selectionTmp.entrySet()){
				List<List<MassBank>> tmp=new ArrayList<List<MassBank>>();
				selectionTmp2.put(e.getKey(), tmp);
				for(int i=0;i<numberSelections;i++)tmp.add(new ArrayList<MassBank>());
				List<MassBank> tmp2=new ArrayList<MassBank>(e.getValue());
				int i=0;
				while(!tmp2.isEmpty()){
					int next=r.nextInt(tmp2.size());				
					tmp.get(i%numberSelections).add(tmp2.get(next));
					tmp2.remove(next);
					i++;
				}			
			}

			Map<String,List<List<MassBank>>> selection=new HashMap<String, List<List<MassBank>>>();
			for(Entry<String, List<List<MassBank>>> e:selectionTmp2.entrySet()){
				List<List<MassBank>> tmp=new ArrayList<List<MassBank>>();
				selection.put(e.getKey(), tmp);
				for(int i=0;i<e.getValue().size();i++)tmp.add(new ArrayList<MassBank>());
				for(int i=0;i<e.getValue().size();i++){
					for(int j=0;j<e.getValue().size();j++){
						if(i!=j||numberSelections==1)tmp.get(i).addAll(e.getValue().get(j));
					}
				}
			}

			for(Entry<String, List<List<MassBank>>> e:selection.entrySet()){
				if(!selectionFinal.containsKey(e.getKey())){
					selectionFinal.put(e.getKey(), e.getValue());
				}else{
					selectionFinal.get(e.getKey()).addAll(e.getValue());
				}
			}
		}

		for (int i=0;i<querys.length;i++){
			String queryLong=querys[i];

			for(String s:decoyMethods){
				if(s.equals("Original"))continue;

				for(String searchMethod:searchMethods){
					for(String masses:new String[]{"HitStatistic_DifferentMassesExcluded_","HitStatistic_DMEHighFDR_"}){
						File qValuesFile=new File(outputFile.getAbsolutePath()+sep+searchMethod);
						if (!qValuesFile.exists())qValuesFile.mkdirs();

						File tmpFile2=new File(inputFile.getAbsolutePath()+sep+searchMethod+sep+queryLong+s);
						List<File> searchresultsFoldersDecoy = new ArrayList<File>();
						if(tmpFile2.exists()){
							for(File f:tmpFile2.listFiles()){
								if(f.getName().startsWith(masses)&&(numberSelections==1||f.getName().endsWith("_1.txt")))searchresultsFoldersDecoy.add(f);
							}
						}
						File searchresultsFoldersOriginal = new File(inputFile.getAbsolutePath()+sep+searchMethod+sep+masses+queryLong+"Original.txt");
						if(!searchresultsFoldersOriginal.exists())continue;

						if(s.equals("RandomDecoy")){

							HitStatistic hsOriginal=new HitStatistic(searchresultsFoldersOriginal,massbankFilesOriginal);
							String v=hsOriginal.hits.get(0).massbankQuery.massbankID.replaceAll("\\d","");
							for(int k=0;k<selectionFinal.get(v).size();k++){
								List<MassBank> mb=selectionFinal.get(v).get(k);
								File file=new File(outputFile.getAbsolutePath()+sep+searchMethod+sep+masses+queryLong+s+sep+"QValues_Selection"+k+"_"+masses+queryLong+s+".qValue");
								if(!file.getParentFile().exists())file.getParentFile().mkdir();
								HitStatistic.writeQValuesDecoyHits(new File(outputFile.getAbsolutePath()+sep+"RandomDecoyValues"), file, hsOriginal, mb, false);

							}
						}else{
							if(searchresultsFoldersDecoy==null||searchresultsFoldersDecoy.isEmpty())continue;

							for(File f:searchresultsFoldersDecoy){
								//							if(!f.getName().equals("HitStatistic_MPTreesOriginal-APTreesRandomPeaks_1.txt")||!searchMethod.equals("MassBank"))continue;
								System.out.print("getQValues: "+f.getAbsolutePath()+"...");
								HitStatistic hsDecoy=new HitStatistic(f,massbankFilesOriginal);
								HitStatistic hsOriginal=new HitStatistic(searchresultsFoldersOriginal,massbankFilesOriginal);

								//							HitStatistic.writeEstimatedQValueVSCalculatedQValueToFile(new File(statisticsFile.getAbsolutePath()+sep+"QValues_"+f.getName()), hsOriginal, hsDecoy, HitStatistic.HitMerge.Add, false, null);
								String v=hsOriginal.hits.get(0).massbankQuery.massbankID.replaceAll("\\d","");
								for(int k=0;k<selectionFinal.get(v).size();k++){
									List<MassBank> mb=selectionFinal.get(v).get(k);																		
									File folder=new File(qValuesFile.getAbsolutePath()+sep+f.getName().substring(0,f.getName().lastIndexOf("_")));
									if(!folder.exists())folder.mkdir();
									HitStatistic.writeEstimatedQValueVSCalculatedQValueToFile(new File(folder.getAbsolutePath()+sep+"QValues_Selection"+k+"_"+f.getName().replaceAll(".txt",".qValue")), new File(folder.getAbsolutePath()+sep+"QValues_Selection"+k+"_"+f.getName().replaceAll(".txt",".qValueMean")), new File(folder.getAbsolutePath()+sep+"QValues_Selection"+k+"_"+f.getName().replaceAll(".txt", ".hitlist")), hsOriginal, hsDecoy, HitStatistic.HitMerge.Add, false, mb);
									HitStatistic.writeEstimatedQValueVSCalculatedQValueToFile(new File(folder.getAbsolutePath()+sep+"QValues_Selection"+k+"_"+f.getName().replaceAll(".txt",".qValueMerge")), new File(folder.getAbsolutePath()+sep+"QValues_Selection"+k+"_"+f.getName().replaceAll(".txt",".qValueMergeMean")), new File(folder.getAbsolutePath()+sep+"QValues_Selection"+k+"_"+f.getName().replaceAll(".txt", ".hitlistMerge")), hsOriginal, hsDecoy, HitStatistic.HitMerge.Merge, false, mb);
									//System.out.println("Selection "+k+" done.");							
								}

								System.out.println("done");

							}
						}
					}
				}
			}
		}

	}

	private static void statisticsOriginal(File inputFile, File outputFile, File outputFileChanges, File graphicsFolder) throws IOException, Exception {
		BufferedWriter bw=new BufferedWriter(new FileWriter(outputFile));
		bw.write("\t");
		bw.write("\t");
		bw.write(HitStatistic.getHeaderString());
		bw.write("\t\t");
		bw.write(HitStatistic.getHeaderString());
		bw.write("\tat least one hit");
		bw.newLine();

		TreeMap<String, MassBank> massbankFilesOriginal=Utils.getAllDBFiles(originalDBFolders, ser, false, false);

		TreeMap<String, Map<MassBank, HitDescription>> stats=new TreeMap<String, Map<MassBank, HitDescription>>();

		for (int i=0;i<querys.length;i++){
			TreeMap<String, Map<MassBank, HitDescription>> statsTmp=new TreeMap<String, Map<MassBank, HitDescription>>();
			TreeMap<String, HitStatistic> allHitstatistics=new TreeMap<String, HitStatistic>();
			TreeMap<String, HitStatistic> allHitstatisticsDifferentMassesExcluded=new TreeMap<String, HitStatistic>();

			String queryLong=querys[i];
			String decoy="Original";

			System.out.println("statisticsOriginal: "+queryLong+decoy);
			bw.write(queryLong+decoy);bw.newLine();			

			for(String searchMethod:searchMethodsStatisticsOriginal){				

				File hitstatistics = new File(inputFile.getAbsolutePath()+sep+searchMethod+sep+"HitStatistic_"+queryLong+decoy+".txt");
				File hitstatisticsDifferentMassesExcluded = new File(inputFile.getAbsolutePath()+sep+searchMethod+sep+"HitStatistic_DifferentMassesExcluded_"+queryLong+decoy+".txt");

				if(hitstatistics.exists()&&hitstatisticsDifferentMassesExcluded.exists()){
					HitStatistic hits=new HitStatistic(hitstatistics, massbankFilesOriginal);					
					bw.write(searchMethod+"\tall\t");
					bw.write(hits.toEntryString());
					Map<MassBank, HitDescription> stat=hits.getHitDescription();
					stats.put(queryLong+decoy+"-"+searchMethod, stat);
					allHitstatistics.put(queryLong+decoy+"-"+searchMethod, hits);

					hits=new HitStatistic(hitstatisticsDifferentMassesExcluded, massbankFilesOriginal);					
					bw.write("\tdifferent masses excluded:\t");
					bw.write(hits.toEntryString());
					bw.write("\t"+hits.getNumberEntriesWithAtLeastOneHit());
					stat=hits.getHitDescription();
					statsTmp.put(queryLong+decoy+"-"+searchMethod+"-SameMass", stat);
					allHitstatisticsDifferentMassesExcluded.put(queryLong+decoy+"-"+searchMethod+"-SameMass", hits);

					bw.newLine();
				}
				bw.flush();
			}
			HitStatistic.writeROCCurves(new File(graphicsFolder.getAbsolutePath()+sep+"ROCCurve_"+queryLong+decoy+".jpg"), allHitstatistics);
			HitStatistic.writeROCCurves(new File(graphicsFolder.getAbsolutePath()+sep+"ROCCurve_DifferentMassesExcluded_"+queryLong+decoy+".jpg"), allHitstatisticsDifferentMassesExcluded);
			stats.putAll(statsTmp);
		}
		bw.close();

		compareSearchMethods(outputFileChanges, stats);
	}
	
	private static void statisticsHighFDR(File inputFile, File outputFile, File outputFileChanges, File graphicsFolder) throws IOException, Exception {
		BufferedWriter bw=new BufferedWriter(new FileWriter(outputFile));
		bw.write("\t\t");
		bw.write(HitStatistic.getHeaderString());
		bw.write("\tat least one hit");
		bw.newLine();

		TreeMap<String, MassBank> massbankFilesOriginal=Utils.getAllDBFiles(originalDBFolders, ser, false, false);

		TreeMap<String, Map<MassBank, HitDescription>> stats=new TreeMap<String, Map<MassBank, HitDescription>>();

		for (int i=0;i<querys.length;i++){
			TreeMap<String, Map<MassBank, HitDescription>> statsTmp=new TreeMap<String, Map<MassBank, HitDescription>>();			
			TreeMap<String, HitStatistic> allHitstatisticsDifferentMassesExcluded=new TreeMap<String, HitStatistic>();

			String queryLong=querys[i];
			String decoy="Original";

			System.out.println("statisticsOriginal: "+queryLong+decoy);
			bw.write(queryLong+decoy);bw.newLine();			

			for(String searchMethod:searchMethodsStatisticsOriginal){				

				File hitstatisticsDifferentMassesExcluded = new File(inputFile.getAbsolutePath()+sep+searchMethod+sep+"HitStatistic_DMEHighFDR_"+queryLong+decoy+".txt");

				if(hitstatisticsDifferentMassesExcluded.exists()){
					HitStatistic hits=new HitStatistic(hitstatisticsDifferentMassesExcluded, massbankFilesOriginal);					
					bw.write("searchMethod\tdifferent masses excluded:\t");
					bw.write(hits.toEntryString());
					bw.write("\t"+hits.getNumberEntriesWithAtLeastOneHit());
					Map<MassBank, HitDescription> stat=hits.getHitDescription();
					statsTmp.put(queryLong+decoy+"-"+searchMethod+"-SameMass", stat);
					allHitstatisticsDifferentMassesExcluded.put(queryLong+decoy+"-"+searchMethod+"-SameMass", hits);

					bw.newLine();
				}
				bw.flush();
			}
			HitStatistic.writeROCCurves(new File(graphicsFolder.getAbsolutePath()+sep+"ROCCurve_DMEHighFDR_"+queryLong+decoy+".jpg"), allHitstatisticsDifferentMassesExcluded);
			stats.putAll(statsTmp);
		}
		bw.close();

		compareSearchMethods(outputFileChanges, stats);
	}

	private static void compareSearchMethods(File outputFileChanges, TreeMap<String, Map<MassBank, HitDescription>> stats) throws IOException {		

		TreeMap<String, Map<String, HitDescription>> tmp=new TreeMap<String, Map<String, HitDescription>>();
		for(Entry<String, Map<MassBank, HitDescription>> e:stats.entrySet()){
			if(!tmp.containsKey(e.getKey()))tmp.put(e.getKey(), new TreeMap<String, HitDescription>());
			Map<String, HitDescription> tmp2=tmp.get(e.getKey());
			for(Entry<MassBank, HitDescription> e2:e.getValue().entrySet()){
				tmp2.put(SimilarityMatrix.getIDDatasetCharge(e2.getKey().massbankID), e2.getValue());
			}
		}

		Map<String, int[][]> countsDetails=new TreeMap<String,int[][]>();
		List<String> allQueries=new ArrayList<String>();
		allQueries.addAll(tmp.keySet());

		for(int i=0;i<allQueries.size();i++){			
			for(int j=0;j<allQueries.size();j++){
				Map<String, HitDescription> stat1=tmp.get(allQueries.get(i));
				Map<String, HitDescription> stat2=tmp.get(allQueries.get(j));
				Set<String> allIDs=new HashSet<String>();
				allIDs.addAll(stat1.keySet());
				allIDs.addAll(stat2.keySet());

				for(String id:allIDs){				
					HitDescription hd1=stat1.get(id);
					HitDescription hd2=stat2.get(id);					
					if(hd1==null)hd1=HitDescription.NotInQuery;			
					if(hd2==null)hd2=HitDescription.NotInQuery;
					Integer c=HitStatistic.compareHitDescription(hd1, hd2);
					String countsQualityKey="_x";
					if(c!=null){
						if(c==0)countsQualityKey="_equal";
						else if(c<0)countsQualityKey="_better";
						else if(c>0)countsQualityKey="_worse";
					}
					if(!countsDetails.containsKey(countsQualityKey))countsDetails.put(countsQualityKey, new int[allQueries.size()][allQueries.size()]);
					countsDetails.get(countsQualityKey)[i][j]++;

					String s1=HitStatistic.getHitDescriptionString(hd1);
					String s2=HitStatistic.getHitDescriptionString(hd2);
					String countsDetailsKey=s1+"->"+s2;
					if(!countsDetails.containsKey(countsDetailsKey))countsDetails.put(countsDetailsKey, new int[allQueries.size()][allQueries.size()]);
					countsDetails.get(countsDetailsKey)[i][j]++;
				}
			}
		}
		BufferedWriter bw2=new BufferedWriter(new FileWriter(outputFileChanges));

		List<String> keys=new ArrayList<String>();
		keys.addAll(countsDetails.keySet());
		for (String k:keys)bw2.write("\t"+k);
		bw2.newLine();
		for(int i=0;i<allQueries.size();i++){				
			for(int j=0;j<allQueries.size();j++){
				String[] query1=allQueries.get(i).split("-");
				String[] query2=allQueries.get(j).split("-");
				int diff=0;
				for(int z=0;z<Math.min(query1.length, query2.length);z++){
					if(!query1[z].equals(query2[z]))diff++;
				}
				diff+=Math.abs(query1.length-query2.length);

				if(diff!=1)continue;
				bw2.write(allQueries.get(i)+" vs. "+allQueries.get(j));
				for(String k:keys){						
					int[][] c=countsDetails.get(k);
					bw2.write("\t"+c[i][j]);
				}
				bw2.newLine();
			}			
		}
		bw2.close();
	}


	private static void getHitLists(File searchResultsFolder, File statisticsFolder) throws IOException, Exception {

		TreeMap<String, MassBank> massbankFilesOriginal=Utils.getAllDBFiles(originalDBFolders, ser, false, false);

		Map<String, Map<MassBank, List<MassBank>>> exclusions70= new HashMap<String, Map<MassBank, List<MassBank>>>();

		for (int i=0;i<querys.length;i++){
			String queryLong=querys[i];


			for(String s:decoyMethods){

				for(String searchMethod:searchMethodsStatisticsOriginal){
					File statisticsFile=new File(statisticsFolder.getAbsolutePath()+sep+searchMethod);
					if (!statisticsFile.exists())statisticsFile.mkdir();

					String ending=(searchMethod.equals("DBSearch")?".txt":".csv");										
					File files=new File(searchResultsFolder.getAbsolutePath()+sep+searchMethod+sep+"pos"+sep+queryLong+s+ending);
					boolean isDirectory=!files.exists();
					List<File> searchresultsFoldersDecoy = getSearchResultsFiles(files);
					if(searchresultsFoldersDecoy.isEmpty())continue;

					for(File f:searchresultsFoldersDecoy){
						SimilarityMatrix mDecoy=new SimilarityMatrix(f, massbankFilesOriginal);
						String name=f.getName();
						if(name.endsWith(".csv"))name=name.replaceAll(".csv",".txt");
						File tmp=new File(statisticsFile.getAbsolutePath()+sep+(isDirectory?f.getParentFile().getName():""));
						if(!tmp.exists())tmp.mkdir();

						System.out.println("getHitLists: " + tmp.getAbsolutePath()+sep+"HitStatistic_"+name);
						
						HitStatistic hits=null;

						mDecoy.clearExclusionList();						
						hits=mDecoy.getHitStatistic();						
						hits.writeHitStatisticToFile(new File(tmp.getAbsolutePath()+sep+"HitStatistic_"+name));


						String[] split=queryLong.split("-");
						boolean isTrees=false;
						if(split[0].contains("Trees")&&split[1].contains("Trees"))isTrees=true;
						mDecoy.clearExclusionList();
						if(isTrees){
							mDecoy.excludeDifferentFormulas();
						}else{
							mDecoy.excludeDifferentMass(ppm, ae);
						}
						hits=mDecoy.getHitStatistic();
						hits.writeHitStatisticToFile(new File(tmp.getAbsolutePath()+sep+"HitStatistic_DifferentMassesExcluded_"+name));

						String data=queryLong;
						if(!exclusions70.containsKey(data))exclusions70.put(data, getExclusions(0.5,mDecoy));
						mDecoy.excludeMassBanksFromDB(exclusions70.get(data));
						hits=mDecoy.getHitStatistic();
						hits.writeHitStatisticToFile(new File(tmp.getAbsolutePath()+sep+"HitStatistic_DMEHighFDR_"+name));

					}
				}

			}
		}
	}


	private static Map<MassBank, List<MassBank>> getExclusions(double percentage, SimilarityMatrix mDecoy) {
		System.out.println("getExclusionList");
		System.out.println("methods DB");
		for(String s:mDecoy.methodsDB){
//			System.out.print(s);
		}
		System.out.println();
		System.out.println("methods query");
		for(String s:mDecoy.methodsQuery){
//			System.out.print(s);
		}
		System.out.println();
		Map<MassBank, List<MassBank>> result=new HashMap<MassBank, List<MassBank>>();
		Map<MassBank, List<MassBank>> hitsTP=new HashMap<MassBank, List<MassBank>>();
		Map<MassBank, List<MassBank>> hitsFP=new HashMap<MassBank, List<MassBank>>();

		for(int i=0;i<mDecoy.massbankQuery.size();i++){
			hitsTP.put(mDecoy.massbankQuery.get(i),new ArrayList<MassBank>());
			hitsFP.put(mDecoy.massbankQuery.get(i),new ArrayList<MassBank>());
			for(int j=0;j<mDecoy.massbankDB.size();j++){
				boolean isEntry=!Double.isNaN(mDecoy.getSimilarityValue(i, j));
				boolean isSameInchi=mDecoy.massbankQuery.get(i).hasEqualInChiKey(mDecoy.massbankDB.get(j));
				if(isEntry){
					if(isSameInchi)hitsTP.get(mDecoy.massbankQuery.get(i)).add(mDecoy.massbankDB.get(j));
					else hitsFP.get(mDecoy.massbankQuery.get(i)).add(mDecoy.massbankDB.get(j));
				}
			}
		}

		List<MassBank> onlyTP=new ArrayList<MassBank>();
		List<MassBank> onlyFP=new ArrayList<MassBank>();
		List<MassBank> both=new ArrayList<MassBank>();
		List<MassBank> none=new ArrayList<MassBank>();
		for(int i=0;i<mDecoy.massbankQuery.size();i++){
			int numberTP=hitsTP.get(mDecoy.massbankQuery.get(i)).size();
			int numberFP=hitsFP.get(mDecoy.massbankQuery.get(i)).size();
			if(numberTP==0&&numberFP==0)none.add(mDecoy.massbankQuery.get(i));
			if(numberTP!=0&&numberFP==0)onlyTP.add(mDecoy.massbankQuery.get(i));
			if(numberTP==0&&numberFP!=0)onlyFP.add(mDecoy.massbankQuery.get(i));
			if(numberTP!=0&&numberFP!=0)both.add(mDecoy.massbankQuery.get(i));
		}
		
		System.out.println("number only TP: "+ onlyTP.size());
		System.out.println("number only FP: "+ onlyFP.size());
		System.out.println("number both: "+ both.size());
		System.out.println("number none: "+ none.size());
		
		double currentFDR=1.0*onlyFP.size()/(onlyFP.size()+onlyTP.size()+both.size());
		
		System.out.println("current number FDR: "+ currentFDR);

		int numberDeletionsBoth=(int)Math.ceil(percentage*(onlyFP.size()+onlyTP.size()+both.size())-onlyFP.size());
		System.out.println("number deletions in both recommended: "+numberDeletionsBoth);
		
		numberDeletionsBoth=Math.min(numberDeletionsBoth, both.size());
		
		System.out.println("number deletions in both possible: "+numberDeletionsBoth);

		Random random=new Random();

//		for(int i=0;i<numberDeletionsBoth;i++){
//			int r=random.nextInt(both.size());
//			MassBank mb=both.get(r);
//			result.put(mb, hitsTP.get(mb));
//			both.remove(mb);
//			onlyFP.add(mb);
//		}
		
		currentFDR=1.0*onlyFP.size()/(onlyFP.size()+onlyTP.size()+both.size());

		System.out.println("number only TP after deletions: "+ onlyTP.size());
		System.out.println("number only FP after deletions: "+ onlyFP.size());
		System.out.println("number both after deletions: "+ both.size());
		System.out.println("number none after deletions: "+ none.size());
		
		System.out.println("FDR after deletions: "+ currentFDR);
		
		if(currentFDR<percentage){
			double numberDeletionsTP=(int)Math.ceil(onlyTP.size()-(onlyFP.size()-percentage*(onlyFP.size()+both.size()))/percentage);
			System.out.println("number deletions in TPs: "+ numberDeletionsTP);
			
			for(int i=0;i<numberDeletionsTP;i++){
				if(onlyTP.size()==0)break;
				int r=random.nextInt(onlyTP.size());
				MassBank mb=onlyTP.get(r);
				result.put(mb, hitsTP.get(mb));
				onlyTP.remove(mb);
				none.add(mb);
			}
		}
		
		System.out.println("number only TP: "+ onlyTP.size());
		System.out.println("number only FP: "+ onlyFP.size());
		System.out.println("number both: "+ both.size());
		System.out.println("number none: "+ none.size());
		
		currentFDR=1.0*onlyFP.size()/(onlyFP.size()+onlyTP.size()+both.size());
		
		System.out.println("current number FDR: "+ currentFDR);
		

		return result;
	}

	private static void getRankStatistic(File searchResultsFolder, File graphicsFolder) throws IOException, Exception {

		TreeMap<String, MassBank> massbankFilesOriginal=Utils.getAllDBFiles(originalDBFolders, ser, false, false);

		for (int i=0;i<querys.length;i++){
			String queryLong=querys[i];			

			for(String searchMethod:searchMethods){
				Map<String, Map<DatabaseType, Double>[]> result=new HashMap<String, Map<DatabaseType, Double>[]>();
				Map<String, Map<DatabaseType, Double>[]> resultDifferentMassesExcluded=new HashMap<String, Map<DatabaseType, Double>[]>();

				for(String s:decoyMethods){
					if(s.equals("Original"))continue;

					File graphicsFile=new File(graphicsFolder.getAbsolutePath()+sep+searchMethod);
					if (!graphicsFile.exists())graphicsFile.mkdirs();

					String ending=(searchMethod.equals("DBSearch")?".txt":".csv");
					File original=new File(searchResultsFolder.getAbsolutePath()+sep+searchMethod+sep+"pos"+sep+queryLong+"Original"+ending);
					if(!original.exists())continue;

					SimilarityMatrix mOriginal=new SimilarityMatrix(original, massbankFilesOriginal);

					File files=new File(searchResultsFolder.getAbsolutePath()+sep+searchMethod+sep+"pos"+sep+queryLong+s+ending);					
					List<File> searchresultsFoldersDecoy = getSearchResultsFiles(files);
					if(searchresultsFoldersDecoy.isEmpty())continue;

					List<SimilarityMatrix> mDecoys=new ArrayList<SimilarityMatrix>();
					for(File f:searchresultsFoldersDecoy){
						mDecoys.add(new SimilarityMatrix(f, massbankFilesOriginal));
					}					
					String name=files.getName()+".jpg";
					name=name.replaceAll(".csv","").replaceAll(".txt","");

					System.out.println("getRankStatistics: " + graphicsFile.getAbsolutePath()+sep+"Hitstatistics_RankDistribution_"+name);

//					for(SimilarityMatrix mDecoy:mDecoys)mDecoy.excludeSameInChi();
//					mOriginal.excludeSameInChi();
					Map<DatabaseType, Double>[] r=SimilarityMatrix.getRankDistribution(new File(graphicsFile.getAbsolutePath()+sep+"Hitstatistics_RankDistribution_"+name), mOriginal, mDecoys,20);
					result.put(name.replaceAll(".jpg", ""), r);

					System.out.println("getRankStatistics: " + graphicsFile.getAbsolutePath()+sep+"Hitstatistics_DifferentMassesExcluded_RankDistribution_"+name);

					String[] split=queryLong.split("-"); 
					boolean isTrees=false;
					if(split[0].contains("Trees")&&split[1].contains("Trees"))isTrees=true;				
					if(isTrees){
						for(SimilarityMatrix mDecoy:mDecoys)mDecoy.excludeDifferentFormulas();
						mOriginal.excludeDifferentFormulas();
					}else{
						for(SimilarityMatrix mDecoy:mDecoys)mDecoy.excludeDifferentMass(ppm, ae);
						mOriginal.excludeDifferentMass(ppm, ae);
					}						

					r=SimilarityMatrix.getRankDistribution(new File(graphicsFile.getAbsolutePath()+sep+"Hitstatistics_DifferentMassesExcluded_RankDistribution_"+name), mOriginal, mDecoys, 20);
					resultDifferentMassesExcluded.put(name.replaceAll(".jpg", ""), r);
					System.out.println();
				}


				for(Map<String, Map<DatabaseType, Double>[]> res:new Map[]{result, resultDifferentMassesExcluded}){
					String add="";
					if(res==resultDifferentMassesExcluded)add="DifferentMassesExcluded_";
					File outputFile=new File(graphicsFolder.getAbsolutePath()+sep+searchMethod+sep+"Hitstatistics_"+add+"RankDistribution_"+queryLong+".jpg");
					if(!outputFile.getParentFile().exists())outputFile.getParentFile().mkdirs();
					final DefaultCategoryDataset dataset = new DefaultCategoryDataset();
					for(Entry<String, Map<DatabaseType, Double>[]> e:res.entrySet()){
						for(int j=0;j<e.getValue().length;j++){
							double decoy=0.0;
							double original=0.0;
							if(e.getValue()[j].containsKey(DatabaseType.Decoy))decoy=e.getValue()[j].get(DatabaseType.Decoy);
							if(e.getValue()[j].containsKey(DatabaseType.Original))original=e.getValue()[j].get(DatabaseType.Original);
							if(decoy!=0.0||original!=0.0)dataset.addValue(decoy/(decoy+original)-0.5,e.getKey(),Integer.toString(j));
						}					
					}
					final JFreeChart chart = ChartFactory.createBarChart("Boxplot","Ranks","Percentage",dataset);
					ChartUtilities.saveChartAsJPEG(outputFile, chart, Math.min(2000,20*100), 1000);
				}
				System.out.println();
			}
		}
	}

	private static void checkEquality(File searchResultsFolder, File graphicsFolder) throws IOException, Exception {

		TreeMap<String, MassBank> massbankFilesOriginal=Utils.getAllDBFiles(originalDBFolders, ser, false, false);			

		//			for(String searchMethod:searchMethods){
		for(String searchMethod:new String[]{"Equality"}){

			for(String s:decoyMethods){
				if(s.equals("Original"))continue;

				File graphicsFile=new File(graphicsFolder.getAbsolutePath()+sep+searchMethod);
				if (!graphicsFile.exists())graphicsFile.mkdirs();
				File outputFile=new File(graphicsFile.getAbsolutePath()+sep+"EuqalityCheck_"+s+".txt");
				BufferedWriter bw=new BufferedWriter(new FileWriter(outputFile));

				for (int i=0;i<querys.length;i++){
					String queryLong=querys[i];

					String[] starts=new String[]{"AP","MP","BP","OP","QP","GP"};
					for(String ending:starts){
						if(queryLong.contains(ending)){
							queryLong=ending+queryLong.substring(2);
						}
					}

					String ending=(searchMethod.equals("DBSearch")?".txt":".csv");

					//					File original=new File(searchResultsFolder.getAbsolutePath()+sep+searchMethod+sep+"pos"+sep+queryLong+"Original"+ending);
					//					if(!original.exists())continue;
					//					SimilarityMatrix mOriginal=new SimilarityMatrix(original, massbankFilesOriginal);

					File files=new File(searchResultsFolder.getAbsolutePath()+sep+searchMethod+sep+"pos"+sep+queryLong+s+ending);					
					List<File> searchresultsFoldersDecoy = getSearchResultsFiles(files);
					if(searchresultsFoldersDecoy.isEmpty())continue;

					for(File f:searchresultsFoldersDecoy){
						bw.write(f.getAbsolutePath()+": ");
						int numberEqual=0;
						String which="(";
						SimilarityMatrix m=new SimilarityMatrix(f, massbankFilesOriginal);
						for(int l=0;l<m.similarity.length;l++){
							for(int n=0;n<m.similarity[l].length;n++){
								if(m.similarity[l][n]!=0){
									numberEqual++;
									which+=m.massbankQuery.get(l).massbankID+"-"+m.massbankDB.get(n).massbankID+" ";
								}
							}
						}
						which+=")";
						bw.write(Double.toString(numberEqual)+" "+which);
						bw.newLine();
					}

					String name=files.getName();
					name=name.replaceAll(".csv","").replaceAll(".txt","");

					System.out.println("checkEquality: " + s+" "+name+".txt");


					bw.newLine();

				}
				bw.close();


			}
		}
	}

	private static List<File> getSearchResultsFiles(File files) {
		List<File> searchresultsFoldersDecoy = new ArrayList<File>();
		if(files.exists()){
			searchresultsFoldersDecoy.add(files);
		}else{
			File allFiles=new File(files.getAbsolutePath().substring(0,files.getAbsolutePath().lastIndexOf(".")));
			if(allFiles.exists()&&allFiles.isDirectory()){
				for(File inFile:allFiles.listFiles()){
					searchresultsFoldersDecoy.add(inFile);
				}
			}
		}
		return searchresultsFoldersDecoy;
	}

}

class FileExtensionFilter implements FilenameFilter {
	private Set<String> exts = new HashSet<String>();

	/**
	 * @param extensions
	 *            a list of allowed extensions, without the dot, e.g.
	 *            <code>"xml","html","rss"</code>
	 */
	public FileExtensionFilter(String... extensions) {
		for (String ext : extensions) {
			exts.add("." + ext.toLowerCase().trim());
		}
	}

	public boolean accept(File dir, String name) {
		final Iterator<String> extList = exts.iterator();
		while (extList.hasNext()) {
			if (name.toLowerCase().endsWith(extList.next())) {
				return true;
			}
		}
		return false;
	}
}
