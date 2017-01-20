package de.unijena.bioinf.em;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Map.Entry;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.GumbelDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.commons.math3.distribution.WeibullDistribution;

import de.unijena.bioinf.em.Sample.Type;
import de.unijena.bioinf.statistics.MainStatistics;

public class EMMain {

	static String R="C:\\Program Files\\R\\R-3.0.2\\bin\\RScript.exe";
	static String EMScript="U:\\MSBlast\\EM.R";
	static String base="U:\\MSBlast\\";

	static InterfaceMLDistribution[] allDistributions=new InterfaceMLDistribution[]{
		new MLNormalDistribution(false, 1E-5),
		//		new MLLogNormalDistribution(false, 1E-5, false),
		//		new MLLogNormalDistribution(true, 1E-5, false),
		//		new MLBetaDistribution(false, 1E-5),
		//		new MLBetaDistribution(true, 1E-5),
		new MLGammaDistribution(false, 1E-5),
		new MLGammaDistribution(true, 1E-5),
		new MLExponentialDistribution(false, 1E-5),
		new MLExponentialDistribution(true, 1E-5),
		new MLGumbelDistribution(false, 1E-5),
		new MLGumbelDistribution(true, 1E-5),
		new MLWeibullDistribution(false, 1E-5),
		new MLWeibullDistribution(true, 1E-5),
		new MLLogisticDistribution(false, 1E-5)
	};

	static InterfaceMLDistribution[] selectedTPDistributions=new InterfaceMLDistribution[]{
		new MLGammaDistribution(true, 1E-5),
		new MLGumbelDistribution(true, 1E-5),
		new MLWeibullDistribution(true, 1E-5)
	};

	static InterfaceMLDistribution[] selectedFPDistributions=new InterfaceMLDistribution[]{
		new MLGammaDistribution(false, 1E-5)
	};

	static InterfaceMLDistribution[] selectedTPDistributions_LOG=new InterfaceMLDistribution[]{
		new MLGammaDistribution(false, 1E-5),
		new MLGumbelDistribution(false, 1E-5),
		new MLWeibullDistribution(false, 1E-5),
		new MLNormalDistribution(false, 1E-5), 
		new MLLogisticDistribution(false, 1E-5)
	};

	static InterfaceMLDistribution[] selectedFPDistributions_LOG=new InterfaceMLDistribution[]{
		new MLGammaDistribution(false, 1E-5),
		new MLWeibullDistribution(false, 1E-5)
	};

	public static void main(String[] args) throws Exception{

		//		mainTest();

//		mainReal("single");
		mainReal("combined");
		//		getAllQValuesAverage("");

		//				mainReal("variable_component_mixture_model");
		//		getAllQValuesAverage("VCMM");
	}


	public static void mainReal(String what) throws Exception{		

		int maxIterations=1000;
		double minChangeLikelihood=1E-20;
//		for(String compounds : new String[]{"Compounds","Selections"}){
		for(String compounds : new String[]{"Compounds"}){
//			for(boolean log:new boolean[]{false,true}){
			for(boolean log:new boolean[]{false}){
				boolean discardZeroEntries=log;
				BufferedWriter bw=new BufferedWriter(new FileWriter(base+"qValues_PEP_All"+compounds+"RPP\\EM_Likelihoods_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+what+".txt"));
				for(String comparisonMethod : new String[]{"MassBank", "CosineDistance"}){
//				for(String comparisonMethod : new String[]{"MassBank"}){
					for(String hs : new String[]{"DifferentMassesExcluded","DMEHighFDR"}){
//					for(String hs : new String[]{"DifferentMassesExcluded"}){
						for(String dataset : new String[]{"APFilesOriginal-GPFiles","APFilesOriginal-GPTrees","APTreesOriginal-GPTrees","OPFilesOriginal-GPFiles","OPFilesOriginal-GPTrees","OPTreesOriginal-GPTrees"}){
//						for(String dataset : new String[]{"APFilesOriginal-GPFiles"}){

							List<File> outputFilesQValuesPEP=new ArrayList<File>();
							List<File> outputFilesQValuesEFDR=new ArrayList<File>();
							int selMax=1;
							if(compounds.equals("Selections")){
								selMax=10;
							}

							for(int sel=0;sel<selMax;sel++){
								String add="";
								if(compounds.equals("Selections")){
									add="Selection"+sel+"_";
								}

								File inputFile=new File(base+"qValues_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+hs+"_"+dataset+"RandomPeaks\\QValues_Selection"+sel+"_HitStatistic_"+hs+"_"+dataset+"RandomPeaks_1.hitlist");
								if(!inputFile.exists())continue;


								List<Sample> sample_sep=EMUtils.getScores(inputFile, log, discardZeroEntries);

								if(what.contains("single")){

									for(int i=0;i<sample_sep.size();i++){
										List<EMResult> results=new ArrayList<EMResult>();

										Sample sample=new Sample(sample_sep.get(i));								

										File outputFileData=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"Single\\SampleData_"+add+(i==1?"TP":"FP")+"_"+hs+"_"+dataset+".txt");
										EMUtils.writeData(outputFileData, sample.sample, null);

										for(InterfaceMLDistribution d:allDistributions){	
											System.out.println("processing "+compounds+" "+comparisonMethod+" "+add+" "+dataset+" "+hs+" "+" "+d.getDescription()+" "+(i==1?"TP":"FP")+ " log:"+log+" discardZeroEntries:"+discardZeroEntries+"...");									
											EM em=new EM(sample, d, maxIterations, minChangeLikelihood);
											try{
												EMResult r=em.doEM(false);

												File outputFileDistribution=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"Single\\Distribution_"+add+d.getDescription()+"_"+(i==1?"TP":"FP")+"_"+hs+"_"+dataset+".txt");
												EMUtils.writeDistribution(outputFileDistribution, sample.sample, r.distributions, r.parameters);

												File outputFilePDF=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"Single\\HistogramDistribution_"+add+d.getDescription()+"_"+(i==1?"TP":"FP")+"_"+hs+"_"+dataset+".pdf");
												if(!outputFilePDF.getParentFile().exists())outputFilePDF.getParentFile().mkdirs();

												String command="\""+R+"\" \""+EMScript+"\" "+outputFilePDF.getAbsolutePath()+" "+outputFileData.getAbsolutePath()+" "+outputFileDistribution.getAbsolutePath()+"\"";
												Process proc=Runtime.getRuntime().exec(command);

												System.out.println(r.likelihood);
												System.out.println();
												if(!Double.isNaN(r.likelihood))results.add(r);
											}catch(Exception e){
												System.err.println(e.getMessage() +" while proceeding "+ compounds+" "+comparisonMethod+" "+add+" "+dataset+" "+hs+" "+d.getDescription()+" "+(i==1?"TP":"FP")+ " log:"+log+" discardZeroEntries:"+discardZeroEntries);
												for(StackTraceElement el:e.getStackTrace()){
													System.err.println(el.toString());
												}
											}
										}
										Collections.sort(results);

										bw.write(compounds+" "+comparisonMethod+" "+dataset+" "+hs+" "+add+(i==1?"TP":"FP")+ " log:"+log+" discardZeroEntries:"+discardZeroEntries);
										bw.newLine();
										for(EMResult r:results){
											bw.write(r.toString());
											bw.newLine();
											bw.flush();
										}
										bw.newLine();
									}
									bw.newLine();									
								}

								if(what.contains("combined")){

									Sample sample=new Sample(sample_sep);
									List<Double> qValuesReal=sample.getQValues();
									List<EMResult> results=new ArrayList<EMResult>();

									File outputFileData=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"Combined\\SampleData_"+add+hs+"_"+dataset+".txt");
									EMUtils.writeData(outputFileData, sample.sample, sample.types);
									
									InterfaceMLDistribution[] selectedFPDistributions_Tmp=(log?selectedFPDistributions_LOG:selectedFPDistributions);
									InterfaceMLDistribution[] selectedTPDistributions_Tmp=(log?selectedTPDistributions_LOG:selectedTPDistributions);

									for(InterfaceMLDistribution d1:selectedFPDistributions_Tmp){
										for(InterfaceMLDistribution d2:selectedTPDistributions_Tmp){
											//										if(d1.getDescription().compareTo(d2.getDescription())>0)continue;
											System.out.println("processing "+compounds+" "+comparisonMethod+" "+add+" "+dataset+" "+hs+" "+d1.getDescription()+" "+d2.getDescription()+ " log:"+log+" discardZeroEntries:"+discardZeroEntries+"...");
											List<InterfaceMLDistribution> distributions=new ArrayList<InterfaceMLDistribution>();
											distributions.add(d1);
											distributions.add(d2);

											EM em=new EM(sample, distributions, maxIterations, minChangeLikelihood);
											try{
												EMResult r=em.doEM(false);

												File outputFileDistribution=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"Combined\\Distribution_"+add+hs+"_"+dataset+"_"+d1.getDescription()+"_"+d2.getDescription()+".txt");
												EMUtils.writeDistribution(outputFileDistribution, sample.sample, r.distributions, r.parameters);

												File outputFilePDF_FP=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"Combined\\HistogramDistribution_"+add+"FP_"+hs+"_"+dataset+"_"+d1.getDescription()+"_"+d2.getDescription()+".pdf");
												File outputFilePDF_TP=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"Combined\\HistogramDistribution_"+add+"TP_"+hs+"_"+dataset+"_"+d1.getDescription()+"_"+d2.getDescription()+".pdf");
												File outputFilePDF=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"Combined\\HistogramDistribution_"+add+hs+"_"+dataset+"_"+d1.getDescription()+"_"+d2.getDescription()+".pdf");
												if(!outputFilePDF.getParentFile().exists())outputFilePDF.getParentFile().mkdirs();

												String command="\""+R+"\" \""+EMScript+"\" "+outputFilePDF.getAbsolutePath()+" "+outputFileData.getAbsolutePath()+" "+outputFileDistribution.getAbsolutePath()+" "+outputFilePDF_FP.getAbsolutePath()+" "+outputFilePDF_TP.getAbsolutePath()+"\"";
												Process proc=Runtime.getRuntime().exec(command);


												List<Double> qValuesPred=r.getQValuesByPEP();									
												r.deviation_QValuesByPEP=EMUtils.deviation(qValuesReal,qValuesPred);
												File outputFileQValues=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"PEP\\QValues_HitStatistic_"+add+hs+"_"+dataset+"_"+d1.getDescription()+"_"+d2.getDescription()+".txt");
												EMUtils.writeQValues(outputFileQValues, r.toString(), qValuesReal, qValuesPred);										

												qValuesPred=r.getQValuesByEDFR();										
												r.deviation_QValuesByEFDR=EMUtils.deviation(qValuesReal,qValuesPred);
												outputFileQValues=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"EFDR\\QValues_HitStatistic_"+add+hs+"_"+dataset+"_"+d1.getDescription()+"_"+d2.getDescription()+".txt");
												EMUtils.writeQValues(outputFileQValues, r.toString(), qValuesReal, qValuesPred);
												
												for(int bins:new int[]{10,20,30}){
												
													Map<Type,List<Double>> pValuesPredMap=r.getPValues(sample);
													List<Double> pValuesPred=new ArrayList<Double>();
													pValuesPred.addAll(pValuesPredMap.get(Type.FalsePositiveMatch));
													pValuesPred.addAll(pValuesPredMap.get(Type.TruePositiveMatch));
													File outputFilePValuesAll=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"PValues\\QValues_HitStatistic_"+add+hs+"_"+dataset+"_"+d1.getDescription()+"_"+d2.getDescription()+"_All_"+bins+"bins.txt");
													EMUtils.writePValues(outputFilePValuesAll, r.toString()+" All", pValuesPred, bins);
													
													File outputFilePValuesFP=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"PValues\\QValues_HitStatistic_"+add+hs+"_"+dataset+"_"+d1.getDescription()+"_"+d2.getDescription()+"_FP_"+bins+"bins.txt");
													EMUtils.writePValues(outputFilePValuesFP, r.toString()+" FP", pValuesPredMap.get(Type.FalsePositiveMatch), bins);
													
													File outputFilePValuesTP=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"PValues\\QValues_HitStatistic_"+add+hs+"_"+dataset+"_"+d1.getDescription()+"_"+d2.getDescription()+"_TP_"+bins+"bins.txt");
													EMUtils.writePValues(outputFilePValuesTP, r.toString()+" TP", pValuesPredMap.get(Type.TruePositiveMatch), bins);
	
												}
												
												System.out.println(r.likelihood+" "+r.deviation_QValuesByPEP+" "+r.deviation_QValuesByEFDR);
												System.out.println();

												if(!Double.isNaN(r.likelihood))results.add(r);																
											}catch(Exception e){	
												System.err.println(e.getMessage() +" while proceeding "+ compounds+" "+comparisonMethod+" "+add+" "+dataset+" "+hs+" "+d1.getDescription()+" "+d2.getDescription()+ " log:"+log+" discardZeroEntries:"+discardZeroEntries);
												for(StackTraceElement el:e.getStackTrace()){
													System.err.println(el.toString());
												}
											}

										}
									}
									Collections.sort(results);

									bw.write(compounds+" "+comparisonMethod+" "+add+" "+dataset+" "+hs+ " log:"+log+" discardZeroEntries:"+discardZeroEntries);
									bw.newLine();
									for(EMResult r:results){
										bw.write(r.toString());
										bw.newLine();
										bw.flush();
									}
									bw.newLine();

									EMResult r=results.get(0);

									File outputFileDistribution=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"Combined\\Distribution_"+add+hs+"_"+dataset+".txt");
									EMUtils.writeDistribution(outputFileDistribution, sample.sample, r.distributions, r.parameters);

									File outputFilePDF_FP=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"Combined\\HistogramDistribution_"+add+"FP_"+hs+"_"+dataset+".pdf");
									File outputFilePDF_TP=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"Combined\\HistogramDistribution_"+add+"TP_"+hs+"_"+dataset+".pdf");
									File outputFilePDF=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"Combined\\HistogramDistribution_"+add+hs+"_"+dataset+".pdf");

									String command="\""+R+"\" \""+EMScript+"\" "+outputFilePDF.getAbsolutePath()+" "+outputFileData.getAbsolutePath()+" "+outputFileDistribution.getAbsolutePath()+" "+outputFilePDF_FP.getAbsolutePath()+" "+outputFilePDF_TP.getAbsolutePath()+"\"";
									Process proc=Runtime.getRuntime().exec(command);

									List<Double> qValuesPred=r.getQValuesByPEP();
									File outputFileQValues=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"PEP\\QValues_HitStatistic_"+add+hs+"_"+dataset+"PEP.qValueMeanAverage");
									outputFilesQValuesPEP.add(outputFileQValues);
									EMUtils.writeQValues(outputFileQValues, r.toString(), qValuesReal, qValuesPred);

									qValuesPred=r.getQValuesByEDFR();
									outputFileQValues=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"EFDR\\QValues_HitStatistic_"+add+hs+"_"+dataset+"EFDR.qValueMeanAverage");
									outputFilesQValuesEFDR.add(outputFileQValues);
									EMUtils.writeQValues(outputFileQValues, r.toString(), qValuesReal, qValuesPred);
									
									for(int bins:new int[]{10,20,30}){
										Map<Type,List<Double>> pValuesPredMap=r.getPValues(sample);
										List<Double> pValuesPred=new ArrayList<Double>();
										pValuesPred.addAll(pValuesPredMap.get(Type.FalsePositiveMatch));
										pValuesPred.addAll(pValuesPredMap.get(Type.TruePositiveMatch));
										File outputFilePValuesAll=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"PValues\\QValues_HitStatistic_"+add+hs+"_"+dataset+"_All_"+bins+"bins.txt");
										EMUtils.writePValues(outputFilePValuesAll, r.toString()+" All", pValuesPred, bins);
										
										File outputFilePValuesFP=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"PValues\\QValues_HitStatistic_"+add+hs+"_"+dataset+"_FP_"+bins+"bins.txt");
										EMUtils.writePValues(outputFilePValuesFP, r.toString()+" FP", pValuesPredMap.get(Type.FalsePositiveMatch), bins);
										
										File outputFilePValuesTP=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"PValues\\QValues_HitStatistic_"+add+hs+"_"+dataset+"_TP_"+bins+"bins.txt");
										EMUtils.writePValues(outputFilePValuesTP, r.toString()+" TP", pValuesPredMap.get(Type.TruePositiveMatch), bins);
									}
								}										

								if(what.contains("variable_component_mixture_model")){
									Sample sample=new Sample(sample_sep);
									List<Double> qValuesReal=sample.getQValues();

									File outputFileData=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+add+hs+"_"+dataset+"CombinedVCMM\\SampleData_"+hs+"_"+dataset+".txt");
									EMUtils.writeData(outputFileData, sample.sample, null);

									System.out.println("processing "+compounds+" "+comparisonMethod+" "+add+" "+dataset+" "+hs+ " log:"+log+" discardZeroEntries:"+discardZeroEntries+" variable component mixture model...");
									List<InterfaceMLDistribution> distributions=new ArrayList<InterfaceMLDistribution>();
									for(int i=0;i<20;i++){
										distributions.add(new MLNormalDistribution(false,0));
									}


									EM em=new EM(sample, distributions, maxIterations, minChangeLikelihood);
									try{
										EMResult r=em.doEM(true);

										File outputFileDistribution=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+add+hs+"_"+dataset+"CombinedVCMM\\Distribution_"+hs+"_"+dataset+".txt");
										EMUtils.writeDistribution(outputFileDistribution, sample.sample, r.distributions, r.parameters);

										File outputFilePDF=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+add+hs+"_"+dataset+"CombinedVCMM\\HistogramDistribution_"+hs+"_"+dataset+".pdf");
										if(!outputFilePDF.getParentFile().exists())outputFilePDF.getParentFile().mkdirs();

										String command="\""+R+"\" \""+EMScript+"\" "+outputFilePDF.getAbsolutePath()+" "+outputFileData.getAbsolutePath()+" "+outputFileDistribution.getAbsolutePath()+"\"";
										Process proc=Runtime.getRuntime().exec(command);

										List<Double> qValuesPred=r.getQValuesByPEP();
										r.deviation_QValuesByPEP=EMUtils.deviation(qValuesReal,qValuesPred);
										File outputFileQValues=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+add+hs+"_"+dataset+"PEPVCMM\\QValues_HitStatistic_"+hs+"_"+dataset+"PEPVCMM.qValueMeanAverage");
										outputFilesQValuesPEP.add(outputFileQValues);
										EMUtils.writeQValues(outputFileQValues, r.toString(), qValuesReal, qValuesPred);								

										qValuesPred=r.getQValuesByEDFR();
										r.deviation_QValuesByEFDR=EMUtils.deviation(qValuesReal,qValuesPred);
										outputFileQValues=new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+add+hs+"_"+dataset+"EFDRVCMM\\QValues_HitStatistic_"+hs+"_"+dataset+"EFDRVCMM.qValueMeanAverage");
										outputFilesQValuesEFDR.add(outputFileQValues);
										EMUtils.writeQValues(outputFileQValues, r.toString(), qValuesReal, qValuesPred);								

										System.out.println(r.likelihood+" "+r.deviation_QValuesByPEP+" "+r.deviation_QValuesByEFDR);
										System.out.println();

										if(!Double.isNaN(r.likelihood)){
											bw.write(compounds+" "+comparisonMethod+" "+add+" "+dataset+" "+hs+ " log:"+log+" discardZeroEntries:"+discardZeroEntries);
											bw.newLine();
											bw.write(r.toString());
											bw.newLine();
											bw.flush();
											bw.newLine();
										}
									}catch(Exception e){	
										System.err.println(e.getMessage() +" while proceeding "+ compounds+" "+comparisonMethod+" "+add+" "+dataset+" "+hs+" variable component mixture model");
										for(StackTraceElement el:e.getStackTrace()){
											System.err.println(el.toString());
										}
									}
								}
							}
							bw.newLine();

							if(compounds.equals("Selections")&&(what.contains("combined")||what.contains("variable_component_mixture_model"))){
								String add="";
								if(what.contains("variable_component_mixture_model"))add="VCMM";
								getQValuesAverage(outputFilesQValuesPEP, new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"PEP"+add+"\\QValues_HitStatistic_"+hs+"_"+dataset+"PEP"+add+".qValueMeanAverage"));
								getQValuesAverage(outputFilesQValuesEFDR, new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+(discardZeroEntries?"DZE_":"")+(log?"Log_":"")+hs+"_"+dataset+"EFDR"+add+"\\QValues_HitStatistic_"+hs+"_"+dataset+"EFDR"+add+".qValueMeanAverage"));
							}
							//							}
							//							bw.newLine();
						}
						bw.newLine();
					}
					bw.newLine();
				}
				bw.newLine();
				bw.close();
			}
		}
	}

	private static void getStandardOutput(Process proc) throws IOException {
		BufferedReader stdInput = new BufferedReader(new 
				InputStreamReader(proc.getInputStream()));

		BufferedReader stdError = new BufferedReader(new 
				InputStreamReader(proc.getErrorStream()));

		// read the output from the command
		System.out.println("Here is the standard output of the command:\n");
		String s = null;
		while ((s = stdInput.readLine()) != null) {
			System.out.println(s);
		}

		// read any errors from the attempted command
		System.out.println("Here is the standard error of the command (if any):\n");
		while ((s = stdError.readLine()) != null) {
			System.out.println(s);
		}
	}

	public static void mainTest() throws Exception{
		int maxIterations=10000;

		List<List<Double>> sample_sep=EMUtils.getSimulatedScores(new int[]{1000,2000},Arrays.asList(new RealDistribution[]{
				new BetaDistribution(2,2),
				new BetaDistribution(5,1)
		}));

		//		List<List<Double>> sample_sep=EMUtils.getSimulatedScores(new int[]{1000,2000},Arrays.asList(new RealDistribution[]{
		//				new GammaDistribution(2,2),
		//				new GammaDistribution(7.5,1)
		//		}));
		List<Double> values=new ArrayList<Double>();
		List<InterfaceMLDistribution> distributions=new ArrayList<InterfaceMLDistribution>();
		for(List<Double> s:sample_sep){
			values.addAll(s);			
			distributions.add(new MLBetaDistribution(false, 0));
			//			distributions.add(new MLGammaDistribution(false, 0));
		}

		Sample sample=new Sample(values, null, null);
		//		EM em=new EM(sample, distributions, maxIterations, minChangeLikelihood);
		EM em=new EM(sample, distributions, maxIterations, Double.NEGATIVE_INFINITY);
		EMResult r=em.doEM(true);
		System.out.println(r.likelihood);
	}


	public static void getAllQValuesAverage(String type) throws Exception{
		for(String compounds : new String[]{"Selections"}){
			for(String comparisonMethod : new String[]{"MassBank", "CosineDistance"}){
				//			for(String comparisonMethod : new String[]{"MassBank"}){
				for(String hs : new String[]{"DifferentMassesExcluded","DMEHighFDR"}){
					//				for(String hs : new String[]{"DifferentMassesExcluded"}){
					for(String dataset : new String[]{"APFilesOriginal-GPFiles","APFilesOriginal-GPTrees","APTreesOriginal-GPTrees","OPFilesOriginal-GPFiles","OPFilesOriginal-GPTrees","OPTreesOriginal-GPTrees"}){
						//					for(String dataset : new String[]{"APFilesOriginal-GPFiles"}){


						List<File> outputFilesQValuesPEP=new ArrayList<File>();
						List<File> outputFilesQValuesEFDR=new ArrayList<File>();
						int selMax=0;
						if(compounds.equals("Selections")){
							selMax=10;
						}
						for(int sel=0;sel<selMax;sel++){
							String add="";
							if(compounds.equals("Selections")){
								add="Selection"+sel+"_";
							}

							outputFilesQValuesPEP.add(new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+hs+"_"+dataset+"PEP"+type+"\\QValues_HitStatistic_"+add+hs+"_"+dataset+"PEP"+type+".qValueMeanAverage"));						
							outputFilesQValuesEFDR.add(new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+hs+"_"+dataset+"EFDR"+type+"\\QValues_HitStatistic_"+add+hs+"_"+dataset+"EFDR"+type+".qValueMeanAverage"));
						}

						getQValuesAverage(outputFilesQValuesPEP, new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+hs+"_"+dataset+"PEP"+type+"\\QValues_HitStatistic_"+hs+"_"+dataset+"PEP"+type+".qValueMeanAverage"));
						getQValuesAverage(outputFilesQValuesEFDR, new File(base+"qValues_PEP_All"+compounds+"RPP\\pos\\"+comparisonMethod+"\\HitStatistic_"+hs+"_"+dataset+"EFDR"+type+"\\QValues_HitStatistic_"+hs+"_"+dataset+"EFDR"+type+".qValueMeanAverage"));

					}
				}
			}
		}
	}

	private static void getQValuesAverage(List<File> inputFiles, File outputFile) throws Exception{

		List<List<double[]>> r=MainStatistics.getPoints(inputFiles, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, true, true, "and", false, false);

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

			double[] as=MainStatistics.getAverageAndStandardDeviation(points);
			p.add(points);
			res.add(new double[]{k,as[0],as[1]});
		}		
		BufferedWriter bw=new BufferedWriter(new FileWriter(outputFile));
		for(int i=0;i<inputFiles.size();i++){				
			File f=inputFiles.get(i);
			BufferedReader br=new BufferedReader(new FileReader(f));
			if(i!=0)bw.write("\t");
			bw.write(f.getName()+": " + br.readLine());
			br.close();
		}
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
	}

}
