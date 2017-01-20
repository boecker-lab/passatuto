package de.unijena.bioinf.statistics;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import de.unijena.bioinf.em.EMUtils;
import de.unijena.bioinf.em.Sample.Type;

public class pValues {

	static String base="U:\\MSBlast\\";

	public static void main(String args[]) throws Exception{
		for(String compounds : new String[]{"Compounds"}){
			//			for(boolean log:new boolean[]{false,true}){
			for(boolean log:new boolean[]{false}){
								for(String comparisonMethod : new String[]{"MassBank", "CosineDistance"}){
//				for(String comparisonMethod : new String[]{"MassBank"}){
										for(String hs : new String[]{"DifferentMassesExcluded","DMEHighFDR"}){
//					for(String hs : new String[]{"DifferentMassesExcluded"}){
						for(String dataset : new String[]{"APFilesOriginal-GPFiles","APFilesOriginal-GPTrees","APTreesOriginal-GPTrees","OPFilesOriginal-GPFiles","OPFilesOriginal-GPTrees","OPTreesOriginal-GPTrees"}){
//						for(String dataset : new String[]{"APFilesOriginal-GPFiles"}){

							for(String decoydataset : new String[]{"ConditionalFast","MixedSpectrum","RandomPeaks","RandomTree","Reroot"}){
//							for(String decoydataset : new String[]{"ConditionalFast"}){
								
								File inputFolder=new File("U:\\MSBlast\\qValues_All"+compounds+"RPP\\pos\\MassBank\\HitStatistic_"+hs+"_"+dataset+decoydataset);
								if(!inputFolder.exists())continue;

								int selMax=1;
								if(compounds.equals("Selections")){
									selMax=10;
								}

								int decoyMax=1;
								if(compounds.equals("Compounds")){
									decoyMax=10;
								}
								
								List<Double> TP=new ArrayList<Double>();
								List<Double> FP=new ArrayList<Double>();
								List<Double> All=new ArrayList<Double>();

								for(int sel=0;sel<selMax;sel++){
									String addSel="_Selection"+sel;
									
									for(int decoy=1;decoy<decoyMax;decoy++){
										String addDecoy="_"+decoy;

										File inputFile=new File("U:\\MSBlast\\qValues_All"+compounds+"RPP\\pos\\MassBank\\HitStatistic_"+hs+"_"+dataset+decoydataset+"\\QValues"+addSel+"_HitStatistic_"+hs+"_"+dataset+decoydataset+addDecoy+".hitlist");

										if(!inputFile.exists())continue;
										
										BufferedReader br=new BufferedReader(new FileReader(inputFile));
										br.readLine();br.readLine();
										String line;
										int allDecoys=0;
										while((line=br.readLine())!=null){
											String l[]=line.split("\t");
											if(l[1].equals("DecoyMatch"))allDecoys++;
										}
										br.close();
										
										br=new BufferedReader(new FileReader(inputFile));
										br.readLine();br.readLine();
										int numberDecoys=0;
										while((line=br.readLine())!=null){
											String l[]=line.split("\t");
											if(l[1].equals("DecoyMatch"))numberDecoys++;
											if(l[1].equals("TruePositiveMatch"))TP.add(1.0*numberDecoys/allDecoys);
											if(l[1].equals("FalsePositiveMatch"))FP.add(1.0*numberDecoys/allDecoys);
											if(!l[1].equals("DecoyMatch"))All.add(1.0*numberDecoys/allDecoys);
										}
										
										br.close();
									}
								}

								for(int bins:new int[]{10,20,30}){

									File outputFilePValuesAll=new File("U:\\MSBlast\\qValues_All"+compounds+"RPP\\pos\\MassBank\\HitStatistic_"+hs+"_"+dataset+decoydataset+"\\QValues_HitStatistic_"+hs+"_"+dataset+decoydataset+"_All_"+bins+"bins.pValues");
									EMUtils.writePValues(outputFilePValuesAll, "All", All, bins);
	
									File outputFilePValuesFP=new File("U:\\MSBlast\\qValues_All"+compounds+"RPP\\pos\\MassBank\\HitStatistic_"+hs+"_"+dataset+decoydataset+"\\QValues_HitStatistic_"+hs+"_"+dataset+decoydataset+"_FP_"+bins+"bins.pValues");
									EMUtils.writePValues(outputFilePValuesFP, "FP", FP, bins);
	
									File outputFilePValuesTP=new File("U:\\MSBlast\\qValues_All"+compounds+"RPP\\pos\\MassBank\\HitStatistic_"+hs+"_"+dataset+decoydataset+"\\QValues_HitStatistic_"+hs+"_"+dataset+decoydataset+"_TP_"+bins+"bins.pValues");
									EMUtils.writePValues(outputFilePValuesTP, "TP", TP, bins);
								}
								
								System.out.print("");
							}
						}
					}
				}
			}
		}


	}
}
