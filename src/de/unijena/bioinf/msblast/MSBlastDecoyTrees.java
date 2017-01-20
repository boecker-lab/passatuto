package de.unijena.bioinf.msblast;

//Generating Decoy DB

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Map.Entry;

import de.unijena.bioinf.babelms.dot.Graph;
import de.unijena.bioinf.decoy.decoytrees.DecoySpectrum;
import de.unijena.bioinf.decoy.model.decoytreeconstructors.ConditionalFastConstructor;
import de.unijena.bioinf.decoy.model.decoytreeconstructors.DecoySpectrumConstructor;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction.CHARGE;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction.DATASET;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction.METHOD;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction.POSTPROCESS;
public class MSBlastDecoyTrees {	

	public static void main(String args[]) throws Exception{
		Locale.setDefault(Locale.US);

		ParametersDecoySpectrumConstruction p=new ParametersDecoySpectrumConstruction(args);
		if(!p.isCorrectCombination()){
			System.err.println("\""+Arrays.toString(args)+"\" is not a correct parameter combination.");
			return;
		}

		Map<File,Graph> origDBMap=p.getAllGraphsInOriginal();		
		Map<Integer,String> inchis=getInchiKeysByID(p.getInchiKeysFile());
		p.getAllPeaksInOriginal();
				
		BufferedWriter bw=null;
		if(p.rt){
			File f=p.getOutputFileRunningTime();
			if(!f.getParentFile().exists())f.getParentFile().mkdirs();
			bw=new BufferedWriter(new FileWriter(f));
		}

		int i=0;
		long allTimes=0;
		if(p.m == ParametersDecoySpectrumConstruction.METHOD.CONDITIONALFAST){
			ConditionalFastConstructor.setM(p);
		}
		for(Entry<File,Graph> e:origDBMap.entrySet()){
			i++;
			String name=e.getKey().getName().replaceAll(".dot","").replaceAll(".ms","");
			System.out.println(i + " of "+origDBMap.size() +" processed ("+name+")");
//			if(!name.equals("qpos559"))continue;
			int number=Integer.parseInt(name.replaceAll("\\D", ""));
			String inchi=inchis.get(number);
			String acc=p.getMassbankID(number);


			DecoySpectrum tree=null;
			DecoySpectrumConstructor dtc=p.getDecoyTreeConstructor();				
			dtc.setOriginalTree(e.getValue());
			
			Graph tmpGraph=new Graph(e.getValue());
			
			boolean foundDifferentTree=false;
			int n=0;
			long current=System.currentTimeMillis();
			while(!foundDifferentTree&&n<100){
				tree=dtc.getDecoySpectrum();
				if(!tree.equalsGraph(tmpGraph)||p.m == ParametersDecoySpectrumConstruction.METHOD.ORIGINAL){
					foundDifferentTree=true;
				}else{
//					System.out.println("no different tree found, try again.");
				}
				n++;
			}
			long time=System.currentTimeMillis()-current;
			if(bw!=null){
				bw.write("time for processing "+name +" (in milliseconds): "+time);
				bw.newLine();
			}			
			allTimes+=(time);
			if(!foundDifferentTree)System.err.println("no different tree found for "+p.getOutputFileDot(acc).getAbsolutePath());
			


			//				System.out.println(p.getOutputFileDot(acc).getAbsolutePath());			
			if(!p.pp.contains(POSTPROCESS.PPDiff)&&tree.getParent()!=null){
				tree.writeAsDot(p.getOutputFileDot(acc));
			}

			tree.writeAsMassbank(p.getOutputFolderMassbank(), tree.getPeaks(), tree.getParent(), 0, acc, acc+"_"+p.getInstrumentName(), p.getInstrumentName(), inchi, p.isPositive(),p.pp.contains(POSTPROCESS.PPDiff));

		}
		if(bw!=null){
			bw.write("time for processing all compounds (in milliseconds): "+ allTimes);
			bw.newLine();
			bw.close();
		}		
		System.out.println(p.getDecoyTreeConstructor().getMethodNameLong()+" finished.");

		p.deleteAllGraphsInOriginalDB();		

	}

	public static Map<String,List<Integer>> getInchiKeysByInchi(File inchiKeysFile) throws Exception{
		Map<String,List<Integer>> inchis=new HashMap<String,List<Integer>>();
		BufferedReader br=new BufferedReader(new FileReader(inchiKeysFile));
		String line=br.readLine();
		List<String> header=Arrays.asList(line.split(","));
		while((line=br.readLine())!=null){
			String l[]=line.split("\",\"");
			int id=Integer.parseInt(l[header.indexOf("\"mid\"")].replaceAll("\"",""));
			String inchiTmp=l[header.indexOf("\"inchi\"")].replaceAll("\"","").replaceAll("\"","");
			if(!inchiTmp.isEmpty()){
				String inchi=inchiTmp.split("/")[1]+"/"+inchiTmp.split("/")[2];
				if(!inchis.containsKey(inchi))inchis.put(inchi, new ArrayList<Integer>());
				inchis.get(inchi).add(id);
			}
		}
		br.close();
		return inchis;
	}

	public static Map<Integer,String> getInchiKeysByID(File inchiKeysFile) throws Exception{
		Map<Integer,String> inchis=new HashMap<Integer,String>();
		BufferedReader br=new BufferedReader(new FileReader(inchiKeysFile));
		String line=br.readLine();
		List<String> header=Arrays.asList(line.split("\t",Integer.MAX_VALUE));
		while((line=br.readLine())!=null){
			String l[]=line.replaceAll("\"","").split("\t",Integer.MAX_VALUE);
			if(l.length!=header.size())l=line.split(",");
			int id=Integer.parseInt(l[header.indexOf("id")]);
			String inchiTmp=l[header.indexOf("inchi")].replaceAll("\"","");
			if(!inchiTmp.isEmpty()){
				String[] isplit=inchiTmp.split("/");
				if(isplit.length>=3&&!isplit[2].startsWith("c")){
					System.out.print(inchiTmp);
					System.out.println();
				}
				if(isplit.length>=4&&!isplit[3].startsWith("h")){
					System.out.print(inchiTmp);
					System.out.println();
				}
				String inchi=isplit[0]+"/"+isplit[1]+(isplit.length>=3&&(isplit[2].startsWith("c")||isplit[2].startsWith("h"))?("/"+isplit[2]):"")+(isplit.length>=4&&isplit[3].startsWith("h")?("/"+isplit[3]):"");
				inchis.put(id, inchi);
			}
		}
		br.close();
		return inchis;
	}
}
