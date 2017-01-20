package de.unijena.bioinf.treecomparison;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.nio.file.Files;

//generating tree alignments

public class Main {
	
	static String queryLong="Agilent";
	static String DBLong="Metlin";


	public static void main(String[] args) throws Exception{
		
		String queryShort=queryLong.substring(0,1);
		String DBShort=DBLong.substring(0,1);
		
		String[] withs=new String[]{
				DBLong+"PositiveOriginal",
				DBLong+"PositiveReroot"};
		String[] outs=new String[]{
				queryShort+"POriginal-"+DBShort+"POriginal",
				queryShort+"POriginal-"+DBShort+"PReroot"};
		String aligns=queryLong+"PositiveOriginal";

		String dataset="D:\\metabo_tandem_ms\\MSBlast\\";


		for(int i=0;i<withs.length;i++){

			String align=dataset+"DB\\DataBases\\"+queryLong+"\\pos\\"+aligns+"\\dot\\";
			String with=dataset+"DB\\DataBases\\"+DBLong+"\\pos\\"+withs[i]+"\\dot\\";
			String output=dataset+"searchResults\\TreeAlignment\\pos\\"+outs[i];
			
			File tmpAlignmentMatrix=new File(System.getProperty("java.io.tmpdir")+"\\tmp.csv");

			System.out.println("starting... ");    	
			String s1="-z -x -j " +
					"-m "+tmpAlignmentMatrix.getAbsolutePath()+" " +
					"--align "+align+" " +
					"--with "+with;
			args=s1.split(" ");
			new de.unijena.bioinf.ftalign.Main().run(args);
			
			
			BufferedReader br=new BufferedReader(new FileReader(tmpAlignmentMatrix));
			BufferedWriter bw=new BufferedWriter(new FileWriter(output+".csv"));
			
			bw.write("Query: "+align.replaceAll("dot","massbank")+";");
			bw.newLine();
			bw.write("DB: "+with.replaceAll("dot","massbank")+";");
			bw.newLine();
			bw.newLine();
			String line;
			while((line=br.readLine())!=null){
				bw.write(line.replaceAll("\"","").replaceAll("scores", ""));
				bw.newLine();
			}
			br.close();
			bw.close();

			//			System.out.println("starting... ");    	
			//			s1="-z -f -x -j " +
			//					"-m "+output+"FP.csv " +
			//					"--align "+align+" " +
			//					"--with "+with;
			//			args=s1.split(" ");
			//			new de.unijena.bioinf.ftalign.Main().run(args); 	
			//
			//			System.out.println("...finished");

//			File tmpDirQueryAndDB=new File(System.getProperty("java.io.tmpdir")+"\\QueryAndDB");
//			if(tmpDirQueryAndDB.exists()){
//				for (File f:tmpDirQueryAndDB.listFiles())f.delete();
//				tmpDirQueryAndDB.delete();
//			}
//			tmpDirQueryAndDB.mkdirs();
//			for(File f:new File(align).listFiles()){
//				if(f.getName().endsWith(".dot"))Files.copy(f.toPath(), new File(tmpDirQueryAndDB.getAbsolutePath()+"\\"+f.getName()).toPath());
//			}
//			for(File f:new File(with).listFiles()){
//				if(f.getName().endsWith(".dot"))Files.copy(f.toPath(), new File(tmpDirQueryAndDB.getAbsolutePath()+"\\"+f.getName()).toPath());
//			}
//
//			tmpAlignmentMatrix=new File(System.getProperty("java.io.tmpdir")+"\\tmp.csv");
//
//			System.out.println("starting... ");    	
//			s1="-z -f -x -j " +
//					"-m "+tmpAlignmentMatrix.getAbsolutePath()+" " +
//					"--align "+tmpDirQueryAndDB.getAbsolutePath()+" " +
//					"--with "+align;
//			args=s1.split(" ");
//			new de.unijena.bioinf.ftalign.Main().run(args); 
//
//			br=new BufferedReader(new FileReader(tmpAlignmentMatrix));
//			bw=new BufferedWriter(new FileWriter(output+"FP.csv"));
//
//			File[] left=new File(align).listFiles();
//			File[] right=new File(with).listFiles();
//
//			double[][] scores=new double[left.length][right.length];
//			String[] header=br.readLine().split(",");
//			int[] headerIndexes=new int[header.length];
//			for(int j=0;j<header.length;j++)headerIndexes[j]=Integer.MIN_VALUE;
//			for(int j=0;j<header.length;j++){
//				for(int k=0;k<left.length;k++){
//					if(header[j].replaceAll("\"","").equals(left[k].getName().replaceAll(".dot",""))){
//						headerIndexes[j]=k;
//						continue;
//					}
//				}
//			}			
//			while((line=br.readLine())!=null){
//				String[] l=line.split(",");
//				String row=l[0];
//				int indexRight=Integer.MIN_VALUE;
//				for(int j=0;j<right.length;j++)if(row.replaceAll("\"","").equals(right[j].getName().replaceAll(".dot","")))indexRight=j;
//				for(int j=1;j<l.length;j++){
//					int indexLeft=headerIndexes[j];
//					if(indexLeft!=Integer.MIN_VALUE&&indexRight!=Integer.MIN_VALUE){
//						scores[indexLeft][indexRight]=Double.parseDouble(l[j]);
//					}
//				}
//			}
//
//			bw.write("Query: "+align.replaceAll("dot","massbank")+";");
//			bw.newLine();
//			bw.write("DB: "+with.replaceAll("dot","massbank")+";");
//			bw.newLine();
//			bw.newLine();
//
//			for(File h:right)bw.write(","+h.getName().replaceAll(".dot", ""));
//			bw.newLine();
//			for(int k=0;k<scores.length;k++){
//				bw.write(left[k].getName().replaceAll(".dot", ""));
//				for(double s:scores[k])bw.write(","+s);
//				bw.newLine();
//			}
//			bw.close();
//			br.close();



			System.out.println("...finished");


			//			tmpDirQueryAndDB.delete();

		}
	}
}
