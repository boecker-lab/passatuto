package de.unijena.bioinf.statistics;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Locale;

public class qValueSummaryMain {
	
	static File statisticsFolder=new File("D:\\metabo_tandem_ms\\MSBlast\\statistics\\posOld\\");

	public static void main(String args[]) throws Exception{		
		List<String> comparisonMethods=new ArrayList<String>();
		List<String> datasets=new ArrayList<String>();
		
		for(File f:statisticsFolder.listFiles()){
			if(f.isDirectory()){
				comparisonMethods.add(f.getName());
				for(File f2:f.listFiles()){
					if(f2.getName().contains("qValuesComparisonMultiple")&&f2.getName().endsWith("txt")){
						String dataset=f2.getName().replaceAll("qValuesComparisonMultiple_", "").replaceAll(".txt", "");
						if(!datasets.contains(dataset))datasets.add(dataset);
					}
				}
			}
		}
		
		Collections.sort(comparisonMethods);
		Collections.sort(datasets);
		
		double[][] values=new double[comparisonMethods.size()][datasets.size()];
		for(int i=0;i<values.length;i++){		
			for(int j=0;j<values[i].length;j++){
				values[i][j]=Double.NaN;
			}
		}
		
		for(File f:statisticsFolder.listFiles()){
			if(f.isDirectory()){
//				System.out.println(f.getName());
				int i=comparisonMethods.indexOf(f.getName());
				for(File f2:f.listFiles()){
					if(f2.getName().contains("qValuesComparisonMultiple")&&f2.getName().endsWith("txt")){
						String dataset=f2.getName().replaceAll("qValuesComparisonMultiple_", "").replaceAll(".txt", "");
//						System.out.println(dataset);
						int j=datasets.indexOf(dataset);
						double value=getAverageValue(f2);
						values[i][j]=value;
					}
				}
			}
		}
		
		Locale.setDefault(Locale.US);
		DecimalFormat df=new DecimalFormat("0.000000");
		
		for(String d:datasets){
			System.out.print("\t"+d);
		}
		System.out.println();
		for(int i=0;i<values.length;i++){
			System.out.print(comparisonMethods.get(i));
			for(int j=0;j<values[i].length;j++){
				System.out.print("\t"+df.format(values[i][j]));
			}
			System.out.println();
		}
		
	}
	
	public static double getAverageValue(File f) throws Exception{
		BufferedReader br=new BufferedReader(new FileReader(f));
		String line="";
		while((line=br.readLine())!=null){
			if(line.startsWith("average sum")){
				return Double.parseDouble(line.split("\t")[3]);
			}
		}
		return Double.NaN;
	}
	
}
