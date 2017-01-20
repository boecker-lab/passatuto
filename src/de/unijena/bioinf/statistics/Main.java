package de.unijena.bioinf.statistics;

public class Main {

	public static void main(String args[]) throws Exception{
		String a="-base U:\\MSBlast\\ -ppm 10 -ae 2 -addFolder RPP "
		+"-queries "
				+"APFilesOriginal-GPTrees "
		+"-searchMethodsStatisticsOriginal MassBank "
		+"-searchMethods MassBank "
		+"-decoyMethods Reroot "		
//		+"-getHitLists "
//		+"-statisticsOriginal "
////		+"-getRankStatistic "
////		+"-checkEquality "
//		+"-getQValues "
//		+"-getQValuesAverage "
		+"-getQValueSummary "
////		+"-getQValueSummaryPerQValue "
//		+"-getQValueSummaryPerQValueEstimated "
//		+"-getQValueSummaryPerPosition "
		;
		 
		
		MainStatistics.main(a.split(" "));
		
	}
}
