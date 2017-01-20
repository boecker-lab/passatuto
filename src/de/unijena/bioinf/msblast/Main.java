package de.unijena.bioinf.msblast;

public class Main {
	public static void main(String args[]) throws Exception{

		//		for(String db:new String[]{"agilent","massbankorbi","massbankqtof","massbank"}){
		//			for(String d:new String[]{"Trees","Files"}){
		//				for(String method:new String[]{"Original"}){
		//					MSBlastDecoyTrees.main(
		//							("-base U:\\MSBlast\\Data\\ "
		//									+ "-c pos "
		//									+ "-db "+db+" "
		//									+ "-d "+d+" "
		//									+ "-method "+method+" "
		//									+ "-ppm 10 "
		//									+ "-ae 2").split(" "));
		//				}
		//			}
		//		}

		//		for(int i=1;i<=10;i++){
		//			for(String db:new String[]{/*"agilent",*/"massbankorbi","massbankqtof"}){
		//				for(String d:new String[]{/*"Trees",*/"Files"}){
		//					for(String method:new String[]{"RandomPeaks","MixedSpectrum","ConditionalFast","Reroot","RandomTree"}){
		//						MSBlastDecoyTrees.main(
		//								("-base U:\\MSBlast\\Data\\ "
		//										+ "-c pos "
		//										+ "-db "+db+" "
		//										+ "-d "+d+" "
		//										+ "-method "+method+" "
		//										+ "-addFolder _"+i+" "
		//										+ "-ppm 10 "
		//										+ "-ae 2").split(" "));
		//					}
		//				}
		//			}
		//		}

				for(int i=1;i<=1;i++){
					for(String db:new String[]{"massbankqtof"}){
						for(String d:new String[]{"Trees"}){
							for(String method:new String[]{"Reroot"}){
								MSBlastDecoyTrees.main(
										("-base U:\\MSBlast\\Data\\ "
												+ "-c pos "
												+ "-db "+db+" "
												+ "-d "+d+" "
												+ "-method "+method+" "
												+ "-addFolder _"+i+" "
												+ "-ppm 10 "
												+ "-ae 2").split(" "));
							}
						}
					}
				}

	}
}