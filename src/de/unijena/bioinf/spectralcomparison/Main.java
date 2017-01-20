package de.unijena.bioinf.spectralcomparison;

public class Main {
	public static void main(String args[]) throws Exception{

		for (String db1:new String[]{"massbankorbi","massbankqtof","agilent"}){
			for (String db2:new String[]{"massbankorbi","massbankqtof","agilent"}){
				if(db1.equals(db2))continue;
				for (String methodSearch:new String[]{/*"CosineDistance","MassBank", "OberacherWithIntensities","OberacherWithIntensitiesAndMasses", */"TreeAlignment"}){			
					for (String d1:new String[]{"Trees"/*, "Files"*/}){
						for (String d2:new String[]{"Trees"/*, "Files"*/}){
							for (String methodDecoy:new String[]{"Original"}){
								MSBlastSpectralComparison.main(
										("-base U:\\MSBlast\\searchResults\\ "
												+ "-method "+methodSearch+" "
												+ "-ppm 10 "
												+ "-ae 2 "
												+ "-postprocess Simple "

									+ "-ds1 "
									+"-base U:\\MSBlast\\Data\\ "
									+ "-c pos "
									+ "-db "+db1+" "
									+ "-d "+d1+" "
									+ "-method Original "
									+ "-ppm 10 "
									+ "-ae 2 " 

									+ "-ds2 "
									+ "-base U:\\MSBlast\\Data\\ "
									+ "-c pos "
									+ "-db "+db2+" "
									+ "-d "+d2+" "
									+ "-method "+methodDecoy+" "
									+ "-ppm 10 "
									+ "-ae 2").split(" "));	
							}
						}
					}
				}
			}
		}
		
	}
}
