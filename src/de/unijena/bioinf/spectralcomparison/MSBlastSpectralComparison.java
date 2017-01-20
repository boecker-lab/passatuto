package de.unijena.bioinf.spectralcomparison;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

import de.unijena.bioinf.decoy.model.MassBank;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction.CHARGE;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction.DATA;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction.DATASET;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction.METHOD;
import de.unijena.bioinf.msblast.ParametersDecoySpectrumConstruction.POSTPROCESS;
import de.unijena.bioinf.spectralcomparison.ParametersSpectralComparison.COMPARISONMETHOD;
import de.unijena.bioinf.spectralcomparison.ParametersSpectralComparison.COMPARISONPOSTPROCESS;

//generating spectral alignments

public class MSBlastSpectralComparison {

	public static void main (String args[]) throws Exception{
		Locale.setDefault(Locale.US);

		ParametersSpectralComparison p=new ParametersSpectralComparison(args);
		
		if(!p.isCorrectCombination()){
			System.err.println("\""+Arrays.toString(args)+"\" is not a correct parameter combination.");
			return;
		}
		System.out.println(p.getOutputFileComparison());

		List<MassBank> left=new ArrayList<MassBank>();
		List<MassBank> right=new ArrayList<MassBank>();
		List<MassBank> both=new ArrayList<MassBank>();

		for(File f:p.p1.getOutputFolderMassbank().listFiles()){
			MassBank m=new MassBank(f);
			if(p.removePP)m.removeParentPeaks(p.ppm, p.ae);
			left.add(m);
			both.add(m);

		}
		System.out.println(p.p2.getOutputFolderMassbank().getAbsolutePath());
		for(File f:p.p2.getOutputFolderMassbank().listFiles()){
			MassBank m=new MassBank(f);
			if(p.removePP)m.removeParentPeaks(p.ppm, p.ae);
			right.add(m);
			both.add(m);
		}

		ComparisonMethod comparison=p.getComparisonMethod();
		comparison.setSpectra(left,  right);
		AlignmentMatrix matrix=comparison.getAlignmentMatrix();
		matrix.writeToCSV(p.getOutputFileComparison());
		System.out.println("... finished");
		
	}

}
