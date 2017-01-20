folder MSBlast/Data:
contains the data

OriginalFiles: contains the original ms data in .ms format
OriginalTrees: contains the computed trees in .dot format

Agilent: Agilent dataset (2120 files)
GNPS: GNPS dataset (4139 files)
MassBank: complete MassBank dataset (1000 files) 
MassBankOrbi: MassBank dataset restricted to Orbitrap measurements (458 files)
MassBankQTof: MassBank dataset restricted to QTof measurements (542 files)


folder MSBlast/MSBlast:
contains the source files and the libraries

usage:

 1. generation of decoy database
de.unijena.bioinf.msblast.MSBlastDecoyTrees: -base <folder of the original files, "/Data"> -c pos -db [gnps|agilent|massbankorbi|massbankqtof|massbank] -d [Trees|Files] -method [Original|RandomPeaks|MixedSpectrum|Reroot|RandomTree|ConditionalFast] -addFolder <string added to folder for repeated generation of decoy spectra> -ppm 10 -ae 2

2. spectral comparison
de.unijena.bioinf.spectralcomparison.MSBlastSpectralComparison <output folder, "/searchResultsRPP"> -rpp -method [MassBank|CosineDistance|OberacherWithIntensities|OberacherWithIntensitiesAndMasses|TreeAlignment] -ppm 10 -ae 2 -postprocess Simple -ds1 [first dataset, same syntax as for generation of decoy trees] -ds2 [second dataset, same syntax as for generation of decoy trees]

3. statistics for target-decoy-approach (very confusing, better to rewrite)
de.unijena.bioinf.statistics.MainStatistics -base "MSBlast/" -addFolder RPP -ppm 10 -ae 2 -queries APFilesOriginal-GPFiles,APFilesOriginal-GPTrees,APTreesOriginal-GPTrees,OPFilesOriginal-GPFiles,OPFilesOriginal-GPTrees,OPTreesOriginal-GPTrees -searchMethodsStatisticsOriginal CosineDistance,MassBank,OberacherWithIntensities,OberacherWithIntensitiesAndMasses,TreeAlignment -searchMethods CosineDistance,MassBank,OberacherWithIntensities,OberacherWithIntensitiesAndMasses,TreeAlignment -decoyMethods Original,RandomPeaks,MixedSpectrum,ConditionalFast,Reroot,RandomTree -getHitLists -statisticsOriginal -getRankStatistic -getQValues -getQValuesAverage -getQValueSummary

4. statistics for emperirical bayes approach
de.unijena.bioinf.em.EMMain (not parameterized, only works if the given folder structure is used)