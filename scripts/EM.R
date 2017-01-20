outputFile="U:\\MSBlast\\qValues_PEP_AllCompoundsRPP\\tmp\\pos\\MassBank\\HitStatistic_DifferentMassesExcluded_APFilesOriginal-GPFilesCombined\\test.pdf"
inputFileData="U:\\MSBlast\\qValues_PEP_AllCompoundsRPP\\pos\\MassBank\\HitStatistic_DifferentMassesExcluded_APFilesOriginal-GPFilesCombined\\SampleData_DifferentMassesExcluded_APFilesOriginal-GPFiles.txt"
inputFileDistribution="U:\\MSBlast\\qValues_PEP_AllCompoundsRPP\\tmp\\pos\\MassBank\\HitStatistic_DifferentMassesExcluded_APFilesOriginal-GPFilesCombined\\Distribution_DifferentMassesExcluded_APFilesOriginal-GPFiles_GammaDistribution_WeibullDistribution.txt"
outputFileFP="U:\\MSBlast\\qValues_PEP_AllCompoundsRPP\\tmp\\pos\\MassBank\\HitStatistic_DifferentMassesExcluded_APFilesOriginal-GPFilesCombined\\test_FP.pdf"
outputFileTP="U:\\MSBlast\\qValues_PEP_AllCompoundsRPP\\tmp\\pos\\MassBank\\HitStatistic_DifferentMassesExcluded_APFilesOriginal-GPFilesCombined\\test_TP.pdf"

args<-commandArgs(TRUE)


if(length(args)!=0){ 
  outputFile=args[1]
  inputFileData=args[2]
  inputFileDistribution=args[3]
}

data=read.csv(inputFileData, header=F, sep="\t")

pdf(file=outputFile)
hist(data[,1],breaks=20, freq=FALSE, main="", xlab="score", xlim = c(min(min(data[,1]),0),max(max(data[,1]),1)))
distribution=read.csv(inputFileDistribution, header=F, sep="\t")
palette=colorRampPalette(c("red", "green"))(n = dim(distribution)[2]-1);
o=order(colMeans(distribution)[2:dim(distribution)[2]])
for(i in 2:dim(distribution)[2]){
  lines(distribution[,1],distribution[,i], col=palette[o[i-1]], lwd=3)
}

sum=rowSums(distribution[2:dim(distribution)[2]])
lines(distribution[,1],sum, col="black", lwd=2)

dev.off()


sep=(dim(data)[2]>1)
if(sep){
  FP=data[which(data[,2]=="FalsePositiveMatch"),]
  TP=data[which(data[,2]=="TruePositiveMatch"),]
  if(length(args)!=0){ 
    outputFileFP=args[4]
    outputFileTP=args[5]
  }

  pdf(file=outputFileFP)
  hist(FP[,1],breaks=20, freq=FALSE, main="", xlab="score", xlim = c(min(min(FP[,1]),0),max(max(FP[,1]),1)))
  lines(distribution[,1],distribution[,which(o==1)+1]/length(FP[,1])*(length(TP[,1])+length(FP[,1])), col=palette[1], lwd=3)
  dev.off()
  
  pdf(file=outputFileTP)
  hist(TP[,1],breaks=20, freq=FALSE, main="", xlab="score", xlim = c(min(min(TP[,1]),0),max(max(TP[,1]),1)))
  lines(distribution[,1],distribution[,which(o==2)+1]/length(TP[,1])*(length(TP[,1])+length(FP[,1])), col=palette[2], lwd=3)
  dev.off()
}

