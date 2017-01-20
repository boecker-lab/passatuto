
usearcTan=FALSE;

#folder="D:\\metabo_tandem_ms\\fdr_db_search\\results\\EstimatingPi\\mixdist\\R\\";
#for(f in list.files(folder)){
#  source(paste(folder,f,sep=""));
#}
library(mixdist)
weibull <- function(x, k, lambda, scale){
  z <- k/lambda*((x/lambda)^(k-1))*exp(-(x/lambda)^k);
  return(scale*z);
}

weibullDens <- function(x, k, lambda, scale){
  z <- exp(-(x/lambda)^k);
  return(scale*z);
}

weibullDensCutOff <- function(x, k, lambda, scale, maxValue){
  z <- (weibullDens(x, k, lambda, 1)-weibullDens(maxValue, k, lambda, 1))/(weibullDens(0, k, lambda, 1)-weibullDens(maxValue, k, lambda, 1))
  return(scale*z);
}

for(compounds in c("Compounds","Selections")){
#for(compounds in c("Compounds")){
  for(comparisonMethod in c("MassBank","CosineDistance")){
  #for(comparisonMethod in c("MassBank")){
     for(dataset in c("APFilesOriginal-GPFiles","APFilesOriginal-GPTrees","APTreesOriginal-GPTrees","OPFilesOriginal-GPFiles","OPFilesOriginal-GPTrees","OPTreesOriginal-GPTrees")){
     #for(dataset in c("APFilesOriginal-GPFiles")){
       for(hs in c("DifferentMassesExcluded","DMEHighFDR")){
       #for(hs in c("DifferentMassesExcluded")){

          tryCatch({
            inputFile=paste0("U:\\MSBlast\\qValues_All",compounds,"RPP\\pos\\",comparisonMethod,"\\HitStatistic_",hs,"_",dataset,"RandomPeaks\\QValues_Selection0_HitStatistic_",hs,"_",dataset,"RandomPeaks_1.hitlist")
            
            #args<-commandArgs(TRUE)
            #inputFile=args[1]
            #inputFile="U:\\MSBlast\\qValues_AllCompoundsRPP\\pos\\MassBank\\HitStatistic_DifferentMassesExcluded_APFilesOriginal-GPTreesConditionalFast\\QValues_Selection0_HitStatistic_DifferentMassesExcluded_APFilesOriginal-GPTreesConditionalFast_1.hitlist";
            #inputFile="U:\\MSBlast\\qValues_AllCompoundsRPP\\pos\\CosineDistance\\HitStatistic_DifferentMassesExcluded_OPFilesOriginal-GPFilesRandomPeaks\\QValues_Selection0_HitStatistic_DifferentMassesExcluded_OPFilesOriginal-GPFilesRandomPeaks_1.hitlist";
            print(inputFile)
            
            if(file.exists(inputFile)){   
              
              AllDataFull=read.csv(inputFile, header=T, skip=1, sep="\t") 
              if(usearcTan)AllDataFull[,"score"]=atan(AllDataFull[,"score"])
              AllData=AllDataFull[which(AllDataFull[,"match.type"]!="DecoyMatch"),]
              
              AllData=AllData[rev(order(AllData["score"])),]
              
              
              # i=dim(AllData)[1]
              # s=AllData[i,"score"]
              # tmp=AllData[which(AllData[,"score"]==s),]
              # firstVektor=as.integer(rownames(tmp))
              # secondVektor=as.integer(as.vector(rownames(tmp[order(tmp[,"match.type"],decreasing = TRUE),])))
              # for(j in 1:length(firstVektor)){
              #   print(AllData[firstVektor[j],])
              # }
              
              #for(i in 1:dim(AllData)[1]){
              #  s=AllData[i,"score"]
              #  tmp=AllData[which(AllData[,"score"]==s),]
              #  firstVektor=order(as.numeric(rownames(tmp)))
              #  secondVektor=as.numeric(rownames(tmp[order(tmp[,"match.type"],decreasing = TRUE),]))
              #  AllData[firstVektor,]=AllData[secondVektor,]
              #}
              
              TP=AllData[which(AllData[,"match.type"]=="TruePositiveMatch"),"score"]
              FP=AllData[which(AllData[,"match.type"]=="FalsePositiveMatch"),"score"]
              
              all <- c(TP, FP)  
              
              breaks <- 30  
              
              his <- hist(all, breaks=breaks)  
              df <- data.frame(mid=his$mids, cou=his$counts)   
              
              guemea <- c(0.001,0.999)  
              guesig <- c(1, 1)  
              guedis <- "weibull"              
              
              fitpro <- mix(as.mixdata(df), mixparam(mu=guemea, sigma=guesig), dist=guedis)
              #head(df) 
              
              #length(FP)/(length(TP)+length(FP))
              #fitpro$parameters
              
              add=""
              if(usearcTan)add="ATan"
              outputBase=paste0("U:\\MSBlast\\qValues_PEP_All",compounds,"RPP\\pos\\",comparisonMethod,"\\HitStatistic_",hs,"_",dataset,"CombinedR",add,"\\")
              outputBaseResult=paste0("U:\\MSBlast\\qValues_PEP_All",compounds,"RPP\\pos\\",comparisonMethod,"\\HitStatistic_",hs,"_",dataset,"PEPR",add,"\\")
              dir.create(outputBase, showWarnings = FALSE)
              dir.create(outputBaseResult, showWarnings = FALSE)
              #pdf(file=paste0("outputBase,"qValues_",hs,"_",dataset,"PEP.pdf"))
              #plot(fitpro, main=paste0("Weibull fitting\n",dataset,"\n",comparisonMethod))
              #dev.off()
              
              wFP=weibullpar(fitpro$parameters[1,"mu"], fitpro$parameters[1,"sigma"])
              kFP=wFP[1,"shape"]
              lambdaFP=wFP[1,"scale"]
              scaleFP=weibullDens(0, kFP, lambdaFP, fitpro$parameters[1,"pi"])-weibullDens(max(AllData["score"]), kFP, lambdaFP, fitpro$parameters[1,"pi"])
              
              wTP=weibullpar(fitpro$parameters[2,"mu"], fitpro$parameters[2,"sigma"])
              kTP=wTP[1,"shape"]
              lambdaTP=wTP[1,"scale"]
              scaleTP=weibullDens(0, kTP, lambdaTP, fitpro$parameters[2,"pi"])-weibullDens(max(AllData["score"]), kTP, lambdaTP, fitpro$parameters[2,"pi"])
              
              scaleFP=scaleFP/(scaleFP+scaleTP)
              scaleTP=1-scaleFP
              
              #print(scaleFP)
              
              xs=seq(0, 1, 0.01)
              
              hisFP <- hist(FP, breaks=20)              
              hisTP <- hist(TP, breaks=20)
              countsFP=hisFP$counts/(sum(hisFP$counts)+sum(hisTP$counts))*20
              countsTP=hisTP$counts/(sum(hisFP$counts)+sum(hisTP$counts))*20
              hisFP$counts=countsFP
              hisTP$counts=countsTP
              
              #pdf(file=paste0(ouputBase,"qValues_",hs,"_",dataset,"PEP_TP.pdf"))
              #plot( hisTP, xlim=c(0,1), main=paste0("Weibull fitting for TPs\n",dataset,"\n",comparisonMethod))  # first histogram              
              #lines(xs, weibull(x=xs,k=kTP,lambda=lambdaTP, scale=scaleTP), col="green", lwd=3)
              #dev.off()
              
              #pdf(file=paste0(outputBase,"qValues_",hs,"_",dataset,"PEP_FP.pdf"))
              #plot(hisFP, xlim=c(0,1), main=paste0("Weibull fitting for FPs\n",dataset,"\n",comparisonMethod))
              #lines(xs, weibull(x=xs,k=kFP,lambda=lambdaFP, scale=scaleFP), col="red", lwd=3)              
              #dev.off()
              
              his <- hist(c(TP,FP), breaks=20)
              counts=his$counts/(sum(his$counts))*20
              his$counts=counts
              pdf(file=paste0(outputBase,"HistogramDistribution_",hs,"_",dataset,".pdf"))
              plot(his, xlim=c(0,1), main=paste0("Weibull fitting for TPs and FPs\n",dataset,"\n",comparisonMethod))
              lines(xs, weibull(x=xs,k=kFP,lambda=lambdaFP, scale=scaleFP), col="red", lwd=3)   
              lines(xs, weibull(x=xs,k=kTP,lambda=lambdaTP, scale=scaleTP), col="green", lwd=3)
              lines(xs, weibull(x=xs,k=kFP,lambda=lambdaFP, scale=scaleFP)+weibull(x=xs,k=kTP,lambda=lambdaTP, scale=scaleTP), col="black", lwd=2)
              dev.off()
              
              AllData["FDR_true"]=0
              AllData["FDR_est_dist"]=0
              AllData["FDR_est_decoy"]=0
              for(i in 1:dim(AllData)[1]){
                s=AllData[i,"score"];
                FP_loc=length(which(AllData[,"score"]>=s&AllData[,"match.type"]=="FalsePositiveMatch"));
                TP_loc=length(which(AllData[,"score"]>=s&AllData[,"match.type"]=="TruePositiveMatch"));
                AllData[i,"FDR_true"]=FP_loc/(FP_loc+TP_loc);
                
                DP=length(which(AllDataFull[,"score"]>=s&AllDataFull[,"match.type"]=="DecoyMatch"));
                AllData[i,"FDR_est_decoy"]=DP/(FP_loc+TP_loc)*length(FP)/(length(TP)+length(FP));
                
                FP_est=weibullDensCutOff(s,k=kFP,lambda=lambdaFP, scale=scaleFP, max(AllData["score"]+0.0001));
                TP_est=weibullDensCutOff(s,k=kTP,lambda=lambdaTP, scale=scaleTP, max(AllData["score"]+0.0001));
                AllData[i,"FDR_est_dist"]=FP_est/(FP_est+TP_est);
              }
              
              AllData["qValue_true"]=0
              AllData["qValue_est_dist"]=0
              AllData["qValue_est_decoy"]=0
              min_qValue_true=1;
              min_qValue_est_dist=1;
              min_qValue_est_decoy=1;
              for(i in dim(AllData)[1]:1){
                min_qValue_true=min(min_qValue_true,AllData[i,"FDR_true"]);
                AllData[i,"qValue_true"]=min_qValue_true;
                
                min_qValue_est_dist=min(min_qValue_est_dist,AllData[i,"FDR_est_dist"]);
                AllData[i,"qValue_est_dist"]=min_qValue_est_dist;
                
                min_qValue_est_decoy=min(min_qValue_est_decoy,AllData[i,"FDR_est_decoy"]);
                AllData[i,"qValue_est_decoy"]=min_qValue_est_decoy;
              }
              
              print("result")
              
              resultTmp=AllData[,c("qValue_true","qValue_est_dist")]
              
              x=unique(resultTmp[,"qValue_true"])
              result=matrix(nrow=length(x),ncol=2)
              result[,1]=x
              
              for(i in 1:length(x)){
                result[i,2]=mean(resultTmp[which(resultTmp[,1]==x[i]),2])
              }
              
              diff=mean((result[,1]-result[,2])^2)
              parametersTP=paste0("WeibullDistribution: lambda=",lambdaTP," alpha=",kTP," pi=",scaleTP);
              parametersFP=paste0("WeibullDistribution: lambda=",lambdaFP," alpha=",kFP," pi=",scaleFP);
              
              oF=paste0(outputBaseResult,"QValues_HitStatistic_",hs,"_",dataset,"PEPR",add,".qValueMeanAverage")
              write(paste0(diff, " ",parametersTP, " ",parametersFP,"\n","calculated qValue\testimated qValue"), oF)              
              write.table(result,oF, row.names=FALSE, sep = "\t", col.names=FALSE, append=TRUE)
              
              #m=max(max(max(AllData[,"FDR_true"])),max(AllData[,"FDR_est_dist"]), max(AllData[,"FDR_est_decoy"]))
              #plot(seq(0,m,0.1),seq(0,m, 0.1))
              #points(AllData[,"FDR_true"], AllData[,"FDR_est_dist"], col="blue")
              #points(AllData[,"FDR_true"], AllData[,"FDR_est_decoy"], col="red")
              
            }
          }, error = function(e) {
            print(paste0(inputFile, " failed ", e))
          }, finally = {
            
          })
        }
      }
    }
}

