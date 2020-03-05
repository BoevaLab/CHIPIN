library(stringr)
library(rtracklayer)
library(ggplot2)
library(matrixStats)
library(viridis)
library(preprocessCore)
library(tidyr)
library(plotrix)
library(rlist)
library(RColorBrewer)
library(tiger)
library(pracma)
library(GenomicFeatures)
library(dplyr)
library(gridExtra)
library(biomaRt)
library(rminer)
library(ggpubr)


###
'%ni%' <- Negate('%in%')
###

#######
# Takes as input the originel bw, duplicates it and transformes the score using a and b from the regression
#######
correctbigwigscore <-function(bw, a, b){
  #####
  #Load bigwig
  #####
  bwNonNormalized=import.wig(format="BigWig", bw)
  #########
  #Normalize with regression coefficient a and b
  #Also remove the values that are <0
  #########
  bw1_corr=bwNonNormalized
  score(bw1_corr)=(as.numeric(unlist(score(bwNonNormalized)))-b)/(a*1.0)
  tmp=which(score(bw1_corr)<0)
  if (length(tmp)!=0){
    score(bw1_corr[tmp])=0
  }
  return(bw1_corr)
}



#######
#Take as input a gtf file of the exons for the current organism
#######
create_exons <- function(gtf_file){
  exons_mm10_Jan19 <- read.delim(gtf_file)
  #
  #Correct the name field (remove the number after the point at the end)
  exons_mm10_Jan19$name=str_split(exons_mm10_Jan19$name, "[.]", 2, simplify=TRUE)[,1]
  colnames(exons_mm10_Jan19)=c("ENSEMBL_ID", colnames(exons_mm10_Jan19)[2:8])
  exons_mm10_Jan19_geneName=merge(exons_mm10_Jan19, corres, by="ENSEMBL_ID", all.x=F, all.y=F)
  #

  exons_mm10_Jan19_geneName=exons_mm10_Jan19_geneName[which(exons_mm10_Jan19_geneName$geneName != ""),]

  #Compute exons length VERSION2
  exonsLength=c()
  geneName=c()
  for (k in 1:nrow(exons_mm10_Jan19_geneName)){
    # #for (k in 1:10){
    total=0
    if (as.character(unlist(exons_mm10_Jan19_geneName[k,]$geneName)) %ni% geneName){
      tmp=which(exons_mm10_Jan19_geneName$geneName == (exons_mm10_Jan19_geneName[k,]$geneName))
      tmp_dataframe=exons_mm10_Jan19_geneName[tmp,]
      start=c()
      stop=c()
      for (i in 1:length(tmp)){
        start=c(start, noquote(str_split(exons_mm10_Jan19_geneName[tmp[i],]$exonStarts, "[,]", exons_mm10_Jan19_geneName[tmp[i],]$exonCount+1, simplify=TRUE)))
        start=start[which(start != "")]
        stop=c(stop, noquote(str_split(exons_mm10_Jan19_geneName[tmp[i],]$exonEnds, "[,]", exons_mm10_Jan19_geneName[tmp[i],]$exonCount+1, simplify=TRUE)))
        stop=stop[which(stop != "")]
      }
      df_tmp=data.frame()
      df_tmp=data.frame(rep(exons_mm10_Jan19_geneName[tmp[1],]$chrom, length(start)), start, stop, rep(exons_mm10_Jan19_geneName[tmp[1],]$strand, length(start)), rep(as.character(unlist(exons_mm10_Jan19_geneName[tmp[1],]$geneName)), length(start)))
      colnames(df_tmp)=c("chrom", "start", "stop", "strand", "geneName")
      tmp_Gr=makeGRangesFromDataFrame(df_tmp, keep.extra.columns=TRUE, start.field="start",  end.field="stop")
      exonsLength=c(exonsLength, sum(width(reduce(tmp_Gr))))
      geneName=c(geneName, as.character(unlist(exons_mm10_Jan19_geneName[tmp[1],]$geneName)))
    }
  }

  A=data.frame(geneName, exonsLength)
  return(A)
  #write.table(A, "/home/caramelito/NAS_Public/data/annotations/Mouse/mm10/GENCODE/exons_mm10_Jan19_length.txt", quote=F, row.names=F, col.names=F, sep="\t")

}

#######
#RPKMvalues should have two columns: one with geneName the second with RPKM values.
#exonsLength should have two columns: one with geneName the second with gene length
#######
get_readCounts<-function(RPKMvalues, exonsLength){
  colnames(RPKMvalues)=c("geneName", "value")
  colnames(exonsLength)=c("geneName", "value")
  exonsnames=t(exonsLength$geneName)
  readCounts=c()
  geneName=c()
  for (k in 1:nrow(RPKMvalues)){
    if (as.character(unlist(RPKMvalues[k,]$geneName)) %in% exonsnames){
      tmp=which(as.character(unlist(exonsLength$geneName)) == as.character(unlist(RPKMvalues[k,]$geneName)))
      geneName=c(geneName, as.character(unlist(RPKMvalues[k,]$geneName)))
      readcount=(RPKMvalues[k,]$value*exonsLength[tmp,]$value)
      readCounts=c(readCounts, readcount)
    }
  }
  rep=data.frame(geneName, readCounts)
  colnames(rep)=c("geneName", "CPM")
  return(rep)
}


#######
#raw read counts values should have two columns: one with geneName the second with raw read count values.
#exonsLength should have two columns: one with geneName the second with gene length
#######
get_RPKM<-function(raw_read_count, exonsLength){
  colnames(raw_read_count)=c("geneName", "value")
  colnames(exonsLength)=c("geneName", "value")
  RPKM_final=c()
  geneName=c()

  for (k in 1:nrow(raw_read_count)){
    if (as.character(unlist(raw_read_count[k,]$geneName)) %in% as.character(unlist(exonsLength$geneName))){
      tmp=which(as.character(unlist(exonsLength$geneName)) == as.character(unlist(raw_read_count[k,]$geneName)))
      geneName=c(geneName, as.character(unlist(raw_read_count[k,]$geneName)))
      RPKM_val=ceiling(raw_read_count[k,]$value/exonsLength[tmp,]$value)*((10^9)/sum(raw_read_count$value))

      RPKM_final=c(RPKM_final, RPKM_val)

    }
  }

  rep=data.frame(geneName, RPKM_final)
  colnames(rep)=c("geneName", "RPKM")
  return(rep)
}



#######
#Use for build density profile: prepare the input for compute_matrixOne
#######
prepareRegionToPlot_strandInDataP <- function(dataP, radius, step){
  regionsToPlot=makeGRangesFromDataFrame(dataP, keep.extra.columns=TRUE, start.field="start", end.field="stop", strand.field = "strand")
  myListOfCentralRegionsToPlot=GenomicRangesList()
  centers=as.numeric(as.character(unlist(dataP$peakmax)))
  add=dataP$geneName
  myListOfCentralRegionsToPlot=GRanges(seqnames(regionsToPlot), IRanges(centers-radius,centers+radius)[1:length(regionsToPlot)])
  A=data.frame(myListOfCentralRegionsToPlot)

  tmp=cbind(A, regionsToPlot$geneName, strand(regionsToPlot))

  colnames(tmp)=c("seqnames", "start", "end", "width", "*", "geneName", "strand")

  rep=makeGRangesFromDataFrame(tmp, keep.extra.columns = TRUE, ignore.strand=FALSE, start.field="start", end.field="end", strand.field = "strand")
  return(rep)
}



#######
#Compute matrix for one sample
#step is usually 50
#######
compute_matrixOne <- function(Gr_sample, myListOfCentralRegionsToPlot, radius, step){
  cpt=0
  overlapsDenSample=findOverlaps(myListOfCentralRegionsToPlot, Gr_sample)
  matrixSample=matrix(0,nrow=length(myListOfCentralRegionsToPlot),ncol=(2*radius/step+1))

  skippedRegions=NULL
  for(i in 1:length(myListOfCentralRegionsToPlot)){
    denRegionSample=Gr_sample[subjectHits(overlapsDenSample[which(queryHits(overlapsDenSample)==i)])]

    if ((as.character(unlist(strand(myListOfCentralRegionsToPlot[i])))=="+") | (as.character(unlist(strand(myListOfCentralRegionsToPlot[i])))=="-")){
      if (as.character(unlist(strand(myListOfCentralRegionsToPlot[i])))=="+"){
        newScoresSample=as.numeric(as.character(unlist(denRegionSample$score)))
      }else{
        newScoresSample=rev(denRegionSample$score)
      }
      newScoresSample[newScoresSample<0]=0
      newScoresSample=newScoresSample
      if (length(newScoresSample)==2*radius/step+1) {
        matrixSample[i,]=newScoresSample
        cpt=cpt+1
      }else{
        skippedRegions=c(skippedRegions,i)
      }
    }else{
      skippedRegions=c(skippedRegions,i)
    }
  }
  return(list(matrixSample, skippedRegions))
}


#######
#get bigwig data
#######
getDatabw_woRemoveNoise <- function(bwDox0){
  bwDox0.normalized=import.wig(format="BigWig", bwDox0)
  bwDox0.normalized=GRanges(bwDox0.normalized)
  return(bwDox0.normalized)
}


#######
#take as input either RPKM or raw read count (put None for the one you don't give) and find genes that are constant between samples/conditions
#######
find_gene_not_moving <- function(RPKM, raw_read_count, sample_name, output_dir, exons_length, Allgenes, percentage){
  ######################################
  #Step1: Find genes that do not move
  ######################################
  cat("\n")
  cat("*****************************************")
  cat("\n")
  cat("Will find ConstantGenes")
  cat("\n")
  cat("*****************************************")
  cat("\n")

  ####
  # Use raw read count
  ####
  if (is.null(raw_read_count)==FALSE & is.null(RPKM)==TRUE){
    cat("\n")
    cat("Will use raw read count to find ConstantGenes")
    cat("\n")
    CPM_rep=read.table(raw_read_count, header=TRUE)
    n_samples=ncol(CPM_rep)-1
    #Transform to matrix
    CPMvalues_all=data.frame(CPM_rep[,1:ncol(CPM_rep)])
    CPMvalues_all=na.omit(CPMvalues_all)

    geneNameToGetThrough=unique(CPMvalues_all[,1])
    Df_matrix=data.frame()
    #Do the sum on the duplicates: we sum the values for each replicates in the case of several transcipts
    for (k in 1:length(geneNameToGetThrough)){
      Df_matrix=rbind(Df_matrix, colSums(CPMvalues_all[which(CPMvalues_all$geneName == as.character(unlist(geneNameToGetThrough[k]))),2:ncol(CPM_rep)]))
    }
    Df_matrix_NoDup=cbind(geneNameToGetThrough, Df_matrix)
    #Compute mean and standard deviation
    CPM_rep_matrix_mean=apply(Df_matrix_NoDup[,2:ncol(CPM_rep)], 1, mean)
    CPM_rep_matrix_SD=apply(Df_matrix_NoDup[,2:ncol(CPM_rep)], 1, sd)
    #Plot
    pdf(paste(output_dir, sample_name, "_CPMmeanVSsd.pdf", sep=""), width=5, height=5)
    plot(log2(CPM_rep_matrix_mean+1), log2(CPM_rep_matrix_SD+1), xlim=c(0, max(log2(CPM_rep_matrix_mean+1))+5), ylim=c(0, max(log2(CPM_rep_matrix_SD+1))+5), type="p")

    #Now we have to select genes that do not move
    namesCorres=Df_matrix_NoDup[,1]

    #Order genes according to their mean (increasing order)
    tmp_order=order(CPM_rep_matrix_mean)
    namesCorresOrdered=namesCorres[tmp_order]
    CPMvalues_mean_ordered=CPM_rep_matrix_mean[tmp_order]
    CPMvalues_SD_ordered=CPM_rep_matrix_SD[tmp_order]


    #Define the number of genes to extract from each bin: nbGenesToPick
    ntot=dim(Df_matrix_NoDup)[1]
    nNotMove=round(ntot*(percentage/100))
    nbGenesPerBins=floor(ntot/100)
    nbGenesToPick=round(nNotMove/100)

    genesNotMoving_sd=c()
    genesNotMoving_names=c()
    genesNotMoving_mean=c()
    #Through the 100 bins
    for (j in 1:100){
      start=nbGenesPerBins*(j-1)+1
      stop=(nbGenesPerBins*j)
      #Get genes with the smallest standard deviation
      tmp_val=order(CPMvalues_SD_ordered[start:stop])
      tmp_val=start+tmp_val-1
      tmp_names=as.character(unlist(namesCorresOrdered[tmp_val]))
      genesNotMoving_names=c(genesNotMoving_names, head(tmp_names, nbGenesToPick))
      genesNotMoving_sd=c(genesNotMoving_sd, head(as.numeric(unlist(CPMvalues_SD_ordered[tmp_val])), nbGenesToPick))
      genesNotMoving_mean=c(genesNotMoving_mean, head(as.numeric(unlist(CPMvalues_mean_ordered[tmp_val])), nbGenesToPick))
    }

    tmp_geneToColor=which(namesCorresOrdered %in% genesNotMoving_names)
    #Colors on the plot sd function of mean the genes that are not moving in red.
    CPMvalues_mean_colors=CPMvalues_mean_ordered[tmp_geneToColor]
    CPMvalues_SD_colors=CPMvalues_SD_ordered[tmp_geneToColor]
    points(log2(CPMvalues_mean_colors+1), log2(CPMvalues_SD_colors+1), col="red", pch=16)
    dev.off()
    #####
    #Use RPKM
    #####
  }else if (is.null(RPKM)==FALSE & is.null(raw_read_count)==TRUE){
    cat("\n")
    cat("Will use RPKM to find ConstantGenes")
    cat("\n")
    RPKM_rep=read.table(RPKM, header=FALSE)
    #Remove RPKM with a geneName equal to "--"
    RPKM_rep_NoDup=RPKM_rep[which(RPKM_rep$V1 != "--"),]
    n_samples=ncol(RPKM_rep)-1

    #Get read counts for each sample
    for (k in 1:(n_samples)){
      RPKM_current=data.frame(RPKM_rep[,1], RPKM_rep[,k+1])
      tmp=get_readCounts(RPKM_current, exons_length)
      if (k==1){
        D=data.frame(tmp)
      }else{
        D=data.frame(D, tmp$CPM)
      }
    }
    D_matrix_all=data.frame(D[,1:ncol(D)])
    D_matrix_all=na.omit(D_matrix_all)

    #Do the sum on the duplicates: we sum the values for each replicates in the case of several transcipts
    geneNameToGetThrough=unique(D_matrix_all$geneName)
    Df_matrix=data.frame()
    for (k in 1:length(geneNameToGetThrough)){
      Df_matrix=rbind(Df_matrix, colSums(D_matrix_all[which(D_matrix_all$geneName == as.character(unlist(geneNameToGetThrough[k]))),2:ncol(D_matrix_all)]))
    }
    Df_matrix_NoDup=cbind(geneNameToGetThrough, Df_matrix)

    #Compute standard deviation and mean

    D_matrix_mean=apply(Df_matrix_NoDup[, 2:ncol(Df_matrix_NoDup)], 1, mean)

    D_matrix_SD=apply(Df_matrix_NoDup[, 2:ncol(Df_matrix_NoDup)], 1, sd)

    pdf(paste(output_dir, sample_name, "CPMmeanVSsd.pdf", sep=""), width=5, height=5)
    plot(log2(D_matrix_mean+1), log2(D_matrix_SD+1), xlim=c(0, max(log2(D_matrix_mean+1))+5), ylim=c(0, max(log2(D_matrix_SD+1))+5), type="p")

    #Now we have to select genes that do not move
    namesCorres=Df_matrix_NoDup[,1]
    #Order genes according to their mean
    tmp_order=order(D_matrix_mean)
    namesCorresOrdered=namesCorres[tmp_order]
    RPMvalues_mean_ordered=D_matrix_mean[tmp_order]
    RPMvalues_SD_ordered=D_matrix_SD[tmp_order]

    ntot=dim(Df_matrix_NoDup)[1]
    nNotMove=round(ntot*(percentage/100))
    nbGenesPerBins=floor(ntot/100)
    nbGenesToPick=round(nNotMove/100)

    genesNotMoving=c()
    genesNotMoving_names=c()
    genesNotMoving_sd=c()
    genesNotMoving_mean=c()
    for (j in 1:100){
      start=nbGenesPerBins*(j-1)+1
      stop=(nbGenesPerBins*j)
      #Get the genes with the smallest standard deviation in each bin
      tmp_val=order(RPMvalues_SD_ordered[start:stop])
      tmp_val=start+tmp_val-1

      tmp_names=as.character(unlist(namesCorresOrdered[tmp_val]))

      genesNotMoving_names=c(genesNotMoving_names, head(tmp_names, nbGenesToPick))
      genesNotMoving_sd=c(genesNotMoving_sd, head(as.numeric(unlist(RPMvalues_SD_ordered[tmp_val])), nbGenesToPick))
      genesNotMoving_mean=c(genesNotMoving_mean, head(as.numeric(unlist(RPMvalues_mean_ordered[tmp_val])), nbGenesToPick))
    }
    tmp_geneToColor=which(namesCorres %in% genesNotMoving_names)
    #Colors on the plot sd function of mean the genes that are not moving in red.
    D_matrix_mean_colors=D_matrix_mean[tmp_geneToColor]
    D_matrix_SD_colors=D_matrix_SD[tmp_geneToColor]
    points(log2(D_matrix_mean_colors+1), log2(D_matrix_SD_colors+1), col="red", pch=16)
    dev.off()
    ######
    # Use all TSS
    ######
  }else if (is.null(raw_read_count)==TRUE & is.null(RPKM)==TRUE){
    cat("\n")
    cat("Will use all genes as ConstantGenes")
    cat("\n")
    write.table(Allgenes, paste(output_dir, "ConstantGenes.bed", sep=""), quote=F, sep="\t", row.names=F, col.names=F)
    return(paste(output_dir, "ConstantGenes.bed", sep=""))
  }
  Df_genesNotMoving_names=data.frame(genesNotMoving_names)
  colnames(Df_genesNotMoving_names)=c("geneName")
  colnames(Allgenes)=c("chrom", "start", "stop", "geneName", "score", "strand")




  Genes_NotMoving=merge(Df_genesNotMoving_names, Allgenes, by="geneName", all.x=F, all.y=F)
  Genes_NotMovingNEW=data.frame(Genes_NotMoving$chr, Genes_NotMoving$start, Genes_NotMoving$stop, Genes_NotMoving$geneName, "1", Genes_NotMoving$strand)
  colnames(Genes_NotMovingNEW)=c("chrom", "start", "stop", "geneName", "score", "strand")


  cat(nrow(Genes_NotMovingNEW))
  cat(" constantGenes have been determined.")
  cat("\n")

  write.table(Genes_NotMovingNEW, paste(output_dir, sample_name, "GenesNotMoving.bed", sep=""), quote=F, sep="\t", row.names=F, col.names=F)
  return(paste(output_dir, sample_name, "GenesNotMoving.bed", sep=""))
}



#######
#run deeptools using genes not moving determined by find_gene_not_moving or supplied by the user then perform linear normalization with non zero intercept
#######
linear_normalization <- function(genesNotMoving_file, fileWithPaths, beforeRegionStartLength, afterRegionStartLength, regionBodyLength, binSize, output_dir){
  cat("\n")
  cat("*****************************************")
  cat("\n")
  cat("Will perform linear normalization with non-zero intercept")
  cat("\n")
  cat("*****************************************")
  cat("\n")

  #Initialize the list that will contain the plots before_after
  plist=list()
  listrenorm=c()
  #Read the file with inputs
  paths=read.table(fileWithPaths, header=FALSE)
  nSamples=nrow(paths)

  if (nSamples>2){
    #Get the median sample
    valMediane=c()
    nameMediane=c()
    for (k in 1:nSamples){
      bw1=as.character(unlist(paths[[1]][k]))
      current_bw=getDatabw_woRemoveNoise(bw1)
      current_score=score(current_bw)
      valMediane=c(valMediane, mean(current_score))
      nameMediane=c(nameMediane, bw1)
    }
    median_val=median(valMediane)
    #get the sample which is the closest to the median value
    ref=which.min((valMediane-median_val)^2)
    cat("\n")
    cat("Sample of reference is: ")
    cat("\n")
    cat(nameMediane[ref])
    cat("\n")
    bw1=nameMediane[ref]
    sample_nametmp1=strsplit(bw1, "/")[[1]]
    sample_nametmp2=sample_nametmp1[length(sample_nametmp1)]
    sample_ref=paste(strsplit(sample_nametmp2, "[.]")[[1]][1:length(strsplit(sample_nametmp2, "[.]")[[1]])-1], collapse = ".")
  }else{
    ref=1
    bw1=as.character(unlist(paths[[1]][ref]))
    sample_nametmp1=strsplit(bw1, "/")[[1]]
    sample_nametmp2=sample_nametmp1[length(sample_nametmp1)]
    sample_ref=paste(strsplit(sample_nametmp2, "[.]")[[1]][1:length(strsplit(sample_nametmp2, "[.]")[[1]])-1], collapse = ".")
  }

  #Work on the reference sample
  sample_nametmp1=strsplit(bw1, "/")[[1]]
  sample_nametmp2=sample_nametmp1[length(sample_nametmp1)]
  sample_name=paste(strsplit(sample_nametmp2, "[.]")[[1]][1:length(strsplit(sample_nametmp2, "[.]")[[1]])-1], collapse = ".")
  output_name=paste(output_dir, sample_name, ".mat.gz", sep="")

  bw_reference=bw1
  listrenorm=c(listrenorm, bw_reference)
  #
  cat("\n")
  cat("compute deeptools for reference sample")
  cat("\n")
  system(paste("computeMatrix scale-regions -S", bw1, "-R", genesNotMoving_file, "--beforeRegionStartLength", beforeRegionStartLength, "--afterRegionStartLength", afterRegionStartLength, "--regionBodyLength", regionBodyLength, "--binSize", binSize, "-o", output_name))
  matrix <- read.delim(output_name, header=FALSE, comment.char="@")
  val=data.frame(matrix[,7:ncol(matrix)])
  moyC1=apply(val, 2, function(x){return(mean(na.rm=TRUE, as.numeric(as.character(unlist(x)))))})
  genesNotMoving=read.table(genesNotMoving_file, header=F)
  colnames(genesNotMoving)=c("chrom", "start", "stop", "geneName", "score", "strand")
  genesNotMoving2=data.frame(genesNotMoving, genesNotMoving$start)
  colnames(genesNotMoving2)=c(colnames(genesNotMoving), "peakmax")
  cpt_sample=0
  for (i in 1:nSamples){
    cpt_sample=cpt_sample+1
    if(i==ref) next
    #Compute matrix using computeMatrix from deeptools
    bw1=as.character(unlist(paths[[1]][i]))
    sample_nametmp1=strsplit(bw1, "/")[[1]]
    sample_nametmp2=sample_nametmp1[length(sample_nametmp1)]
    sample_name=paste(strsplit(sample_nametmp2, "[.]")[[1]][1:length(strsplit(sample_nametmp2, "[.]")[[1]])-1], collapse = ".")
    output_name=paste(output_dir, sample_name, ".mat.gz", sep="")
    cat("\n")
    cat("compute deeptools for sample: ")
    cat(sample_name)
    cat("\n")

    system(paste("computeMatrix scale-regions -S", bw1, "-R", genesNotMoving_file, "--beforeRegionStartLength", beforeRegionStartLength, "--afterRegionStartLength", afterRegionStartLength, "--regionBodyLength", regionBodyLength, "--binSize", binSize, "-o", output_name))
    cat("deeptools done")
    #Now all matrix are computed.

    #Define a reference to perform the linear regression with non zero intercept

    #Perform linear regression with non-zero intercept
    #Get matrix
    matrixCondition1 <- read.delim(output_name, header=FALSE, comment.char="@")
    val=data.frame(matrixCondition1[,7:ncol(matrixCondition1)])
    moyCsample=apply(val, 2, function(x){return(mean(na.rm=TRUE, as.numeric(as.character(unlist(x)))))})

    reg=lm(moyCsample ~ moyC1)
    pdf(paste(output_dir, sample_name, "_regression.pdf", sep=""), width=5, height=5)
    #par(mfrow=c(2,2))
    plot(moyC1, moyCsample, pch=3, xlim=c(0, max(moyC1)+1), ylim=c(0, max(moyCsample)+2), xlab=sample_ref, ylab=sample_name)
    abline(0, 1, col="blue")
    abline(lm(moyCsample ~ moyC1), col="red")
    dev.off()
    cat("\n")
    cat("Coefficents of Regression")
    cat("\n")
    print(lm(moyCsample ~ moyC1))
    #a and b are the coefficients we will to correct the signal

    a=as.numeric(unlist(reg$coefficients[[2]]))
    b=as.numeric(unlist(reg$coefficients[[1]]))


    #Correct the signal
    filename_tmp=str_split(bw1, "[.]")[[1]][1:length(str_split(bw1, "[.]")[[1]])-1]
    M=c()
    for (j in 1:length(filename_tmp)){
      M=c(M, paste(filename_tmp[j], sep=""))
    }
    filename=M

    sample_nametmp1=strsplit(bw1, "/")[[1]]
    sample_nametmp2=sample_nametmp1[length(sample_nametmp1)]
    sample_name=paste(strsplit(sample_nametmp2, "[.]")[[1]][1:length(strsplit(sample_nametmp2, "[.]")[[1]])-1], collapse = ".")

    bw1_corr=correctbigwigscore(bw1, a, b)
    #output the new bigwig
    export(bw1_corr, paste(output_dir, sample_name, ".renormReg.bw", sep=""))
    #Import bigwig before
    bw_before=getDatabw_woRemoveNoise(bw1)
    #Import bigwig reference
    bw_ref=getDatabw_woRemoveNoise(bw_reference)

    listrenorm=c(listrenorm, paste(output_dir, sample_name, ".renormReg.bw", sep=""))
  }
  write.table(listrenorm, paste(output_dir, "pathRenorm.txt", sep=""), sep="\t", row.names = F, col.names = F, quote=F)
  #ggsave(paste(output_dir, "profiles_samples.pdf", sep=""), gridExtra::marrangeGrob(grobs = plist, nrow=2, ncol=2))
  return(paste(output_dir, "pathRenorm.txt", sep=""))
}


#######
#run deeptools using genes not moving determined by find_gene_not_moving or supplied by the user then perform quantile normalization
#######
quantile_norm <- function(fileWithPaths, genesNotMoving_file, nGroup, output_folder, beforeRegionStartLength, afterRegionStartLength, regionBodyLength, binSize){
  cat("\n")
  cat("*****************************************")
  cat("\n")
  cat("Will perform quantile normalization")
  cat("\n")
  cat("*****************************************")
  cat("\n")
  cat("\n")
  plist=list()
  #Read the file with inputs
  paths=read.table(fileWithPaths, header=FALSE)
  nSamples=nrow(paths)
  dataToNormHUGE=data.frame()
  for (i in 1:nSamples){
    current_bw=as.character(unlist(paths[[1]][i]))
    sample_nametmp1=strsplit(current_bw, "/")[[1]]
    sample_nametmp2=sample_nametmp1[length(sample_nametmp1)]
    sample_name=paste(strsplit(sample_nametmp2, "[.]")[[1]][1:length(strsplit(sample_nametmp2, "[.]")[[1]])-1], collapse = ".")
    output_name=paste(output_folder, sample_name, ".mat.gz", sep="")
    cat("Compute deeptools for sample:")
    cat("\n")
    cat(sample_name)
    cat("\n")
    system(paste("computeMatrix scale-regions -S", current_bw, "-R", genesNotMoving_file, "--beforeRegionStartLength", beforeRegionStartLength, "--afterRegionStartLength", afterRegionStartLength, "--regionBodyLength", regionBodyLength, "--binSize", binSize, "-o", output_name))
    matrix <- read.delim(output_name, header=FALSE, comment.char="@")
    current_values=unlist(matrix[,7:ncol(matrix)], use.names=FALSE)
    current_values[which(is.nan(current_values)==TRUE)]=0
    if ((i!=1) & (length(current_values) != dim(dataToNormHUGE)[1])){
      if (length(current_values)<nrow(dataToNormHUGE)){
        lengthNeeded=dim(dataToNormHUGE)[1]
        current_values=c(current_values, rep(0, lengthNeeded-length(current_values)))
      }else{
        current_values=current_values[1:dim(dataToNormHUGE)[1]]
      }
    }
    if (i==1){
      dataToNormHUGE=data.frame(current_values)
    }else{
      dataToNormHUGE=cbind(dataToNormHUGE, current_values)
    }
  }
  #Do the sum across bins (result: each gene has a mean value)
  dataSum=apply(dataToNormHUGE, 1, sum)
  #Order the genes
  tmp_order=order(dataSum)
  dataToNormHUGE_ordered=dataToNormHUGE[tmp_order,]
  nlignes=dim(dataToNormHUGE_ordered)[1]
  if (nGroup>1){
    #Define the number of genes per group
    nValPerGroup=floor(nlignes/nGroup)
    dataToNormHUGE_ordered_grouped=matrix(nrow=nGroup, ncol=nSamples)
    for (k in 1:nGroup){
      if (k!=nGroup){
        for (j in 1:nSamples){
          current_val=dataToNormHUGE_ordered[((k-1)*nValPerGroup)+1:(k*nValPerGroup), j]
          #current_val=dataToNormHUGE_ordered[((k-1)*nValPerGroup)+1:(k*nValPerGroup), j]
          dataToNormHUGE_ordered_grouped[k,j]=mean(current_val, na.rm=TRUE)
        }
      }else{
        for (j in 1:nSamples){
          #current_val=dataToNormHUGE_ordered[((k-1)*nValPerGroup):nlignes, j]
          current_val=dataToNormHUGE_ordered[((k-1)*nValPerGroup)+1:nlignes, j]
          dataToNormHUGE_ordered_grouped[k,j]=mean(current_val, na.rm=TRUE)
        }
      }
    }
    #Perform the quantile normalization on the groups
    dataNormHUGE_ordered_grouped=normalize.quantiles(dataToNormHUGE_ordered_grouped, copy=TRUE)
    dataToNormHUGE=dataToNormHUGE_ordered_grouped
    dataNormHUGE=dataNormHUGE_ordered_grouped
  }else{
    #Perform the quantile normalization
    dataNormHUGE_ordered=normalize.quantiles(as.matrix(dataToNormHUGE), copy=TRUE)
    dataNormHUGE=dataNormHUGE_ordered

  }
  ##
  #Get the orginal bigwig
  ##
  table_nom=data.frame()
  for (s in 1:nSamples){
    #Get the original bigwig
    name_bw_old=as.character(unlist(paths[[1]][s]))
    bw_old=getDatabw_woRemoveNoise(name_bw_old)
    bw_old_toModify=bw_old
    #Learn the transformation to apply
    smoothingSpline_current=smooth.spline(dataToNormHUGE[,s], dataNormHUGE[,s], spar=0.01)
    smoothingSpline_current_predict=predict(smoothingSpline_current, bw_old_toModify$score, deriv = 0)$y

    #Apply the transformation
    bw_old_toModify$score=smoothingSpline_current_predict

    sample_nametmp1=strsplit(name_bw_old, "/")[[1]]
    sample_nametmp2=sample_nametmp1[length(sample_nametmp1)]
    sample_name=paste(strsplit(sample_nametmp2, "[.]")[[1]][1:length(strsplit(sample_nametmp2, "[.]")[[1]])-1], collapse = ".")
    table_nom=c(table_nom, paste(output_folder, "QuantileNorm_", sample_name, ".bw", sep=""))
    export(bw_old_toModify, paste(output_folder, "QuantileNorm_", sample_name, ".bw", sep=""))

  }
  output_paths=paste(output_folder, "pathRenorm.txt", sep="")
  write.table(table_nom, output_paths, row.names=F, col.names=F, quote=F, sep="\n")
  return(output_paths)
}




##################
######
#Plots
#x should be: x=seq(-radius, radius, by=step)
#color should be a vector with 3 values
##c("indianred4", "sandybrown", "steelblue4")
######
##################

#OK
plot_after_quantile<-function(fileWithPaths, output_folder, genesNotMoving, step, DF_after, histoneMark){
  NotMoving=read.table(genesNotMoving, header=F)
  colnames(NotMoving)=c("chrom", "start", "stop", "geneName", "score", "strand")
  peakmax=c()
  for (k in 1:nrow(NotMoving)){
    peakmax=c(peakmax, NotMoving[k,]$start)
  }
  NotMoving$peakmax=peakmax
  paths=read.table(fileWithPaths, header=FALSE)

  #Define the number of plots: There will be 5 samples per plot including one "reference" sample which will be present on each plot
  nSamples=nrow(paths)
  if (nSamples!=5){
    #nSamples_woPremierGraph=nSamples-5
    #nbGraph=(nSamples_woPremierGraph %/% 4)+2
    nbGraph=2+((nSamples-5)%/%4)
  }else{
    nbGraph=1
  }
  plist=list()
  if (nbGraph>1){
    for (i in 1:(nbGraph-1)){
      bw_current=c()
      if (i==1){
        #The five first samples
        bw_current=paths[[1]][((i-1)*5+2):(i*5)]
        ref="popo"
      }else{
        #The four next samples + the reference one (i e the first)
        bw_current=paths[[1]][(((i-1)*4)+2):((i*4)+1)]
        ref="keke"
      }
      legend=c()
      if (i==1){
        #Define the reference
        bw_ref=paths[[1]][1]
        sample_nametmp1_ref=strsplit(as.character(unlist(bw_ref)), "/")
        sample_nametmp2_ref=sample_nametmp1_ref[[1]][length(sample_nametmp1_ref[[1]])]
        sample_name_ref_tmp=strsplit(sample_nametmp2_ref, ".bw")[[1]]
        sample_name_ref_tpm=strsplit(sample_name_ref_tmp, "_")[[1]][2]
        #sample_name_ref=paste(sample_name_ref_tpm[1], sample_name_ref_tpm[2], sep="_")
        legend=c(legend, sample_name_ref_tpm)
      }
      for (k in 1:length(bw_current)){
        if (i!= nbGraph){
          if ((i != 1) & (k==1)){
            #The first sample in legend is always the reference
            legend=c(legend, sample_name_ref_tpm)
          }
          sample_nametmp1=strsplit(as.character(unlist(bw_current[[k]])), "/")
          sample_nametmp2=sample_nametmp1[[1]][length(sample_nametmp1[[1]])]

          sample_name=strsplit(sample_nametmp2, ".bw")[[1]]

          sample_name_final=strsplit(sample_name, "_")[[1]][2]
          legend=c(legend, sample_name_final)
          #sample_name_final=paste(sample_name_tmp[1], sample_name_tmp[2], sep="_")

          # if (k!=5 | i==1){
          #   legend=c(legend, sample_name)
          # }
        }

      }
      #Open bw
      if (i != 1){
        ref="keke"
        bw1=getDatabw_woRemoveNoise(as.character(unlist(bw_ref)))
        bw2=getDatabw_woRemoveNoise(as.character(unlist(bw_current[1])))
        bw3=getDatabw_woRemoveNoise(as.character(unlist(bw_current[2])))
        bw4=getDatabw_woRemoveNoise(as.character(unlist(bw_current[3])))
        bw5=getDatabw_woRemoveNoise(as.character(unlist(bw_current[4])))
      }else{
        ref="popo"
        bw1=getDatabw_woRemoveNoise(as.character(unlist(bw_ref)))
        bw2=getDatabw_woRemoveNoise(as.character(unlist(bw_current[1])))
        bw3=getDatabw_woRemoveNoise(as.character(unlist(bw_current[2])))
        bw4=getDatabw_woRemoveNoise(as.character(unlist(bw_current[3])))
        bw5=getDatabw_woRemoveNoise(as.character(unlist(bw_current[4])))
        # bw1=getDatabw_woRemoveNoise(as.character(unlist(bw_current[1])))
        # bw2=getDatabw_woRemoveNoise(as.character(unlist(bw_current[2])))
        # bw3=getDatabw_woRemoveNoise(as.character(unlist(bw_current[3])))
        # bw4=getDatabw_woRemoveNoise(as.character(unlist(bw_current[4])))
        # bw5=getDatabw_woRemoveNoise(as.character(unlist(bw_current[5])))
      }

      p=plot_before_afternorm_several(bw1, bw2, bw3, bw4, bw5, NotMoving, legend, seq(-4000, 4000, by=step), "After normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen", "lightsalmon2"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
      DF_after=list.append(fill_statsAfter(bw1, bw2, bw3, bw4, bw5, NotMoving, legend, seq(-4000, 4000, by=step), "After normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen", "lightsalmon2"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after, ref))
      plist=list.append(plist, p)
    }

    #Samples that left: those that were not multiple of 5
    ref="keke"

    nResteToPlot=nSamples-(4*(nbGraph-1))-1
    #nPlotted=5+(nbGraph-1)*4
    nPlotted=length(DF_after)
    while (nResteToPlot>0){
      listToPlot=c()
      legend=c()
      if ((nResteToPlot-4)<=0){
        listToPlot=c(listToPlot, getDatabw_woRemoveNoise(as.character(unlist(bw_ref))))
        legend=c(legend, sample_name_ref_tpm)
        for (j in 1:nResteToPlot){
          listToPlot=c(listToPlot, getDatabw_woRemoveNoise(as.character(unlist(paths[[1]][(nPlotted+j)]))))

          sample_nametmp1=strsplit(as.character(unlist(as.character(unlist(paths[[1]][(nPlotted+j)])))), "/")
          sample_nametmp2=sample_nametmp1[[1]][length(sample_nametmp1[[1]])]

          sample_name=strsplit(sample_nametmp2, ".bw")[[1]]
          sample_name_final=strsplit(sample_name, "_")[[1]][2]
          legend=c(legend, sample_name_final)
          print(legend)
        }
        #listToPlot=c(listToPlot, getDatabw_woRemoveNoise(as.character(unlist(paths[[1]][(nPlotted+1):nSamples]))))
      }else if (nResteToPlot-4>0){
        listToPlot=c(listToPlot, getDatabw_woRemoveNoise(as.character(unlist(bw_ref))))
        legend=c(legend, sample_name_ref_tpm)
        for (j in 1:nResteToPlot){
          listToPlot=c(listToPlot, getDatabw_woRemoveNoise(as.character(unlist(paths[[1]][(nPlotted+j)]))))

          sample_nametmp1=strsplit(as.character(unlist(as.character(unlist(paths[[1]][(nPlotted+j)])))), "/")
          sample_nametmp2=sample_nametmp1[[1]][length(sample_nametmp1[[1]])]

          sample_name=strsplit(sample_nametmp2, ".bw")[[1]]
          sample_name_final=strsplit(sample_name, "_")[[1]][2]
          #sample_name_final=paste(sample_name_tmp[1], sample_name_tmp[2], sep="_")
          legend=c(legend, sample_name_final)
        }
      }

      if (length(listToPlot)==2){
        p=plot_before_afternorm_2profiles(listToPlot[[1]], listToPlot[[2]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
        DF_after=list.append(fill_statsAfter_1profiles(listToPlot[[2]], NotMoving, legend[2], seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after))
      }else if(length(listToPlot)==3){
        p=plot_before_afternorm_3profiles(listToPlot[[1]], listToPlot[[2]], listToPlot[[3]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
        DF_after=list.append(fill_statsAfter_2profiles(listToPlot[[2]], listToPlot[[3]], NotMoving, legend[2:length(legend)], seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after))
      }else if(length(listToPlot)==4){
        p=plot_before_afternorm_4profiles(listToPlot[[1]], listToPlot[[2]], listToPlot[[3]], listToPlot[[4]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
        DF_after=list.append(fill_statsAfter_3profiles(listToPlot[[2]], listToPlot[[3]], listToPlot[[4]], NotMoving, legend[2:length(legend)], seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after))
      }else if(length(listToPlot)==5){
        p=plot_before_afternorm_5profiles(listToPlot[[1]], listToPlot[[2]], listToPlot[[3]], listToPlot[[4]], listToPlot[[5]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
        DF_after=list.append(fill_statsAfter_4profiles(listToPlot[[2]], listToPlot[[3]], listToPlot[[4]], listToPlot[[5]][2:length(legend)], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after))
      }
      plist=list.append(plist, p)
      nPlotted=nPlotted+4
      nResteToPlot=nResteToPlot-4
    }

    main=marrangeGrob(grobs=plist,ncol=2, nrow=2)
    ggsave(paste(output_folder, "After_Normalisation.pdf", sep=""), main, width=10, height=10)

    #ggsave(paste(output_folder, "Before_Normalisation.pdf", sep=""), gridExtra::marrangeGrob(grobs = plist, nrow=2, ncol=2))
  }else if (nbGraph==1){
    ref="popo"
    legend=c()
    bw_current=paths[[1]][1:nSamples]
    list_bw_open=list()
    for (k in 1:nSamples){
      list_bw_open=list.append(list_bw_open, getDatabw_woRemoveNoise(as.character(unlist(bw_current[k]))))
      sample_nametmp1=strsplit(as.character(unlist(bw_current[[k]])), "/")
      sample_nametmp2=sample_nametmp1[[1]][length(sample_nametmp1[[1]])]
      sample_name=strsplit(sample_nametmp2, ".bw")[[1]]
      sample_name_final=strsplit(sample_name, "_")[[1]][2]
      #sample_name_final=paste(sample_name_tmp[1], sample_name_tmp[2], sep="_")
      legend=c(legend, sample_name_final)
    }
    #Given the number of samples to plot we choose the appropriate function
    if (nSamples==2){
      pdf(paste(output_folder, "After_Normalisation.pdf", sep=""), width=5, height=5)
      plot_before_afternorm_2profiles(list_bw_open[[1]], list_bw_open[[2]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
      dev.off()
      DF_after=fill_statsAfter_2profiles(list_bw_open[[1]], list_bw_open[[2]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after)
    }else if(nSamples==3){
      pdf(paste(output_folder, "After_Normalisation.pdf", sep=""), width=5, height=5)
      plot_before_afternorm_3profiles(list_bw_open[[1]], list_bw_open[[2]], list_bw_open[[3]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
      dev.off()
      DF_after=fill_statsAfter_3profiles(list_bw_open[[1]], list_bw_open[[2]], list_bw_open[[3]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after)
    }else if(nSamples==4){
      pdf(paste(output_folder, "After_Normalisation.pdf", sep=""), width=5, height=5)
      plot_before_afternorm_4profiles(list_bw_open[[1]], list_bw_open[[2]], list_bw_open[[3]], list_bw_open[[4]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
      dev.off()
      DF_after=fill_statsAfter_4profiles(list_bw_open[[1]], list_bw_open[[2]], list_bw_open[[3]], list_bw_open[[4]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after)
    }else if(nSamples==5){
      pdf(paste(output_folder, "After_Normalisation.pdf", sep=""), width=5, height=5)
      plot_before_afternorm_5profiles(list_bw_open[[1]], list_bw_open[[2]], list_bw_open[[3]], list_bw_open[[4]],  list_bw_open[[5]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen", "mediumvioletred"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
      dev.off()
      DF_after=fill_statsAfter_5profiles(list_bw_open[[1]], list_bw_open[[2]], list_bw_open[[3]], list_bw_open[[4]], list_bw_open[[5]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen", "mediumvioletred"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after)
    }

  }
  return(DF_after)
}

#OK
plot_after_linear<-function(fileWithPaths, output_folder, genesNotMoving, step, DF_after, histoneMark){
  NotMoving=read.table(genesNotMoving, header=F)
  colnames(NotMoving)=c("chrom", "start", "stop", "geneName", "score", "strand")
  peakmax=c()
  for (k in 1:nrow(NotMoving)){
    peakmax=c(peakmax, NotMoving[k,]$start)
  }
  NotMoving$peakmax=peakmax
  paths=read.table(fileWithPaths, header=FALSE)

  #Define the number of plots: There will be 5 samples per plot including one "reference" sample which will be present on each plot
  nSamples=nrow(paths)
  if (nSamples!=5){
    #nSamples_woPremierGraph=nSamples-5
    #nbGraph=(nSamples_woPremierGraph %/% 4)+2
    nbGraph=2+((nSamples-5)%/%4)
  }else{
    nbGraph=1
  }
  plist=list()
  if (nbGraph>1){
    for (i in 1:(nbGraph-1)){
      bw_current=c()
      if (i==1){
        #The five first samples
        bw_current=paths[[1]][((i-1)*5+2):(i*5)]
        ref="popo"
      }else{
        #The four next samples + the reference one (i e the first)
        bw_current=paths[[1]][(((i-1)*4)+2):((i*4)+1)]
        ref="keke"
      }
      legend=c()
      if (i==1){
        #Define the reference
        bw_ref=paths[[1]][1]
        sample_nametmp1_ref=strsplit(as.character(unlist(bw_ref)), "/")
        sample_nametmp2_ref=sample_nametmp1_ref[[1]][length(sample_nametmp1_ref[[1]])]
        sample_name_ref_tmp=strsplit(sample_nametmp2_ref, ".bw")[[1]]
        sample_name_ref=sample_name_ref_tmp
        #sample_name_ref_tpm=strsplit(sample_name_ref_tmp, "_")[[1]][2]
        #sample_name_ref=paste(sample_name_ref_tpm[1], sample_name_ref_tpm[2], sep="_")
        legend=c(legend, sample_name_ref_tmp)
      }
      for (k in 1:length(bw_current)){
        if (i!= nbGraph){
          if ((i != 1) & (k==1)){
            #The first sample in legend is always the reference
            legend=c(legend, sample_name_ref)
          }
          sample_nametmp1=strsplit(as.character(unlist(bw_current[[k]])), "/")
          sample_nametmp2=sample_nametmp1[[1]][length(sample_nametmp1[[1]])]

          sample_name=strsplit(sample_nametmp2, ".bw")[[1]]
          legend=c(legend, sample_name)
          #sample_name_final=strsplit(sample_name, "_")[[1]][2]
          #sample_name_final=paste(sample_name_tmp[1], sample_name_tmp[2], sep="_")

          # if (k!=5 | i==1){
          #   legend=c(legend, sample_name)
          # }
        }

      }
      #Open bw
      if (i != 1){
        ref="keke"
        bw1=getDatabw_woRemoveNoise(as.character(unlist(bw_ref)))
        bw2=getDatabw_woRemoveNoise(as.character(unlist(bw_current[1])))
        bw3=getDatabw_woRemoveNoise(as.character(unlist(bw_current[2])))
        bw4=getDatabw_woRemoveNoise(as.character(unlist(bw_current[3])))
        bw5=getDatabw_woRemoveNoise(as.character(unlist(bw_current[4])))
      }else{
        ref="popo"
        bw1=getDatabw_woRemoveNoise(as.character(unlist(bw_ref)))
        bw2=getDatabw_woRemoveNoise(as.character(unlist(bw_current[1])))
        bw3=getDatabw_woRemoveNoise(as.character(unlist(bw_current[2])))
        bw4=getDatabw_woRemoveNoise(as.character(unlist(bw_current[3])))
        bw5=getDatabw_woRemoveNoise(as.character(unlist(bw_current[4])))
        # bw1=getDatabw_woRemoveNoise(as.character(unlist(bw_current[1])))
        # bw2=getDatabw_woRemoveNoise(as.character(unlist(bw_current[2])))
        # bw3=getDatabw_woRemoveNoise(as.character(unlist(bw_current[3])))
        # bw4=getDatabw_woRemoveNoise(as.character(unlist(bw_current[4])))
        # bw5=getDatabw_woRemoveNoise(as.character(unlist(bw_current[5])))
      }

      p=plot_before_afternorm_several(bw1, bw2, bw3, bw4, bw5, NotMoving, legend, seq(-4000, 4000, by=step), "After normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen", "lightsalmon2"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
      DF_after=list.append(fill_statsAfter(bw1, bw2, bw3, bw4, bw5, NotMoving, legend, seq(-4000, 4000, by=step), "After normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen", "lightsalmon2"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after, ref))
      plist=list.append(plist, p)
    }

    #Samples that left: those that were not multiple of 5
    ref="keke"

    nResteToPlot=nSamples-(4*(nbGraph-1))-1
    #nPlotted=5+(nbGraph-1)*4
    nPlotted=length(DF_after)
    while (nResteToPlot>0){
      listToPlot=c()
      legend=c()
      if ((nResteToPlot-4)<=0){
        listToPlot=c(listToPlot, getDatabw_woRemoveNoise(as.character(unlist(bw_ref))))
        legend=c(legend, sample_name_ref)
        for (j in 1:nResteToPlot){
          listToPlot=c(listToPlot, getDatabw_woRemoveNoise(as.character(unlist(paths[[1]][(nPlotted+j)]))))

          sample_nametmp1=strsplit(as.character(unlist(as.character(unlist(paths[[1]][(nPlotted+j)])))), "/")
          sample_nametmp2=sample_nametmp1[[1]][length(sample_nametmp1[[1]])]

          sample_name=strsplit(sample_nametmp2, ".bw")[[1]]
          legend=c(legend, sample_name)
          print(legend)
        }
        #listToPlot=c(listToPlot, getDatabw_woRemoveNoise(as.character(unlist(paths[[1]][(nPlotted+1):nSamples]))))
      }else if (nResteToPlot-4>0){
        listToPlot=c(listToPlot, getDatabw_woRemoveNoise(as.character(unlist(bw_ref))))
        legend=c(legend, sample_name_ref)
        for (j in 1:nResteToPlot){
          listToPlot=c(listToPlot, getDatabw_woRemoveNoise(as.character(unlist(paths[[1]][(nPlotted+j)]))))

          sample_nametmp1=strsplit(as.character(unlist(as.character(unlist(paths[[1]][(nPlotted+j)])))), "/")
          sample_nametmp2=sample_nametmp1[[1]][length(sample_nametmp1[[1]])]

          sample_name=strsplit(sample_nametmp2, ".bw")[[1]]
          legend=c(legend, sample_name)
        }
      }

      if (length(listToPlot)==2){
        p=plot_before_afternorm_2profiles(listToPlot[[1]], listToPlot[[2]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
        DF_after=list.append(fill_statsAfter_1profiles(listToPlot[[2]], NotMoving, legend[2], seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after))
      }else if(length(listToPlot)==3){
        p=plot_before_afternorm_3profiles(listToPlot[[1]], listToPlot[[2]], listToPlot[[3]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
        DF_after=list.append(fill_statsAfter_2profiles(listToPlot[[2]], listToPlot[[3]], NotMoving, legend[2:length(legend)], seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after))
      }else if(length(listToPlot)==4){
        p=plot_before_afternorm_4profiles(listToPlot[[1]], listToPlot[[2]], listToPlot[[3]], listToPlot[[4]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
        DF_after=list.append(fill_statsAfter_3profiles(listToPlot[[2]], listToPlot[[3]], listToPlot[[4]], NotMoving, legend[2:length(legend)], seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after))
      }else if(length(listToPlot)==5){
        p=plot_before_afternorm_5profiles(listToPlot[[1]], listToPlot[[2]], listToPlot[[3]], listToPlot[[4]], listToPlot[[5]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
        DF_after=list.append(fill_statsAfter_4profiles(listToPlot[[2]], listToPlot[[3]], listToPlot[[4]], listToPlot[[5]][2:length(legend)], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after))
      }
      plist=list.append(plist, p)
      nPlotted=nPlotted+4
      nResteToPlot=nResteToPlot-4
    }

    main=marrangeGrob(grobs=plist,ncol=2, nrow=2)
    ggsave(paste(output_folder, "After_Normalisation.pdf", sep=""), main, width=10, height=10)

    #ggsave(paste(output_folder, "Before_Normalisation.pdf", sep=""), gridExtra::marrangeGrob(grobs = plist, nrow=2, ncol=2))
  }else if (nbGraph==1){
    ref="popo"
    legend=c()
    bw_current=paths[[1]][1:nSamples]
    list_bw_open=list()
    for (k in 1:nSamples){
      list_bw_open=list.append(list_bw_open, getDatabw_woRemoveNoise(as.character(unlist(bw_current[k]))))
      sample_nametmp1=strsplit(as.character(unlist(bw_current[[k]])), "/")
      sample_nametmp2=sample_nametmp1[[1]][length(sample_nametmp1[[1]])]
      sample_name=strsplit(sample_nametmp2, ".bw")[[1]]
      #sample_name_final=strsplit(sample_name, "_")[[1]][2]
      #sample_name_final=paste(sample_name_tmp[1], sample_name_tmp[2], sep="_")
      legend=c(legend, sample_name)
    }
    #Given the number of samples to plot we choose the appropriate function
    if (nSamples==2){
      pdf(paste(output_folder, "After_Normalisation.pdf", sep=""), width=5, height=5)
      plot_before_afternorm_2profiles(list_bw_open[[1]], list_bw_open[[2]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
      dev.off()
      DF_after=fill_statsAfter_2profiles(list_bw_open[[1]], list_bw_open[[2]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after)
    }else if(nSamples==3){
      pdf(paste(output_folder, "After_Normalisation.pdf", sep=""), width=5, height=5)
      plot_before_afternorm_3profiles(list_bw_open[[1]], list_bw_open[[2]], list_bw_open[[3]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
      dev.off()
      DF_after=fill_statsAfter_3profiles(list_bw_open[[1]], list_bw_open[[2]], list_bw_open[[3]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after)
    }else if(nSamples==4){
      pdf(paste(output_folder, "After_Normalisation.pdf", sep=""), width=5, height=5)
      plot_before_afternorm_4profiles(list_bw_open[[1]], list_bw_open[[2]], list_bw_open[[3]], list_bw_open[[4]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
      dev.off()
      DF_after=fill_statsAfter_4profiles(list_bw_open[[1]], list_bw_open[[2]], list_bw_open[[3]], list_bw_open[[4]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after)
    }else if(nSamples==5){
      pdf(paste(output_folder, "After_Normalisation.pdf", sep=""), width=5, height=5)
      plot_before_afternorm_5profiles(list_bw_open[[1]], list_bw_open[[2]], list_bw_open[[3]], list_bw_open[[4]],  list_bw_open[[5]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen", "mediumvioletred"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
      dev.off()
      DF_after=fill_statsAfter_5profiles(list_bw_open[[1]], list_bw_open[[2]], list_bw_open[[3]], list_bw_open[[4]], list_bw_open[[5]], NotMoving, legend, seq(-4000, 4000, by=step), "After Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen", "mediumvioletred"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after)
    }

  }
  return(DF_after)
}

#OK also works for before linear
plot_before_quantile<-function(fileWithPaths, output_folder, genesNotMoving, step, DF_after, histoneMark){
  NotMoving=read.table(genesNotMoving, header=F)
  colnames(NotMoving)=c("chrom", "start", "stop", "geneName", "score", "strand")
  peakmax=c()
  for (k in 1:nrow(NotMoving)){
    peakmax=c(peakmax, NotMoving[k,]$start)
  }
  NotMoving$peakmax=peakmax
  paths=read.table(fileWithPaths, header=FALSE)

  #Define the number of plots: There will be 5 samples per plot including one "reference" sample which will be present on each plot
  nSamples=nrow(paths)
  if (nSamples!=5){
    #nSamples_woPremierGraph=nSamples-5
    #nbGraph=(nSamples_woPremierGraph %/% 4)+2
    nbGraph=2+((nSamples-5)%/%4)
  }else{
    nbGraph=1
  }
  plist=list()
  if (nbGraph>1){
    for (i in 1:(nbGraph-1)){
      bw_current=c()
      if (i==1){
        #The five first samples
        bw_current=paths[[1]][((i-1)*5+2):(i*5)]
        ref="popo"
      }else{
        #The four next samples + the reference one (i e the first)
        bw_current=paths[[1]][(((i-1)*4)+2):((i*4)+1)]
        ref="keke"
      }
      legend=c()
      if (i==1){
        #Define the reference
        bw_ref=paths[[1]][1]
        sample_nametmp1_ref=strsplit(as.character(unlist(bw_ref)), "/")
        sample_nametmp2_ref=sample_nametmp1_ref[[1]][length(sample_nametmp1_ref[[1]])]
        sample_name_ref_tmp=strsplit(sample_nametmp2_ref, ".bw")[[1]]
        sample_name_ref=sample_name_ref_tmp
        #sample_name_ref_tpm=strsplit(sample_name_ref_tmp, "_")[[1]][2]
        #sample_name_ref=paste(sample_name_ref_tpm[1], sample_name_ref_tpm[2], sep="_")
        legend=c(legend, sample_name_ref_tmp)
      }
      for (k in 1:length(bw_current)){
        if (i!= nbGraph){
          if ((i != 1) & (k==1)){
            #The first sample in legend is always the reference
            legend=c(legend, sample_name_ref)
          }
          sample_nametmp1=strsplit(as.character(unlist(bw_current[[k]])), "/")
          sample_nametmp2=sample_nametmp1[[1]][length(sample_nametmp1[[1]])]

          sample_name=strsplit(sample_nametmp2, ".bw")[[1]]
          legend=c(legend, sample_name)
          #sample_name_final=strsplit(sample_name, "_")[[1]][2]
          #sample_name_final=paste(sample_name_tmp[1], sample_name_tmp[2], sep="_")

          # if (k!=5 | i==1){
          #   legend=c(legend, sample_name)
          # }
        }

      }
      #Open bw
      if (i != 1){
        ref="keke"
        bw1=getDatabw_woRemoveNoise(as.character(unlist(bw_ref)))
        bw2=getDatabw_woRemoveNoise(as.character(unlist(bw_current[1])))
        bw3=getDatabw_woRemoveNoise(as.character(unlist(bw_current[2])))
        bw4=getDatabw_woRemoveNoise(as.character(unlist(bw_current[3])))
        bw5=getDatabw_woRemoveNoise(as.character(unlist(bw_current[4])))
      }else{
        ref="popo"
        bw1=getDatabw_woRemoveNoise(as.character(unlist(bw_ref)))
        bw2=getDatabw_woRemoveNoise(as.character(unlist(bw_current[1])))
        bw3=getDatabw_woRemoveNoise(as.character(unlist(bw_current[2])))
        bw4=getDatabw_woRemoveNoise(as.character(unlist(bw_current[3])))
        bw5=getDatabw_woRemoveNoise(as.character(unlist(bw_current[4])))
        # bw1=getDatabw_woRemoveNoise(as.character(unlist(bw_current[1])))
        # bw2=getDatabw_woRemoveNoise(as.character(unlist(bw_current[2])))
        # bw3=getDatabw_woRemoveNoise(as.character(unlist(bw_current[3])))
        # bw4=getDatabw_woRemoveNoise(as.character(unlist(bw_current[4])))
        # bw5=getDatabw_woRemoveNoise(as.character(unlist(bw_current[5])))
      }

      p=plot_before_afternorm_several(bw1, bw2, bw3, bw4, bw5, NotMoving, legend, seq(-4000, 4000, by=step), "Before normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen", "lightsalmon2"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
      DF_after=list.append(fill_statsAfter(bw1, bw2, bw3, bw4, bw5, NotMoving, legend, seq(-4000, 4000, by=step), "Before normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen", "lightsalmon2"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after, ref))
      plist=list.append(plist, p)
    }

    #Samples that left: those that were not multiple of 5
    ref="keke"

    nResteToPlot=nSamples-(4*(nbGraph-1))-1
    #nPlotted=5+(nbGraph-1)*4
    nPlotted=length(DF_after)
    while (nResteToPlot>0){
      listToPlot=c()
      legend=c()
      if ((nResteToPlot-4)<=0){
        listToPlot=c(listToPlot, getDatabw_woRemoveNoise(as.character(unlist(bw_ref))))
        legend=c(legend, sample_name_ref)
        for (j in 1:nResteToPlot){
          listToPlot=c(listToPlot, getDatabw_woRemoveNoise(as.character(unlist(paths[[1]][(nPlotted+j)]))))

          sample_nametmp1=strsplit(as.character(unlist(as.character(unlist(paths[[1]][(nPlotted+j)])))), "/")
          sample_nametmp2=sample_nametmp1[[1]][length(sample_nametmp1[[1]])]

          sample_name=strsplit(sample_nametmp2, ".bw")[[1]]
          legend=c(legend, sample_name)
          print(legend)
        }
        #listToPlot=c(listToPlot, getDatabw_woRemoveNoise(as.character(unlist(paths[[1]][(nPlotted+1):nSamples]))))
      }else if (nResteToPlot-4>0){
        listToPlot=c(listToPlot, getDatabw_woRemoveNoise(as.character(unlist(bw_ref))))
        legend=c(legend, sample_name_ref)
        for (j in 1:nResteToPlot){
          listToPlot=c(listToPlot, getDatabw_woRemoveNoise(as.character(unlist(paths[[1]][(nPlotted+j)]))))

          sample_nametmp1=strsplit(as.character(unlist(as.character(unlist(paths[[1]][(nPlotted+j)])))), "/")
          sample_nametmp2=sample_nametmp1[[1]][length(sample_nametmp1[[1]])]

          sample_name=strsplit(sample_nametmp2, ".bw")[[1]]
          legend=c(legend, sample_name)
        }
      }

      if (length(listToPlot)==2){
        p=plot_before_afternorm_2profiles(listToPlot[[1]], listToPlot[[2]], NotMoving, legend, seq(-4000, 4000, by=step), "Before Normalization", 4000, step, c("indianred4", "steelblue4"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
        DF_after=list.append(fill_statsAfter_1profiles(listToPlot[[2]], NotMoving, legend[2], seq(-4000, 4000, by=step), "Before Normalization", 4000, step, c("indianred4"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after))
      }else if(length(listToPlot)==3){
        p=plot_before_afternorm_3profiles(listToPlot[[1]], listToPlot[[2]], listToPlot[[3]], NotMoving, legend, seq(-4000, 4000, by=step), "Before Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
        DF_after=list.append(fill_statsAfter_2profiles(listToPlot[[2]], listToPlot[[3]], NotMoving, legend[2:length(legend)], seq(-4000, 4000, by=step), "Before Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after))
      }else if(length(listToPlot)==4){
        p=plot_before_afternorm_4profiles(listToPlot[[1]], listToPlot[[2]], listToPlot[[3]], listToPlot[[4]], NotMoving, legend, seq(-4000, 4000, by=step), "Before Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
        DF_after=list.append(fill_statsAfter_3profiles(listToPlot[[2]], listToPlot[[3]], listToPlot[[4]], NotMoving, legend[2:length(legend)], seq(-4000, 4000, by=step), "Before Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after))
      }else if(length(listToPlot)==5){
        p=plot_before_afternorm_5profiles(listToPlot[[1]], listToPlot[[2]], listToPlot[[3]], listToPlot[[4]], listToPlot[[5]], NotMoving, legend, seq(-4000, 4000, by=step), "Before Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
        DF_after=list.append(fill_statsAfter_4profiles(listToPlot[[2]], listToPlot[[3]], listToPlot[[4]], listToPlot[[5]][2:length(legend)], NotMoving, legend, seq(-4000, 4000, by=step), "Before Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after))
      }
      plist=list.append(plist, p)
      nPlotted=nPlotted+4
      nResteToPlot=nResteToPlot-4
    }

    main=marrangeGrob(grobs=plist,ncol=2, nrow=2)
    ggsave(paste(output_folder, "Before_Normalisation.pdf", sep=""), main, width=10, height=10)

    #ggsave(paste(output_folder, "Before_Normalisation.pdf", sep=""), gridExtra::marrangeGrob(grobs = plist, nrow=2, ncol=2))
  }else if (nbGraph==1){
    ref="popo"
    legend=c()
    bw_current=paths[[1]][1:nSamples]
    list_bw_open=list()
    for (k in 1:nSamples){
      list_bw_open=list.append(list_bw_open, getDatabw_woRemoveNoise(as.character(unlist(bw_current[k]))))
      sample_nametmp1=strsplit(as.character(unlist(bw_current[[k]])), "/")
      sample_nametmp2=sample_nametmp1[[1]][length(sample_nametmp1[[1]])]
      sample_name=strsplit(sample_nametmp2, ".bw")[[1]]
      #sample_name_final=strsplit(sample_name, "_")[[1]][2]
      #sample_name_final=paste(sample_name_tmp[1], sample_name_tmp[2], sep="_")
      legend=c(legend, sample_name)
    }
    #Given the number of samples to plot we choose the appropriate function
    if (nSamples==2){
      pdf(paste(output_folder, "Before_Normalisation.pdf", sep=""), width=5, height=5)
      plot_before_afternorm_2profiles(list_bw_open[[1]], list_bw_open[[2]], NotMoving, legend, seq(-4000, 4000, by=step), "Before Normalization", 4000, step, c("indianred4", "steelblue4"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
      dev.off()
      DF_after=fill_statsAfter_2profiles(list_bw_open[[1]], list_bw_open[[2]], NotMoving, legend, seq(-4000, 4000, by=step), "Before Normalization", 4000, step, c("indianred4", "steelblue4"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after)
    }else if(nSamples==3){
      pdf(paste(output_folder, "Before_Normalisation.pdf", sep=""), width=5, height=5)
      plot_before_afternorm_3profiles(list_bw_open[[1]], list_bw_open[[2]], list_bw_open[[3]], NotMoving, legend, seq(-4000, 4000, by=step), "Before Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
      dev.off()
      DF_after=fill_statsAfter_3profiles(list_bw_open[[1]], list_bw_open[[2]], list_bw_open[[3]], NotMoving, legend, seq(-4000, 4000, by=step), "Before Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after)
    }else if(nSamples==4){
      pdf(paste(output_folder, "Before_Normalisation.pdf", sep=""), width=5, height=5)
      plot_before_afternorm_4profiles(list_bw_open[[1]], list_bw_open[[2]], list_bw_open[[3]], list_bw_open[[4]], NotMoving, legend, seq(-4000, 4000, by=step), "Before Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
      dev.off()
      DF_after=fill_statsAfter_4profiles(list_bw_open[[1]], list_bw_open[[2]], list_bw_open[[3]], list_bw_open[[4]], NotMoving, legend, seq(-4000, 4000, by=step), "Before Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after)
    }else if(nSamples==5){
      pdf(paste(output_folder, "Before_Normalisation.pdf", sep=""), width=5, height=5)
      plot_before_afternorm_5profiles(list_bw_open[[1]], list_bw_open[[2]], list_bw_open[[3]], list_bw_open[[4]],  list_bw_open[[5]], NotMoving, legend, seq(-4000, 4000, by=step), "Before Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen", "mediumvioletred"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "))
      dev.off()
      DF_after=fill_statsAfter_5profiles(list_bw_open[[1]], list_bw_open[[2]], list_bw_open[[3]], list_bw_open[[4]], list_bw_open[[5]], NotMoving, legend, seq(-4000, 4000, by=step), "Before Normalization", 4000, step, c("indianred4", "steelblue4", "darkorchid3", "forestgreen", "mediumvioletred"), "Distance from TSS [bp]", paste("Average density of", histoneMark, sep=" "), DF_after)
    }

  }
  return(DF_after)
}





plot_ggplot_4curves <- function(matrixC1, matrixC2, matrixC3, matrixC4, legend, x, title, color, ylab, xlab, xlimite){
  G1=data.frame(x, apply(matrixC1, 2, mean), apply(matrixC2, 2, mean), apply(matrixC3, 2, mean), apply(matrixC4, 2, mean))
  colnames(G1)=c("abscisse", legend)
  plot_data=gather(G1, condition, valeur, legend, factor_key=TRUE)
  return(ggplot(plot_data) + geom_line(data=plot_data, aes(x=abscisse, y=valeur, group=condition, colour=condition )) +
           geom_ribbon(aes(ymin = plot_data$valeur - std.error(matrixC1), ymax = plot_data$valeur + std.error(matrixC1), x=plot_data$abscisse, fill=condition), alpha = 0.2, stat="identity") + scale_color_manual(values=color) +
           labs(x=xlab, y=ylab, title=title) + theme_bw() + ylim(0, xlimite))
}

plot_ggplot_3curves <- function(matrixC1, matrixC2, matrixC3, legend, x, title, color, ylab, xlab, xlimite){
  G1=data.frame(x, apply(matrixC1, 2, mean), apply(matrixC2, 2, mean), apply(matrixC3, 2, mean))
  colnames(G1)=c("abscisse", legend)
  plot_data=gather(G1, condition, valeur, legend, factor_key=TRUE)
  return(ggplot(plot_data) + geom_line(data=plot_data, aes(x=abscisse, y=valeur, group=condition, colour=condition )) +
           geom_ribbon(aes(ymin = plot_data$valeur - std.error(matrixC1), ymax = plot_data$valeur + std.error(matrixC1), x=plot_data$abscisse, fill=condition), alpha = 0.2, stat="identity") + scale_color_manual(values=color) +
           labs(x=xlab, y=ylab, title=title) + theme_bw() + ylim(0, xlimite))
}

plot_ggplot_2curves <- function(matrixC1, matrixC2, legend, x, title, color, ylab, xlab, xlimite){
  G1=data.frame(x, apply(matrixC1, 2, mean), apply(matrixC2, 2, mean))
  colnames(G1)=c("abscisse", legend)
  plot_data=gather(G1, condition, valeur, legend, factor_key=TRUE)
  return(ggplot(plot_data) + geom_line(data=plot_data, aes(x=abscisse, y=valeur, group=condition, colour=condition )) +
           geom_ribbon(aes(ymin = plot_data$valeur - std.error(matrixC1), ymax = plot_data$valeur + std.error(matrixC1), x=plot_data$abscisse, fill=condition), alpha = 0.2, stat="identity") + scale_color_manual(values=color) +
           labs(x=xlab, y=ylab, title=title) + theme_bw() + ylim(0, xlimite))
}

plot_ggplot_fivecurves <- function(matrixC1, matrixC2, matrixC3, matrixC4, matrixC5, legend, x, title, color, ylab, xlab, xlimite){
  G1=data.frame(x, apply(matrixC1, 2, mean), apply(matrixC2, 2, mean), apply(matrixC3, 2, mean), apply(matrixC4, 2, mean), apply(matrixC5, 2, mean))
  colnames(G1)=c("abscisse", legend)
  plot_data=gather(G1, condition, valeur, legend, factor_key=TRUE)
  return(ggplot(plot_data) + geom_line(data=plot_data, aes(x=abscisse, y=valeur, group=condition, colour=condition )) +
           geom_ribbon(aes(ymin = plot_data$valeur - std.error(matrixC1), ymax = plot_data$valeur + std.error(matrixC1), x=plot_data$abscisse, fill=condition), alpha = 0.2, stat="identity") + scale_color_manual(values=color) +
           labs(x=xlab, y=ylab, title=title) + theme_bw() + ylim(0, xlimite))
}


plot_before_afternorm_2profiles <- function(sample1, sample2, D_TSS, legend, x, title, radius, step, color, xlab, ylab){
  regionD_TSS=prepareRegionToPlot_strandInDataP(D_TSS, radius, step)
  repD_sample1=compute_matrixOne(sample1, regionD_TSS, radius, step)
  repD_sample2=compute_matrixOne(sample2, regionD_TSS, radius, step)

  valmax=max(apply(repD_sample1[[1]], 2, mean), apply(repD_sample2[[1]], 2, mean))+5
  print(plot_ggplot_2curves(repD_sample1[[1]], repD_sample2[[1]], legend, x, title, color, ylab, xlab, valmax))
}

plot_before_afternorm_3profiles <- function(sample1, sample2, sample3, D_TSS, legend, x, title, radius, step, color, xlab, ylab){

  regionD_TSS=prepareRegionToPlot_strandInDataP(D_TSS, radius, step)
  print(regionD_TSS)
  repD_sample1=compute_matrixOne(sample1, regionD_TSS, radius, step)
  repD_sample2=compute_matrixOne(sample2, regionD_TSS, radius, step)
  repD_sample3=compute_matrixOne(sample3, regionD_TSS, radius, step)


  valmax=max(apply(repD_sample1[[1]], 2, mean), apply(repD_sample2[[1]], 2, mean), apply(repD_sample3[[1]], 2, mean))+5

  print(plot_ggplot_3curves(repD_sample1[[1]], repD_sample2[[1]], repD_sample3[[1]], legend, x, title, color, ylab, xlab, valmax))
}

plot_before_afternorm_4profiles <- function(sample1, sample2, sample3, sample4, D_TSS, legend, x, title, radius, step, color, xlab, ylab){
  regionD_TSS=prepareRegionToPlot_strandInDataP(D_TSS, radius, step)
  repD_sample1=compute_matrixOne(sample1, regionD_TSS, radius, step)
  repD_sample2=compute_matrixOne(sample2, regionD_TSS, radius, step)
  repD_sample3=compute_matrixOne(sample3, regionD_TSS, radius, step)
  repD_sample4=compute_matrixOne(sample4, regionD_TSS, radius, step)


  valmax=max(apply(repD_sample1[[1]], 2, mean), apply(repD_sample2[[1]], 2, mean), apply(repD_sample3[[1]], 2, mean), apply(repD_sample4[[1]], 2, mean))+5
  print(plot_ggplot_4curves(repD_sample1[[1]], repD_sample2[[1]], repD_sample3[[1]], repD_sample4[[1]], legend, x, title, color, ylab, xlab, valmax))

}

plot_before_afternorm_5profiles <- function(sample1, sample2, sample3, sample4, sample5, D_TSS, legend, x, title, radius, step, color, xlab, ylab){
  regionD_TSS=prepareRegionToPlot_strandInDataP(D_TSS, radius, step)
  repD_sample1=compute_matrixOne(sample1, regionD_TSS, radius, step)
  repD_sample2=compute_matrixOne(sample2, regionD_TSS, radius, step)
  repD_sample3=compute_matrixOne(sample3, regionD_TSS, radius, step)
  repD_sample4=compute_matrixOne(sample4, regionD_TSS, radius, step)
  repD_sample5=compute_matrixOne(sample5, regionD_TSS, radius, step)

  valmax=max(apply(repD_sample1[[1]], 2, mean), apply(repD_sample2[[1]], 2, mean), apply(repD_sample3[[1]], 2, mean), apply(repD_sample4[[1]], 2, mean), apply(repD_sample5[[1]], 2, mean))+5

  print(plot_ggplot_fivecurves(repD_sample1[[1]], repD_sample2[[1]], repD_sample3[[1]], repD_sample4[[1]], repD_sample5[[1]], legend, x, title, color, ylab, xlab, valmax))
}


plot_before_afternorm_several <- function(sample1, sample2, sample3, sample4, sample5, D_TSS, legend, x, title, radius, step, color, xlab, ylab){
  regionD_TSS=prepareRegionToPlot_strandInDataP(D_TSS, radius, step)
  repD_sample1=compute_matrixOne(sample1, regionD_TSS, radius, step)
  repD_sample2=compute_matrixOne(sample2, regionD_TSS, radius, step)
  repD_sample3=compute_matrixOne(sample3, regionD_TSS, radius, step)
  repD_sample4=compute_matrixOne(sample4, regionD_TSS, radius, step)
  repD_sample5=compute_matrixOne(sample5, regionD_TSS, radius, step)
  valmax=max(apply(repD_sample1[[1]], 2, mean), apply(repD_sample2[[1]], 2, mean), apply(repD_sample3[[1]], 2, mean), apply(repD_sample4[[1]], 2, mean), apply(repD_sample5[[1]], 2, mean))+5
  print(plot_ggplot_fivecurves(repD_sample1[[1]], repD_sample2[[1]], repD_sample3[[1]], repD_sample4[[1]], repD_sample5[[1]], legend, x, title, color, ylab, xlab, valmax))
}


fill_statsAfter_2profiles <- function(sample1, sample2, D_TSS, legend, x, title, radius, step, color, xlab, ylab, DF_after){
  regionD_TSS=prepareRegionToPlot_strandInDataP(D_TSS, radius, step)
  repD_sample1=compute_matrixOne(sample1, regionD_TSS, radius, step)
  repD_sample2=compute_matrixOne(sample2, regionD_TSS, radius, step)


  repD_sample1_av=apply(repD_sample1[[1]], 2, mean)
  repD_sample2_av=apply(repD_sample2[[1]], 2, mean)

  DF_after[[legend[1]]]=repD_sample1_av
  DF_after[[legend[2]]]=repD_sample2_av

  return(DF_after)
}

fill_statsAfter_1profiles <- function(sample1, D_TSS, legend, x, title, radius, step, color, xlab, ylab, DF_after){
  regionD_TSS=prepareRegionToPlot_strandInDataP(D_TSS, radius, step)
  repD_sample1=compute_matrixOne(sample1, regionD_TSS, radius, step)


  repD_sample1_av=apply(repD_sample1[[1]], 2, mean)

  DF_after[[legend[1]]]=repD_sample1_av

  return(DF_after)
}


fill_statsAfter_3profiles <- function(sample1, sample2, sample3, D_TSS, legend, x, title, radius, step, color, xlab, ylab, DF_after){
  regionD_TSS=prepareRegionToPlot_strandInDataP(D_TSS, radius, step)
  repD_sample1=compute_matrixOne(sample1, regionD_TSS, radius, step)
  repD_sample2=compute_matrixOne(sample2, regionD_TSS, radius, step)
  repD_sample3=compute_matrixOne(sample3, regionD_TSS, radius, step)


  repD_sample1_av=apply(repD_sample1[[1]], 2, mean)
  repD_sample2_av=apply(repD_sample2[[1]], 2, mean)
  repD_sample3_av=apply(repD_sample3[[1]], 2, mean)

  DF_after[[legend[1]]]=repD_sample1_av
  DF_after[[legend[2]]]=repD_sample2_av
  DF_after[[legend[3]]]=repD_sample3_av
  return(DF_after)
}



fill_statsAfter_4profiles <- function(sample1, sample2, sample3, sample4, D_TSS, legend, x, title, radius, step, color, xlab, ylab, DF_after){
  regionD_TSS=prepareRegionToPlot_strandInDataP(D_TSS, radius, step)
  repD_sample1=compute_matrixOne(sample1, regionD_TSS, radius, step)
  repD_sample2=compute_matrixOne(sample2, regionD_TSS, radius, step)
  repD_sample3=compute_matrixOne(sample3, regionD_TSS, radius, step)
  repD_sample4=compute_matrixOne(sample4, regionD_TSS, radius, step)

  repD_sample1_av=apply(repD_sample1[[1]], 2, mean)
  repD_sample2_av=apply(repD_sample2[[1]], 2, mean)
  repD_sample3_av=apply(repD_sample3[[1]], 2, mean)
  repD_sample4_av=apply(repD_sample4[[1]], 2, mean)

  DF_after[[legend[1]]]=repD_sample1_av
  DF_after[[legend[2]]]=repD_sample2_av
  DF_after[[legend[3]]]=repD_sample3_av
  DF_after[[legend[4]]]=repD_sample4_av

  return(DF_after)
}

fill_statsAfter_5profiles <- function(sample1, sample2, sample3, sample4, sample5, D_TSS, legend, x, title, radius, step, color, xlab, ylab, DF_after){
  regionD_TSS=prepareRegionToPlot_strandInDataP(D_TSS, radius, step)
  repD_sample1=compute_matrixOne(sample1, regionD_TSS, radius, step)
  repD_sample2=compute_matrixOne(sample2, regionD_TSS, radius, step)
  repD_sample3=compute_matrixOne(sample3, regionD_TSS, radius, step)
  repD_sample4=compute_matrixOne(sample4, regionD_TSS, radius, step)
  repD_sample5=compute_matrixOne(sample5, regionD_TSS, radius, step)

  repD_sample1_av=apply(repD_sample1[[1]], 2, mean)
  repD_sample2_av=apply(repD_sample2[[1]], 2, mean)
  repD_sample3_av=apply(repD_sample3[[1]], 2, mean)
  repD_sample4_av=apply(repD_sample4[[1]], 2, mean)
  repD_sample5_av=apply(repD_sample5[[1]], 2, mean)

  DF_after[[legend[1]]]=repD_sample1_av
  DF_after[[legend[2]]]=repD_sample2_av
  DF_after[[legend[3]]]=repD_sample3_av
  DF_after[[legend[4]]]=repD_sample4_av
  DF_after[[legend[5]]]=repD_sample5_av

  return(DF_after)
}


#######
#Create cluster to plot density profiles.
#The expressionvalues should contain 2 columns: geneName and values
#######
build_clusters_expression <- function(Expressionvalues){
  #Perform kmeans to determine the three clusters according to gene expression: genes high expressed, medium expressed, low expressed
  colnames(Expressionvalues)=c("geneName", "values")
  Expressionvalues$values=log2(Expressionvalues$values+1)
  Clusters=kmeans(Expressionvalues$values, 3)
  Valuesclusters=c()

  for (i in 1:length(Clusters[1])){
    Valuesclusters=c(Valuesclusters, Clusters[1][i])
  }

  data_RPKM_clustered=data.frame(Expressionvalues, Valuesclusters)

  data_RPKM_cluster1=data_RPKM_clustered[which(data_RPKM_clustered$cluster == "1"),]
  data_RPKM_cluster2=data_RPKM_clustered[which(data_RPKM_clustered$cluster == "2"),]
  data_RPKM_cluster3=data_RPKM_clustered[which(data_RPKM_clustered$cluster == "3"),]
  cat("\n")
  cat("Number of genes per cluster: ")
  cat("\n")
  cat("Cluster 1")
  cat("\n")
  cat(nrow(data_RPKM_cluster1))
  cat("\n")
  cat("Cluster 2")
  cat("\n")
  cat(nrow(data_RPKM_cluster2))
  cat("\n")
  cat("Cluster 3")
  cat("\n")
  cat(nrow(data_RPKM_cluster3))
  cat("\n")

  mini=which.min(c(mean(data_RPKM_cluster1$values), mean(data_RPKM_cluster2$values), mean(data_RPKM_cluster3$values)))
  maxi=which.max(c(mean(data_RPKM_cluster1$values), mean(data_RPKM_cluster2$values), mean(data_RPKM_cluster3$values)))

  cat("The mean values per cluster are: ")
  cat("\n")
  cat("Cluster 1: ")
  cat("\n")
  cat(mean(data_RPKM_cluster1$values))
  cat("\n")
  cat("Cluster 2: ")
  cat("\n")
  cat(mean(data_RPKM_cluster2$values))
  cat("\n")
  cat("Cluster 3: ")
  cat("\n")
  cat(mean(data_RPKM_cluster3$values))
  cat("\n")

  if (mini==1 & maxi==2){
    G1=data_RPKM_cluster1
    G2=data_RPKM_cluster3
    G3=data_RPKM_cluster2
  }else if(mini==1 & maxi==3){
    G1=data_RPKM_cluster1
    G2=data_RPKM_cluster2
    G3=data_RPKM_cluster3
  }else if (mini==2 & maxi==3){
    G1=data_RPKM_cluster2
    G2=data_RPKM_cluster1
    G3=data_RPKM_cluster3
  }else if (mini==2 & maxi==1){
    G1=data_RPKM_cluster2
    G2=data_RPKM_cluster3
    G3=data_RPKM_cluster1
  }else if (mini==3 & maxi==1){
    G1=data_RPKM_cluster3
    G2=data_RPKM_cluster2
    G3=data_RPKM_cluster1
  }else if (mini==3 & maxi==2){
    G1=data_RPKM_cluster3
    G2=data_RPKM_cluster1
    G3=data_RPKM_cluster2
  }
  return(list(G1, G2, G3))
}


#######
#Plot three cruves for plot_expression
#######

plot_ggplot_threecurves <- function(matrixC1, matrixC2, matrixC3, x, title, color, ylab, xlab){
  G1=data.frame(x, apply(matrixC1, 2, mean), apply(matrixC2, 2, mean), apply(matrixC3, 2, mean))
  colnames(G1)=c("abscisse", "Before Normalization", "After Normalization", "Reference")
  plot_data=gather(G1, condition, valeur, c("Before Normalization", "After Normalization", "Reference"), factor_key=TRUE)
  return(ggplot(plot_data) + geom_line(data=plot_data, aes(x=abscisse, y=valeur, group=condition, colour=condition )) +
           geom_ribbon(aes(ymin = plot_data$valeur - std.error(matrixC1), ymax = plot_data$valeur + std.error(matrixC1), x=plot_data$abscisse, fill=condition), alpha = 0.2, stat="identity") + scale_color_manual(values=color) +
           labs(x=xlab, y=ylab, title=title) + theme_bw())
}

plot_ggplot_threecurves_expression <- function(matrixC1, matrixC2, matrixC3, x, title, color, ylab, xlab){
  G1=data.frame(x, apply(matrixC1, 2, mean), apply(matrixC2, 2, mean), apply(matrixC3, 2, mean))
  colnames(G1)=c("abscisse", "Low expressed", "Medium expressed", "High expressed")
  plot_data=gather(G1, condition, valeur, c("Low expressed", "Medium expressed", "High expressed"), factor_key=TRUE)
  return(ggplot(plot_data) + geom_line(data=plot_data, aes(x=abscisse, y=valeur, group=condition, colour=condition )) +
           geom_ribbon(aes(ymin = plot_data$valeur - std.error(matrixC1), ymax = plot_data$valeur + std.error(matrixC1), x=plot_data$abscisse, fill=condition), alpha = 0.2, stat="identity") + scale_color_manual(values=color) +
           labs(x=xlab, y=ylab, title=title) + theme_bw())
}




#######
#Create density profiles according to gene expression with the output of build_clusters_expression
#######
profiles_gene_expression <- function(clusters, bw, D_TSS, radius, step){
  D_TSS_C1=merge(D_TSS, clusters[[1]], by="geneName", all.x=F, all.y=F)
  regionD_TSS_C1=prepareRegionToPlot_strandInDataP(D_TSS_C1, radius, step)
  repD_TSS_C1=compute_matrixOne(bw, regionD_TSS_C1, radius, step)
  D_TSS_C2=merge(D_TSS, clusters[[2]], by="geneName", all.x=F, all.y=F)
  regionD_TSS_C2=prepareRegionToPlot_strandInDataP(D_TSS_C2, radius, step)
  repD_TSS_C2=compute_matrixOne(bw, regionD_TSS_C2, radius, step)
  D_TSS_C3=merge(D_TSS, clusters[[3]], by="geneName", all.x=F, all.y=F)
  regionD_TSS_C3=prepareRegionToPlot_strandInDataP(D_TSS_C3, radius, step)
  repD_TSS_C3=compute_matrixOne(bw, regionD_TSS_C3, radius, step)
  return(list(repD_TSS_C1[[1]], repD_TSS_C2[[1]], repD_TSS_C3[[1]]))
}



####
#Evaluate the quality of the normalization
###

fill_statsAfter <- function(sample1, sample2, sample3, sample4, sample5, D_TSS, legend, x, title, radius, step, color, xlab, ylab, DF_after, ref){
  regionD_TSS=prepareRegionToPlot_strandInDataP(D_TSS, radius, step)
  repD_sample1=compute_matrixOne(sample1, regionD_TSS, radius, step)
  repD_sample2=compute_matrixOne(sample2, regionD_TSS, radius, step)
  repD_sample3=compute_matrixOne(sample3, regionD_TSS, radius, step)
  repD_sample4=compute_matrixOne(sample4, regionD_TSS, radius, step)
  repD_sample5=compute_matrixOne(sample5, regionD_TSS, radius, step)

  repD_sample1_av=apply(repD_sample1[[1]], 2, mean)
  repD_sample2_av=apply(repD_sample2[[1]], 2, mean)
  repD_sample3_av=apply(repD_sample3[[1]], 2, mean)
  repD_sample4_av=apply(repD_sample4[[1]], 2, mean)
  repD_sample5_av=apply(repD_sample5[[1]], 2, mean)
  if (ref=="popo"){
    DF_after[[legend[1]]]=repD_sample1_av
    DF_after[[legend[2]]]=repD_sample2_av
    DF_after[[legend[3]]]=repD_sample3_av
    DF_after[[legend[4]]]=repD_sample4_av
    DF_after[[legend[5]]]=repD_sample5_av
  }else{
    DF_after[[legend[2]]]=repD_sample2_av
    DF_after[[legend[3]]]=repD_sample3_av
    DF_after[[legend[4]]]=repD_sample4_av
    DF_after[[legend[5]]]=repD_sample5_av
  }
  print("La longueur apres fill stat after est : ")
  print(length(DF_after))
  return(DF_after)
}

fill_statsBefore <- function(sample1, sample2, sample3, sample4, sample5, D_TSS, legend, x, title, radius, step, color, xlab, ylab, DF_before, ref){
  regionD_TSS=prepareRegionToPlot_strandInDataP(D_TSS, radius, step)
  repD_sample1=compute_matrixOne(sample1, regionD_TSS, radius, step)
  repD_sample2=compute_matrixOne(sample2, regionD_TSS, radius, step)
  repD_sample3=compute_matrixOne(sample3, regionD_TSS, radius, step)
  repD_sample4=compute_matrixOne(sample4, regionD_TSS, radius, step)
  repD_sample5=compute_matrixOne(sample5, regionD_TSS, radius, step)

  repD_sample1_av=apply(repD_sample1[[1]], 2, mean)
  repD_sample2_av=apply(repD_sample2[[1]], 2, mean)
  repD_sample3_av=apply(repD_sample3[[1]], 2, mean)
  repD_sample4_av=apply(repD_sample4[[1]], 2, mean)
  repD_sample5_av=apply(repD_sample5[[1]], 2, mean)
  if (ref=="popo"){
    DF_before[[legend[1]]]=repD_sample1_av
    DF_before[[legend[2]]]=repD_sample2_av
    DF_before[[legend[3]]]=repD_sample3_av
    DF_before[[legend[4]]]=repD_sample4_av
    DF_before[[legend[5]]]=repD_sample5_av
  }else{
    DF_before[[legend[1]]]=repD_sample1_av
    DF_before[[legend[2]]]=repD_sample2_av
    DF_before[[legend[3]]]=repD_sample3_av
    DF_before[[legend[4]]]=repD_sample4_av
  }
  return(DF_before)
}




######
#Compute stats
######

compute_stats_AreaUnderCurves <- function(values, nSamples){
  val=c()
  tmp=c()
  for (j in 1:3){
    if (j==1){
      start=1
      stop=60
      start_integral=-4000
      stop_integral=-1050
    }else if(j==2){
      start=61
      stop=101
      start_integral=-1000
      stop_integral=1000
    }else{
      start=102
      stop=161
      start_integral=1050
      stop_integral=4000
    }
    for (i in 1:nSamples){
      approx=approxfun(seq(start_integral, stop_integral, by=50), values[[i]][start:stop])
      val=c(val, integral(approx, start_integral, stop_integral))
    }
  }
  percentage_values=c()
  for (k in 1:3){
    tmp=c(tmp, which.max(val[((k*nSamples)-(nSamples-1)):(k*nSamples)]))
    #percentage_values=c()
    for (m in 1:nSamples){
      #Compute percentage
      tmp_val=val[((k*nSamples)-(nSamples-1)):(k*nSamples)]
      percentage_values=c(percentage_values, (tmp_val[m]*100)/tmp_val[tmp[k]])
    }
  }
  return(percentage_values)
}


computeStats<-function(fileBefore, fileAfter, nSamples){
  listBefore=list()
  listAfter=list()
  valsBefore=read.table(fileBefore)
  valsAfter=read.table(fileAfter)
  for (j in 1:nSamples){
    listBefore[[j]]=valsBefore[,j]
    listAfter[[j]]=valsAfter[,j]
  }
  cat("\n")
  cat("Before normalization")
  cat("\n")
  cat("Region going from -4kB to -1kB")
  cat("\n")
  cat(compute_stats_AreaUnderCurves(listBefore, nSamples)[1:nSamples])
  cat("\n")
  cat("Region going from -1kB to +1kB")
  cat("\n")
  cat(compute_stats_AreaUnderCurves(listBefore, nSamples)[(nSamples+1):(2*nSamples)])
  cat("\n")
  cat("Region going from +1kB to +4kB")
  cat("\n")
  cat(compute_stats_AreaUnderCurves(listBefore, nSamples)[((2*nSamples)+1):(3*nSamples)])
  cat("\n")
  cat("\n")
  cat("Mean difference on the three regions before normalization: ")
  cat("\n")
  R1B=compute_stats_AreaUnderCurves(listBefore, nSamples)[1:nSamples][order(compute_stats_AreaUnderCurves(listBefore, nSamples)[1:nSamples], decreasing=T)]
  R2B=compute_stats_AreaUnderCurves(listBefore, nSamples)[(nSamples+1):(2*nSamples)][order(compute_stats_AreaUnderCurves(listBefore, nSamples)[(nSamples+1):(2*nSamples)], decreasing=T)]
  R3B=compute_stats_AreaUnderCurves(listBefore, nSamples)[((2*nSamples)+1):(3*nSamples)][order(compute_stats_AreaUnderCurves(listBefore, nSamples)[((2*nSamples)+1):(3*nSamples)], decreasing=T)]
  ddBefore=cbind(R1B, R2B, R3B)
  moyBefore=mean(ddBefore[1,]-ddBefore[2:nrow(ddBefore),])
  cat(moyBefore)
  cat(" %")
  cat("\n")
  cat("\n")

  cat("After normalization")
  cat("\n")
  cat("Region going from -4kB to -1kB")
  cat("\n")
  cat(compute_stats_AreaUnderCurves(listAfter, nSamples)[1:nSamples])
  cat("\n")
  cat("Region going from -1kB to +1kB")
  cat("\n")
  cat(compute_stats_AreaUnderCurves(listAfter, nSamples)[(nSamples+1):(2*nSamples)])
  cat("\n")
  cat("Region going from +1kB to +4kB")
  cat("\n")
  cat(compute_stats_AreaUnderCurves(listAfter, nSamples)[((2*nSamples)+1):(3*nSamples)])
  cat("\n")
  cat("\n")
  cat("Mean difference on the three regions after normalization: ")
  cat("\n")
  R1=compute_stats_AreaUnderCurves(listAfter, nSamples)[1:nSamples][order(compute_stats_AreaUnderCurves(listAfter, nSamples)[1:nSamples], decreasing=T)]
  R2=compute_stats_AreaUnderCurves(listAfter, nSamples)[(nSamples+1):(2*nSamples)][order(compute_stats_AreaUnderCurves(listAfter, nSamples)[(nSamples+1):(2*nSamples)], decreasing=T)]
  R3=compute_stats_AreaUnderCurves(listAfter, nSamples)[((2*nSamples)+1):(3*nSamples)][order(compute_stats_AreaUnderCurves(listAfter, nSamples)[((2*nSamples)+1):(3*nSamples)], decreasing=T)]
  dd_After=cbind(R1, R2, R3)
  moyAfter=mean(dd_After[1,]-dd_After[2:nrow(dd_After),])
  cat(moyAfter)
  cat(" %")
  cat("\n")
  cat("\n")
  cat("The difference between the highest density curve and the others density curves decreases of: ")
  diff=moyBefore-moyAfter
  cat(diff)
  cat(" %")

}





####################################################
# MAIN FUNCTIONS THAT WILL BE RUN VALENTINA VERSION
####################################################

plot_expression <- function(RPKM, raw_read_count, fileWithPaths, output_dir, organism, histoneMark="ChIP-seq signal"){
  if (organism == "hg19"){
    data("A_hg19")
    data("TSS_hg19")
    data("allgenes_hg19")
    exonsLength=A_hg19
    D_TSS=TSS_hg19
    Allgenes=allgenes_hg19
  }else if(organism == "hg38"){
    data("A_hg38")
    data("TSS_hg38")
    data("allgenes_hg38")
    exonsLength=A_hg38
    D_TSS=TSS_hg38
    Allgenes=allgenes_hg38
  }else if(organism == "mm10"){
    data("A_mm10")
    data("TSS_mm10")
    data("allgenes_mm10")
    exonsLength=A_mm10
    D_TSS=TSS_mm10
    Allgenes=allgenes_mm10
  }else if(organism == "mm9"){
    data("A_mm9")
    data("TSS_mm9")
    data("allgenes_mm9")
    exonsLength=A_mm9
    D_TSS=TSS_mm9
    Allgenes=allgenes_mm9
  }

  cat("\n")
  cat("*****************************************")
  cat("\n")
  cat("Will plot density profiles according to gene expression")
  cat("\n")
  cat("*****************************************")
  cat("\n")
  #initialization
  radius=4000
  step=50
  x=seq(-radius, radius, by=step)
  paths=read.table(fileWithPaths, header=FALSE)
  plist_expression=list()
  nSamples=nrow(paths)
  #Check if the output_dir parameter is provided
  if (is.null(output_dir)==TRUE){
    return(cat("You should provide a path to an existing folder for the output files."))
  }else if(is.null(D_TSS)==TRUE){
    return(cat("You should provide a data frame containing TSS information for the organism of interest. The data frame should contain the following columns: chrom start stop peakmax strand geneName. Some data frame are previously computed in the package: TSS_mm10 and TSS_hg19."))
  }else if (is.null(fileWithPaths)==TRUE){
    return(cat("You should provide a path to a text file containing path to .bigWig files."))
  }
  #Check if at least RPKM or raw read count has been provided
  if (is.null(RPKM)==TRUE & is.null(raw_read_count)==TRUE){
    return(cat("You should provide RPKM or raw_read_count"))

  } else if (is.null(RPKM)==TRUE & is.null(raw_read_count)!=TRUE){
    #if raw read count provided -> will compute RPKM
    if (is.null(exonsLength)==TRUE){
      return(cat("You should provide exons length in order to compute RPKM from raw read counts"))
    }

    rawReadCount=read.table(raw_read_count, header=TRUE)
    #Get read counts for each sample
    for (k in 1:(nSamples)){

      rawreadcount_current=data.frame(rawReadCount[,1], rawReadCount[,k+1])

      tmp=get_RPKM(rawreadcount_current, exonsLength)

      if (k==1){
        D=tmp
      }else{
        D=data.frame(D, tmp$RPKM)
      }
    }
    RPKM_rep=D

  } else if (is.null(RPKM)==FALSE & is.null(raw_read_count)==TRUE){
    #If RPKM provided -> will use it directly

    RPKM_rep=read.table(RPKM, header=TRUE)

  }
  tmpb <- RPKM_rep[,2:ncol(RPKM_rep)]
  #rownames(tmpb) <- RPKM_rep$V1
  tmp_apply=apply(data.matrix(tmpb), 1, mean)
  tmp_apply=data.frame(RPKM_rep[,1], tmp_apply)
  colnames(tmp_apply)=c("geneName", "values")
  #Clustering
  rep_clusters=build_clusters_expression(tmp_apply)
  nSamples=nrow(paths)
  for (i in 1:nSamples){
    bw1=as.character(unlist(paths[[1]][i]))

    sample_nametmp1=strsplit(bw1, "/")[[1]]
    sample_nametmp2=sample_nametmp1[length(sample_nametmp1)]
    sample_name=paste(strsplit(sample_nametmp2, "[.]")[[1]][1:length(strsplit(sample_nametmp2, "[.]")[[1]])-1], collapse = ".")

    bw1_read=getDatabw_woRemoveNoise(bw1)

    rep_profiles=profiles_gene_expression(rep_clusters, bw1_read, D_TSS, radius, step)

    color=c("indianred4", "sandybrown", "steelblue4")
    #plot_ggplot_several(matrixC1, matrixC2, matrixC3, x, title, color)
    p=plot_ggplot_threecurves_expression(rep_profiles[[1]], rep_profiles[[2]], rep_profiles[[3]], x, sample_name, color, paste("Average density of", histoneMark, sep=" "), "Distance from TSS [bp]")
    #p=plot_ggplot_several(rep_profiles[[1]], rep_profiles[[2]], rep_profiles[[3]], x, title, color)
    plist_expression=list.append(plist_expression, p)
  }
  #pdf('expression2.pdf',width=5, height=5)
  main=marrangeGrob(grobs=plist_expression,ncol=2, nrow=2)
  ggsave(paste(output_dir, "expression.pdf", sep=""), main, width=10, height=10)


  #grid.arrange(grobs = plist_expression, ncol = 2, nrow=2) ## display plot
  #ggsave(file = paste(output_dir, "expression.pdf", sep=""), arrangeGrob(grobs = plist_expression, ncol = 2, nrow=2), width=10, height=10)  ## save plot
  #grid.arrange(grobs=plist_expression, ncol = 2, nrow=2)
  #dev.off()
  #ggsave(paste(output_dir, "expression.pdf", sep=""), gridExtra::marrangeGrob(grobs = plist_expression, nrow=2, ncol=2))
}


ChIPIN_normalize <- function(RPKM, raw_read_count, pathToFileGenesNotMoving=NULL, sample_name, output_dir, organism, fileWithPaths, beforeRegionStartLength=4000, afterRegionStartLength=4000, regionBodyLength=40000, binSize=10, expression_plot=FALSE, compute_stat=FALSE, percentage=10, typeNorm="linear", nGroup=20, histoneMark="ChIP-seq signal"){
  radius=4000
  step=50
  DF_before=list()
  DF_after=list()

  paths=read.table(fileWithPaths, header=FALSE)
  nSamples=nrow(paths)

  if (organism == "hg19"){
    data("A_hg19")
    data("TSS_hg19")
    data("allgenes_hg19")
    exonsLength=A_hg19
    D_TSS=TSS_hg19
    Allgenes=allgenes_hg19
  }else if(organism == "hg38"){
    data("A_hg38")
    data("TSS_hg38")
    data("allgenes_hg38")
    exonsLength=A_hg38
    D_TSS=TSS_hg38
    Allgenes=allgenes_hg38
  }else if(organism == "mm10"){
    data("A_mm10")
    data("TSS_mm10")
    data("allgenes_mm10")
    exonsLength=A_mm10
    D_TSS=TSS_mm10
    Allgenes=allgenes_mm10
  }else if(organism == "mm9"){
    data("A_mm9")
    data("TSS_mm9")
    data("allgenes_mm9")
    exonsLength=A_mm9
    D_TSS=TSS_mm9
    Allgenes=allgenes_mm9
  }
  if (is.null(pathToFileGenesNotMoving)==TRUE & (is.null(RPKM)==FALSE | is.null(raw_read_count)==FALSE) & is.null(output_dir)==FALSE){
    pathToFileGenesNotMoving=find_gene_not_moving(RPKM, raw_read_count, sample_name, output_dir, exonsLength, Allgenes, percentage)
  }else if (is.null(pathToFileGenesNotMoving)==TRUE & is.null(RPKM)==TRUE & is.null(raw_read_count)==TRUE){
    pathToFileGenesNotMoving=find_gene_not_moving(RPKM, raw_read_count, sample_name, output_dir, exonsLength, Allgenes, percentage)
  }else if (is.null(pathToFileGenesNotMoving)==FALSE){
    cat("\n")
    cat("Constant Genes are provided.")
    cat("\n")
  }else{
    cat("\n")
    cat("One or more argument is not well provided for the function that will determine Constant Genes.")
    cat("\n")
  }
  if (typeNorm=="linear" & is.null(fileWithPaths)==FALSE & is.null(output_dir)==FALSE & is.null(D_TSS)==FALSE){
    pathRenorm=linear_normalization(pathToFileGenesNotMoving, fileWithPaths, beforeRegionStartLength, afterRegionStartLength, regionBodyLength, binSize, output_dir)

    rep_stats_after=plot_after_linear(pathRenorm, output_dir, pathToFileGenesNotMoving, step, DF_after,histoneMark)

    rep_stats_before=plot_before_quantile(fileWithPaths, output_dir, pathToFileGenesNotMoving, step, DF_before, histoneMark)
  }else if (typeNorm=="quantile" & is.null(fileWithPaths)==FALSE & is.null(output_dir)==FALSE & is.null(D_TSS)==FALSE){

    pathRenorm=quantile_norm(fileWithPaths, pathToFileGenesNotMoving, nGroup, output_dir, beforeRegionStartLength, afterRegionStartLength, regionBodyLength, binSize)

    rep_stats_after=plot_after_quantile(pathRenorm, output_dir, pathToFileGenesNotMoving, step, DF_after, histoneMark)

    rep_stats_before=plot_before_quantile(fileWithPaths, output_dir, pathToFileGenesNotMoving, step, DF_before, histoneMark)
  }else{
    cat("One or more argument is not well provided for the normalization process")
  }
  #Write stats
  write.table(data.frame(rep_stats_before), paste(output_dir, "StatsBefore.txt", sep=""), quote=F, sep="\t", col.names=F, row.names=F)
  write.table(data.frame(rep_stats_after), paste(output_dir, "StatsAfter.txt", sep=""), quote=F, sep="\t", col.names=F, row.names=F)


  if (expression_plot==TRUE){

    x=seq(-radius, radius, by=step)
    plot_expression(RPKM, raw_read_count, fileWithPaths, output_dir, organism, histoneMark)
  }
  if (compute_stat==TRUE){
    cat("\n")
    cat("Will compute statistics before and after normalization in % of the highest density curve.")
    cat("\n")
    cat("The order of the printed percentage correspond to the order of the samples in the file with paths.")
    computeStats(paste(output_dir, "StatsBefore.txt", sep=""), paste(output_dir, "StatsAfter.txt",  sep=""), nSamples)
  }
}

