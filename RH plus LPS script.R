#if you have ever done this before on any other script, you can just open the package with the library function
install.packages("devtools", "pander", "tidyverse","imputeTS", "parallel")
BiocManager::install("Rbowtie2")
install.packages("BiocManager") #BiocManager is very useful for genetic computations and has multiple 'sub packages'
BiocManager::install("ShortRead", 'Rsubread', 'Rqc', 'QuasR')

#Once you have the packages downloaded on your computer, you can open them with the library function
library(pander)
library(Rqc)
library(QuasR)
library(devtools)
library (BiocManager)
library (ShortRead)
library (Rsubread)
library (tidyverse)
library(imputeTS)
library(parallel)
library(doSNOW)
library(Rbowtie2)


RHLPSqc=rqc(path = "~/Lima research/Lima/RH plus LPS", pattern = ".fastq.gz", n = 1000, openBrowser=FALSE) 
rqcCycleQualityBoxPlot(RHLPSqc)
rqcCycleBaseCallsLinePlot(RHLSPqc)
rqcReadFrequencyPlot(RHLSPqc)

RHLPSqc_report <- rqcReport(RHLPSqc, outdir = "~/Lima research/Lima/RH plus LPS", file = "RHLPSqc_report")
openFileInOS(RHLPSqc_report)



RHLPStrimmed <- list.files("~/Lima research/Lima/RH plus LPS")[grepl('Galaxy', list.files("~/Lima research/Lima/RH plus LPS"))]
RHLPSalign_txt <- cbind(RHLPStrimmed, paste0('Sample', rep(1:length(RHLPStrimmed))))
colnames(RHLPSalign_txt) <- c('FileName', 'SampleName')
write.table(RHLPSalign_txt, 'RHLPSalign_txt.txt', col.names = T, quote = F, row.names = F, sep = "\t")
RHLPStxt <- "RHLPSalign_txt.txt" 
cl <- makeCluster(4)
FinalRHLPS <- qAlign("RHLPSalign_txt.txt", "ToxoGenome.fasta", clObj = cl)
alignmentStats("RH_Plus_LPS_4_40bc6b13e3d.bam")


