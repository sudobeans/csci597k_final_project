library(APAlyzer)
library(ggplot2)
library(repmis)
library(GenomicRanges)
library(GenomicAlignments)

#RNA-seq BAM files
flsall <- dir(getwd(),".bam")
flsall<-paste0(getwd(),'/',flsall)
names(flsall)<-gsub('.bam','',dir(getwd(),".bam"))

#Genomic reference
URL="https://github.com/RJWANGbioinfo/PAS_reference_RData/blob/master/"
file="mm10_REF.RData"

#Building 3?UTR and intronic PAS
source_data(paste0(URL,file,"?raw=True"))
refUTRraw=refUTRraw
dfIPAraw=dfIPA
dfLEraw=dfLE   
PASREF=REF4PAS(refUTRraw,dfIPAraw,dfLEraw)
UTRdbraw=PASREF$UTRdbraw
dfIPA=PASREF$dfIPA
dfLE=PASREF$dfLE   

#Building aUTR and cUTR references
refUTRraw=refUTRraw
UTRdbraw=REF3UTR(refUTRraw)

#Calculation of relative expression
DFUTRraw=PASEXP_3UTR(UTRdbraw, flsall)

#Significance analysis of APA events
sampleTable1 = data.frame(samplename = c(names(flsall)),
                          condition = c("Control","Flor","Control","Flor","Control","Flor","Control","Flor"))
test_3UTRmuti=APAdiff(sampleTable1,
                     DFUTRraw, 
                     conKET="Control",
                     trtKEY="Flor",
                     PAS='3UTR',
                     CUTreads=0,
                     p_adjust_methods="fdr",
                     MultiTest='unpaired t-test')
table(test_3UTRmuti$APAreg)
write.csv(test_3UTRmuti, 'C:/Users/Densmore/Desktop/TauKOUTR.csv')

#Visualization of results
APAVolcano(test_3UTRmuti, PAS='3UTR', Pcol = "pvalue", top=5)

