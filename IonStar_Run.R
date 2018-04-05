##Generate the annotated frame list
db <- "IonStarPRIDE_database.sdb" ##File name of the SIEVE database
sp <- "IonStarPRIDE_spectrum report.csv" ##File name of the spectrum report 
col_filename <- 4 ##Column number for rawfile name 
col_scannum <- 17 ##Column number for MS2 scan number
col_framelist <- c(6,18) ##Column number for Protein accession number and Peptide sequence
framelist <- "IonStarPRIDE_frame.csv" ##File name of the annotated frame list (output1)
sampleid <-"IonStarPRIDE_sampleid.csv" ##File name of the sample list (output2)
source ("IonStar_FrameGen.R")

##Perform protein quantification
library("IonStarStat")
setwd("File directory")
rawfile <- "Raw_input_IonStar.csv"
condfile <- "Group list.csv"
raw <- read.csv(rawfile)
cond <- read.csv(condfile)
condition <- cond[match(colnames(raw)[-c(1:3)], cond[,1]),2]
pdata <- newProDataSet(proData=raw, condition=condition)
ndata <- pnormalize(pdata, summarize=TRUE, method="TIC")
cdata<-OutlierPeptideRM(ndata,condition,variance=0.7,critM1=1/3,critM2=1/4,ratio=TRUE)
cdata<-SharedPeptideRM(cdata) ##Optional removal of shared peptides
quan <- ProteinQuan(eset=cdata, method="sum")
wrote.csv(quan,"Protein quantitative result.csv")
write.csv(exprs(cdata),"Peptide quantitative result.csv")

##Remove frames with too many missing values for outlier rejection (when error occurs for OutlierPeptideRM)
raw$mv.num <- apply(raw[,4:ncol(raw)],1,function(x) sum(x<1))
mv.tol <- ceiling((ncol(raw)-3)*0.1)
raw_mv <- subset(raw, mv.num<=mv.tol)
raw <- raw_mv[,1:(ncol(raw_mv)-1)]
