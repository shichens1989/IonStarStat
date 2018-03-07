library("IonStarStat")
setwd("File directory")
rawfile <- "Raw_input_IonStar.csv"
condfile <- "Group list.txt"
raw <- read.csv(rawfile)
cond <- read.table(condfile, header=TRUE)
condition <- cond[match(colnames(raw)[-c(1:3)], cond[,1]),2]
pdata <- newProDataSet(proData=raw, condition=condition)
ndata <- pnormalize(pdata, summarize=TRUE, method="TIC")
cdata<-OutlierPeptideRM(ndata,condition,variance=0.7,critM1=1/3,critM2=1/4,ratio=TRUE)
cdata<-SharedPeptideRM(cdata)
quan <- ProteinQuan(eset=cdata, method="sum")
wrote.csv(quan,"Protein quantitative result.csv")
write.csv(exprs(cdata),"Peptide quantitative result.csv")

##Removing frames with too many missing values for outlier rejection
raw$mv.num <- apply(raw[,4:ncol(raw)],1,function(x) sum(x<1))
mv.tol <- ceiling((ncol(raw)-3)*0.1)
raw_mv <- subset(raw, mv.num<=mv.tol)
raw <- raw_mv[,1:(ncol(raw_mv)-1)]