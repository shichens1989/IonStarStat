setwd("File directory")
##File name of SIEVE file
db<-"SIEVE database.sdb"
##File name of Spectrum report
xls<-"Spectrum report.csv"

##Use Line 8~37 if Spectrum report is generated from Scaffold
library(XLConnect)
library(RSQLite)
drv <- dbDriver("SQLite")
con <- dbConnect(drv, dbname=db)
sql <- "select frames.frameid,(frames.timestart+frames.timestop)/2 as Time,(frames.mzstart+frames.mzstop)/2 as MZ, ms2scans.ms2scan, ms2scans.RawFileID, FrameAttrib_1.* from Frames 
INNER JOIN MS2Scans on MS2Scans.FrameID=Frames.frameid 
INNER JOIN FrameAttrib_1 ON FrameAttrib_1.FrameID=Frames.FrameID
order by Frames.Frameid"
frame <- dbGetQuery(con, sql)
RawFiles <- dbGetQuery(con, "select * from RawFiles")
dbDisconnect(con)
raw <- RawFiles[match(frame[,"ms2scans.RawFileID"], RawFiles[,"RFID"]), "RawFile"]
raw <- sub(".raw", "", raw)
frame <- cbind(frame, RAWFILE=raw)
startRow <- grep("Experiment name", readLines(xls, n=200)) - 1
Spectrum <- read.csv(xls, sep="\t", skip=startRow)
frame[,"RAWFILE"]<-toupper(sub(".RAW","",frame[,"RAWFILE"])) #DELETE .RAW
nframe <- toupper(paste(frame[,"RAWFILE"], frame[,"ms2scans.ms2scan"], sep="_"))
nSepc <- toupper(sub("-\\d*-\\d*-\\d*", "", Spectrum[,"Spectrum.name"]))
nSepc <- toupper(sub("-\\d*-\\d*", "", Spectrum[,"Spectrum.name"]))
idx <- match(nSepc, nframe)
paste(sum(is.na(idx)),"spectrum don't have matched frame.")
FS <- cbind(Spectrum[!is.na(idx),], frame[na.omit(idx),])
FS[,"Protein.accession.numbers"]<- toupper(sub("\\,.*", "", FS[,"Protein.accession.numbers"]))
FS$ID<-paste(FS[,"Protein.accession.numbers"],FS[,"Peptide.sequence"],FS[,"FrameID"],sep="|")
FS2<-FS[!duplicated(FS$ID),]
row.names(FS2)<-FS2$ID
##Subset output (containing Protein AC, Peptide sequence, Frame ID, and Frame intensity in each run)
sub_data<-FS[,"columns"]
write.csv(sub_data,"Annotated Frame List.csv",row.names=FALSE)
write.csv(colnames(sub_data),"Sample List.csv")

##Use Line 39~70 if Spectrum report is generated from IonStarSPG.R
library(XLConnect)
library(RSQLite)
drv <- dbDriver("SQLite")
con <- dbConnect(drv, dbname=db)
sql <- "select frames.frameid,(frames.timestart+frames.timestop)/2 as Time,(frames.mzstart+frames.mzstop)/2 as MZ, ms2scans.ms2scan, ms2scans.RawFileID, FrameAttrib_1.* from Frames 
INNER JOIN MS2Scans on MS2Scans.FrameID=Frames.frameid 
INNER JOIN FrameAttrib_1 ON FrameAttrib_1.FrameID=Frames.FrameID
order by Frames.Frameid"
frame <- dbGetQuery(con, sql)
RawFiles <- dbGetQuery(con, "select * from RawFiles")
dbDisconnect(con)
raw <- RawFiles[match(frame[,"ms2scans.RawFileID"], RawFiles[,"RFID"]), "RawFile"]
raw <- sub(".raw", "", raw)
frame <- cbind(frame, RAWFILE=raw)
frame[,"RAWFILE"]<-toupper(sub(".RAW","",frame[,"RAWFILE"])) #DELETE .RAW
nframe <- toupper(paste(frame[,"RAWFILE"], frame[,"ms2scans.ms2scan"], sep="_"))
Spectrum_MSGF<- read.csv(MSGF)
Spectrum_MSGF[,"SpecFile"]<-sub(".raw.mzXML","",Spectrum_MSGF[,"SpecFile"])
MSepc <- toupper(paste(Spectrum_MSGF[,"SpecFile"], Spectrum_MSGF[,"ScanNum"], sep="_"))
idx <- match(MSepc, nframe)
paste(sum(is.na(idx)),"spectrum don't have matched frame.")
FS <- cbind(Spectrum_MSGF[!is.na(idx),], frame[na.omit(idx),])
FS[,"Accession"]<- sub("\\,.*", "", FS[,"Accession"])
FS[,"Accession"]<- toupper(sub("sp\\|", "", FS[,"Accession"]))
FS[,"Accession"]<- toupper(sub("\\|", ":", FS[,"Accession"]))
FS$ID<-paste(FS[,"Accession"],FS[,"Peptide.Sequence"],FS[,"FrameID"],sep="|")
FS2<-FS[!duplicated(FS$ID),]
row.names(FS2)<-FS2$ID
##Subset output (containing Protein AC, Peptide sequence, Frame ID, and Frame intensity in each run)
sub_data<-FS2[,c("Columns")]
write.csv(sub_data,"Annotated Frame List.csv",row.names=FALSE)
write.csv(colnames(sub_data),"Sample List.csv")

##If spectrum report is generated from other software packages, manual modification of the script is needed.