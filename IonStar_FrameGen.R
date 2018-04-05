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
frame[,"RAWFILE"]<-toupper(sub(".RAW","",frame[,"RAWFILE"]))
nframe <- toupper(paste(frame[,"RAWFILE"], frame[,"ms2scans.ms2scan"], sep="_"))

spectrum<- read.csv(sp)
nspec <- toupper(paste(spectrum[,4], spectrum[,17], sep="_"))
idx <- match(nspec, nframe)
message1 <- paste("Total number of spectra with matching frames:",sum(!is.na(idx)))
message2 <- paste("Total number of spectra without matching frames:",sum(is.na(idx)))
print(message1)
print(message2)
FS <- cbind(spectrum[!is.na(idx),], frame[na.omit(idx),])
FS[,col_framelist[1]]<- sub("\\,.*", "", FS[,col_framelist[1]])
FS[,col_framelist[1]]<- toupper(sub("sp\\|", "", FS[,col_framelist[1]]))
FS[,col_framelist[1]]<- toupper(sub("\\|", ":", FS[,col_framelist[1]]))
FS[,col_framelist[2]]<- toupper(FS[,col_framelist[2]])
FS$ID<-paste(FS[,col_framelist[1]],FS[,col_framelist[2]],FS[,"FrameID"],sep="|")
FS2<-FS[!duplicated(FS$ID),]
row.names(FS2)<-FS2$ID
##########output data################################
sub_start <- ncol(spectrum)+12
sub_end <- ncol(spectrum)+11+length(unique(frame[,"RAWFILE"]))
sub_data<-FS2[,c(col_framelist,(ncol(spectrum)+6),sub_start:sub_end)]
write.csv(sub_data,framelist,row.names=FALSE)
print("The annotated frame list has been generated.")
sampleid_list <- colnames(sub_data[,4:ncol(sub_data)])
sampleid_mat <- as.data.frame(matrix(c(sampleid_list,rep("NA",times=length(sampleid_list))),ncol=2,nrow=length(sampleid_list)))
colnames(sampleid_mat) <- c("Rawfiles","GroupID")
write.csv(sampleid_mat,sampleid,row.names=FALSE)
print("The sample list has been generated.")
