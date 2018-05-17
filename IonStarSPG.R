#Import and format the protein list and spectrum list for all files from IDPicker
#Generate an empty data frame for final output
#Generate the list of sample names
setwd("File Directory")
protein <- read.csv("Protein list.tsv",header=TRUE,sep="\t",fileEncoding="windows-1252")
peptide <- read.table("Peptide list.tsv",header=TRUE,sep="\t",fileEncoding="windows-1252")
spectrum <- read.table("Spectra list.tsv",header=TRUE,sep="\t",fileEncoding="windows-1252")
list <- strsplit(as.matrix(spectrum$Spectrum.Rank),"/")
require("plyr")
df <- ldply(list)
df <- df[,2:3]
colnames(df) <- c("Sample.Name","Scan.Num")
spectrum <- data.frame(spectrum,df)
spectrumrepo_total <- data.frame(matrix(, ncol=25))
colnames(spectrumrepo_total) <- c("SpecFile","SpecID","ScanNum","FragMethod","Precursor","IsotopeError","PrecursorError.ppm.","Charge","DeNovoScore",
                                  "MSGFScore","SpecEValue","EValue","QValue","PepQValue","Peptide.Sequence","Accession","Cluster","Count","Coverage",
                                  "Protein.Group","Distinct.Peptides","Distinct.Matches","Filtered.Spectra","Description","Accession.2")
samplename <- unique(spectrum$Sample.Name)
#Subsetting the sample list and use part of the files
#Create a dataframe with individual protein names and group assignment
protein_2 <- protein[1,]
protein_2$Accession.2[1] <- "000000"
for (j in 1:nrow(protein)) {
    if (grepl(",", as.character(protein$Accession[j]))==FALSE) {
        pro_row <- protein[j,]
        pro_row$Accession.2 <- protein[j,1]
        protein_2 <- rbind(protein_2, pro_row)
    } else if (grepl(",", as.character(protein$Accession[j]))==TRUE) {
        proteinname <- unlist(strsplit(as.character(protein$Accession[j]),","))
        pro_row <- protein[j,]
        pro_row <- do.call("rbind", replicate(length(proteinname),pro_row,simplify=FALSE))
        pro_row$Accession.2 <- matrix(proteinname, nrow=length(proteinname))
        protein_2 <- rbind(protein_2, pro_row)
    }
}
protein_2 <- protein_2[2:nrow(protein_2),]

peptide$Sequence <- apply(as.matrix(peptide$Sequence),2, function(x) gsub("\\d|\\W","",x))

#Combining information from both IDPicker and MS-GF+ to generate the spectrum report
for (i in 1:length(samplename)) {
  filename <- paste(samplename[i],".tsv", sep="")
  MSGF <- read.table(filename,header=FALSE,sep="\t",fileEncoding="windows-1252")
  colnames(MSGF) <- c("SpecFile","SpecID","ScanNum","FragMethod","Precursor","IsotopeError","PrecursorError(ppm)","Charge","Peptide","Protein",
                      "DeNovoScore","MSGFScore","SpecEValue","EValue","QValue","PepQValue")
  matchlist <- match(as.matrix(spectrum$Sample.Name), as.character(samplename[i]), nomatch="NA")
  spectrum_2 <- data.frame(spectrum, matchlist)
  spectrum_used <- subset(spectrum_2, matchlist=="1")
  specnum <- spectrum_used$Scan.Num
  matchrow <- pmatch(MSGF$ScanNum,specnum,nomatch="NA",duplicates.ok=TRUE)
  matchlist <- unlist(strsplit(as.character(matchrow), " "))
  MSGF$Matchlist <- matchlist
  MSGF_2 <- na.omit(MSGF)
  MSGF_2$Peptide.Sequence <- apply(as.matrix(MSGF_2$Peptide),2,function(x) gsub("\\d|\\W","",x))
  MSGF_2$Seq <- row.names(MSGF_2)
  MSGF_2 <- MSGF_2[!duplicated(MSGF_2),]
  MSGF_2$Peptide.Sequence <- apply(MSGF_2$Peptide.Sequence, 2, function(x) substr(as.character(x),2,nchar(as.character(x))-1))
  MSGF_2$pepmatch <- match(MSGF_2$Peptide.Sequence, peptide$Sequence)
  MSGF_2 <- na.omit(MSGF_2)
  promatch <- pmatch(as.matrix(MSGF_2$Protein), as.matrix(protein_2$Accession.2), nomatch="NA", duplicates.ok=TRUE)
  protein_2$Number <- 1:nrow(protein_2)
  MSGF_2 <- data.frame(MSGF_2,promatch)
  colnames(MSGF_2)[21] <- "Number"
  spectrumrepo <- merge(MSGF_2, protein_2, by.y="Number")
  spectrumrepo <- spectrumrepo[,c(2:9,12:17,19,22:31)]
  spectrumrepo_total <- rbind(spectrumrepo_total,spectrumrepo)
}
#Delete the first row containing "NA"s and output the spectrum report
spectrumrepo_total <- spectrumrepo_total[2:nrow(spectrumrepo_total),]
write.csv(spectrumrepo_total, "Spectrum report.csv", row.names=FALSE)