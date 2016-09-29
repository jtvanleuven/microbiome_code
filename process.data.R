library(EnvStats)
library(ggplot2)
library(stringr)

setwd("C:/Users/JT/Dropbox/gates/jt/hmo_091816_mm-if/")

#read count data
all <- read.delim("Gates_03102015_genus_for_paper_noFormula_noDups.txt",stringsAsFactors = FALSE)
meta <- read.delim("Gates_Milk_IF_MF_genus_diet_HMO_for_mar.txt",stringsAsFactors=FALSE)

meta_short <- cbind(subset(meta,select = c("Participant_Time","SampleType","Subject","Time.x")),meta[,1071:1083])
meta_short <- meta_short[which(meta_short$SampleType=="Milk"),]
meta_short <- meta_short[!is.na(meta_short$Sum),]

#fix names and stuff
Milk <- all[which(all$SampleType=="Milk"),]
IF <- all[which(all$SampleType=="IF"),]
in_both <- intersect(Milk$Participant_Time, IF$Participant_Time)
#in_both <- intersect(Milk$Participant_Time, meta_short$Participant_Time)
IF <- IF[IF$Participant_Time %in% in_both,]
Milk <- Milk[Milk$Participant_Time %in% in_both,]
all <- all[all$Participant_Time %in% in_both,]
meta_short <- meta_short[meta_short$Participant_Time %in% in_both,]

#add metadata to milk
#Milk <- merge(Milk,meta_short,by="Participant_Time")

#make matrix of taxa names
oxy.taxa <- rbind(colnames(all[,8:1065]),rep("genus",length(colnames(all[,8:1065]))))
colnames(oxy.taxa) <- oxy.taxa[1,]
for(i in 1:ncol(oxy.taxa)){
  name <- str_match(oxy.taxa[1,i], "(.*)_(.*)")[, 2]
  level <- str_match(oxy.taxa[1,i], "(.*)_(.*)")[, 3]
  oxy.taxa[1,i] <- name
  oxy.taxa[2,i] <- level
}
oxy.taxa <- cbind(c("OUT.ID","Level"),oxy.taxa)
oxy.taxa[1,1] <- "Lowest.Taxa"
oxy.taxa <- data.frame(oxy.taxa)

#get matrices in correct format
colnames(Milk)[3] <- "Time.Point"
colnames(IF)[3] <- "Time.Point"
colnames(Milk)[6] <- "Sample"
colnames(IF)[6] <- "Sample"
colnames(Milk)[5] <- "Subject"
colnames(IF)[5] <- "Subject"
colnames(Milk)[1066] <- "Total.Reads"
colnames(IF)[1066] <- "Total.Reads"

oxy.full.X <- cbind(subset(Milk,select=c("Sample","Subject","Time.Point","Total.Reads")),Milk[,8:1065])
theAnimal <- oxy.full.X$Subject
oxy.full.Y <- cbind(subset(IF,select=c("Sample","Subject","Time.Point","Total.Reads")),IF[,8:1065])
oxy.full.X <- oxy.full.X[do.call("order",oxy.full.X[,c("Subject","Time.Point")]),]
oxy.full.Y <- oxy.full.Y[do.call("order",oxy.full.Y[,c("Subject","Time.Point")]),]
oxy.full.X <- oxy.full.X[, ! names(oxy.full.X) %in% c("Sample","Subject","SampleType","Time.Point")]
oxy.full.Y <- oxy.full.Y[, ! names(oxy.full.Y) %in% c("Sample","Subject","SampleType","Time.Point")]
animal.matrix<-as.matrix(model.matrix(~theAnimal-1))

otuReads.X <- colSums(oxy.full.X[,! names(oxy.full.X) %in% c("SampleType.x","SampleType.y", "Sum", "Participant_SampleType_Time", "Side", "Participant_Time", "SampleType")], na.rm = T)
otuReads.Y <- colSums(oxy.full.Y[,! names(oxy.full.Y) %in% c("SampleType.x","SampleType.y", "Sum", "Participant_SampleType_Time", "Side", "Participant_Time", "SampleType")], na.rm = T)

keepOTU.X <- names(otuReads.X)[otuReads.X > quantile(otuReads.X,0.9)]
keepOTU.Y <- names(otuReads.Y)[otuReads.Y > quantile(otuReads.Y,0.9)]

top10.Y <- oxy.full.Y[,keepOTU.Y]
top10.X <- oxy.full.X[,keepOTU.X]
summary(colSums(top10.Y[,! names(top10.Y) %in% c("SampleType", "Total.Reads")],na.rm=T))
hist(log(colSums(top10.Y[,! names(top10.Y) %in% c("SampleType", "Total.Reads")],na.rm=T)), xlab = "Log(OTU Reads)")

yNames <- setdiff(names(top10.Y), "Total.Reads")
xNames <- setdiff(names(top10.X), "Total.Reads")
oxy.full.Y <- top10.Y[1:nrow(top10.Y),]
oxy.full.X <- top10.X[1:nrow(top10.X),]
oxy.full.X[,xNames] <- oxy.full.X[,xNames] / oxy.full.X[,"Total.Reads"]
oxy.full.X <- oxy.full.X[,!names(oxy.full.X) %in% c("Total.Reads")]
oxy.full.X <- cbind(animal.matrix[1:(nrow(animal.matrix)),], oxy.full.X)

#replace the column names
for(i in 2:ncol(oxy.full.Y)){
  names(oxy.full.Y)[i] <- as.character(oxy.taxa[1,names(oxy.full.Y)[i]])
}
for(i in 23:ncol(oxy.full.X)){
  names(oxy.full.X)[i] <- as.character(oxy.taxa[1,names(oxy.full.X)[i]])
}

save(oxy.taxa, file = "oxy.taxa.rda")
save(oxy.full.X, file = "X.MM.IF.rda")
save(oxy.full.Y, file ="Y.MM.IF.rda")
