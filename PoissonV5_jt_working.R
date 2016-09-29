#Notice a huge error in data processing in V2
#Fixing this in V3 (BJR 2016-05-27 10:08:40 -0700)
library(EnvStats)
library(glmnet)
library(ggplot2)

#read data
setwd("~/Dropbox/Shared/CMCI Projects/cmci_oxalate")
mb <- read.csv("ReadCounts.csv") #raw counts with total counts
#still need oxalate data

oxy.taxa <- read.csv("otu_table_NALB_relative_abundance_trans.csv")
oxy.taxa <- oxy.taxa[33:34,-1*(2:7)]
oxy.taxa <- droplevels(oxy.taxa)
#somehow the nuber of denovoXXXXX samples is different in this data frame than in "mb" - this is a problem.

missing <- names(mb)[!names(mb) %in% names(oxy.taxa)] #5 missing
mb[,missing]

metadata <- read.csv("metadata.csv")
names(metadata)[1] <- "Sample"
mb <- merge(metadata,mb)


extraStep <- c(unique(mb$Time.Point),1000)
extraSet2<-expand.grid(extraStep,unique(mb$Animal))
names(extraSet2)<-c("Time.Point","Animal")
oxy.full <-merge(extraSet2,mb,all.x=T)
oxy.full <- oxy.full[do.call("order",oxy.full[,c("Animal","Time.Point")]),]
theAnimal <- oxy.full$Animal
##### STOPPING HERE!!! FIX ME
oxy.full <- oxy.full[, ! names(oxy.full) %in% c("Time.Point","Animal","Sample","Oxalate.Concentration")]
oxy.full.save <- oxy.full
animal.matrix<-as.matrix(model.matrix(~theAnimal-1))


#oxy.full <- oxy.full[! is.na(oxy.full$Total.Reads), ] #Bonehead error

#let's just say keep the top 10% of the sample OTUs
sum(! names(oxy.full) %in% c("Oxalate.Consumed", "Oxalate.Excreted", "Oxalate.Degraded", "Total.Reads")) * 0.1
#keep top 624
otuReads <- colSums(oxy.full[,! names(oxy.full) %in% c("Oxalate.Consumed", "Oxalate.Excreted", "Oxalate.Degraded", "Total.Reads")], na.rm = T)
keepOTU <- names(otuReads)[otuReads > quantile(otuReads,0.9)]

top10 <- oxy.full[,union(keepOTU,c("Oxalate.Consumed", "Oxalate.Excreted", "Oxalate.Degraded", "Total.Reads"))]
summary(colSums(top10[,! names(top10) %in% c("Oxalate.Consumed", "Oxalate.Excreted", "Oxalate.Degraded", "Total.Reads")],na.rm=T))
hist(log(colSums(top10[,! names(top10) %in% c("Oxalate.Consumed", "Oxalate.Excreted", "Oxalate.Degraded", "Total.Reads")],na.rm=T)), xlab = "Log(OTU Reads)")

#I'm changing the way Oxalate.* is part of the model -- BJR 2016-05-27 09:32:32 -0700
yNames <- setdiff(names(top10), c("Oxalate.Consumed"))
xNames <- setdiff(names(top10), c("Oxalate.Consumed","Oxalate.Excreted","Oxalate.Degraded", "Total.Reads"))
oxy.full.Y <- top10[2:nrow(top10), ]
oxy.full.X <- top10[1:(nrow(top10)-1),]
oxy.full.X[,xNames] <- oxy.full.X[,xNames] / oxy.full.X[,"Total.Reads"]
oxy.full.X <- oxy.full.X[,!names(oxy.full.X) %in% c("Oxalate.Excreted","Oxalate.Degraded","Total.Reads")] 
oxy.full.X <- cbind(animal.matrix[1:(nrow(animal.matrix)-1),], oxy.full.X)

useMe<- ! (is.na(oxy.full.Y$Oxalate.Consumed)|is.na(oxy.full.X$Oxalate.Consumed))
oxy.full.X <- oxy.full.X[useMe,]
oxy.full.Y <- oxy.full.Y[useMe,]

#I think that oxalate consumed needs to be shifted like this to refelect
#the amount of oxalate consumed in the previous week.
oxy.full.X$Oxalate.Consumed <- oxy.full.Y$Oxalate.Consumed  


#workhorse for performing glmnet
#old ... args: grouped = FALSE, offset=theOffset, family = theFamily, pmax = nrow(oxy.full.X) - 1
#seq.alpha defines the alpha values to loop over
#nCV is the number of times to run cv.glmnet for a particular alpha values (to determine lambda)
#pSuccess is the percent of the time an alpha value needed to successfully return a model
#nlambda is as defined in glmnet
#family is as defined in glmnet
#findBetas <- function (index, seq.alpha = seq(0.5,1,0.1), nCV = 500, pSuccess = 0.2, nlambda = 100, grouped = TRUE, ...){
findBetas <- function(index){
  seq.alpha = seq(0.5,1,0.1)
  nCV = 500
  pSuccess = 0.2
  nlambda = 100
  nAlpha <- length(seq.alpha)
  maxIter <- nAlpha*nCV #perform nCV runs per alpha level
  #loop across alpha values
  MSE <- foreach(i = 1:maxIter,.packages='glmnet') %dopar% tryCatch(
    #cv.glmnet(as.matrix(oxy.full.X),oxy.full.Y[,index], alpha=seq.alpha[(i-1) %/% nCV + 1], nlambda=nlambda, grouped = grouped, ...)$cvm,
    cv.glmnet(as.matrix(oxy.full.X),oxy.full.Y[,index], alpha=seq.alpha[(i-1) %/% nCV + 1], nlambda=nlambda, grouped = grouped)$cvm,
    error = function(e){
      print(paste("Index",i,"produced no output.")); 
      rep(NA,nlambda)}
    )
  
  #convert the MSE list into a matrix
  MSE.mat <- sapply(MSE,"[",1:nlambda)
  MSE.long <- MSE.mat[,1:nCV]
  for(i in 1:(nAlpha-1)) MSE.long <- rbind(MSE.long,MSE.mat[,i*nCV+1:nCV])
  
  failures <- matrix(rowSums(is.na(MSE.long)) > (1-pSuccess)*nCV,nlambda,nAlpha) #had to have at least pSuccess to keep data from that alpha value
  mMSE <- matrix(rowMeans(MSE.long),nlambda,nAlpha) #get mean MSE values across lambda values (rows) then reshape to 100 (number of lambdas by default) x nAlpha
  sdMSE <- matrix(apply(MSE.long,1,var)^0.5,nlambda,nAlpha)
  mMSE[failures] <- NA
  sdMSE[failures] <- NA
  mMSE <<- mMSE #diagnostic
  indices.min <- apply(mMSE+sdMSE,2,which.min)
  
  findAIC <- function(glmnet.mod, nY,lambda){
    nP <- nrow(summary(coef(glmnet.mod,lambda)))
    theOffset <- eval(as.list(glmnet.mod$call)$offset)
    #this isn't thorough, only going to check for poisson
    #otherwise use normal
    prob <- ifelse(sum(class(glmnet.mod) == "fishnet"),
                   sum(dpois(oxy.full.Y[,nY],predict(glmnet.mod,as.matrix(oxy.full.X),s = lambda,type = "response",offset = theOffset),log=T)),
                   sum(dnorm(oxy.full.Y[,nY],predict(glmnet.mod,as.matrix(oxy.full.X),s = lambda,type = "response",offset = theOffset),log=T))
    )
    2*nP - 2*prob
  }
  
  AIC <- c()
  for(i in 1:nAlpha){
    if(1-length(indices.min[[i]])){AIC[i] <- Inf; next} #need this for cases where no index for the minimum exists
    #cur.mod <- glmnet(as.matrix(oxy.full.X),oxy.full.Y[,index], alpha=seq.alpha[i], nlambda = nlambda, ...)
    cur.mod <- glmnet(as.matrix(oxy.full.X),oxy.full.Y[,index], alpha=seq.alpha[i], nlambda = nlambda)
    AIC[i] <- findAIC(cur.mod,index,cur.mod$lambda[indices.min[[i]]])
  }
  AICg <<- AIC #write this out for diagnostic purposes
  
  if(min(AIC) == Inf) return(data.frame()) #return nothing if every alpha level failed
  
  mod.index <- which.min(AIC)
  final.model <- glmnet(as.matrix(oxy.full.X),oxy.full.Y[,index], alpha=seq.alpha[mod.index], nlambda = nlambda, ...)
  #final.model <- glmnet(as.matrix(oxy.full.X),oxy.full.Y[,index], alpha=seq.alpha[mod.index], nlambda = nlambda)
  final.lambda <- final.model$lambda[indices.min[[mod.index]]]
  
  out <- summary(coef(final.model,final.lambda))
  if(nrow(out) > 0){
    out$j <- index
    out$alpha <- seq.alpha[mod.index]
    out$lambda <- final.lambda
    out$dev.ratio <- final.model$dev.ratio[which(final.lambda == final.model$lambda)]
    out$AIC <- AIC[mod.index]
  }
  as.data.frame(out)
  
}

library(doParallel)
registerDoParallel(cores=7) #24 for big mac!

startT <- Sys.time()

stopIndex <- (ncol(oxy.full.Y)-1) # last column is total reads, don't use it
#stopIndex <- 4
#pf <- c(rep(1,ncol(oxy.full.X)-1),0) #penalty factor vector that forces inclusion of oxalate.consumed
pf <- rep(1,ncol(oxy.full.X)) #do not force oxalate.consumed
out <- c()
usePoisson <- !names(oxy.full.Y) %in% c("Oxalate.Consumed","Oxalate.Excreted","Oxalate.Degraded","Total.Reads") 
for(i in 1:stopIndex){
  print(i)
  try(ifelse(usePoisson[i],
             out <- rbind(out, findBetas(i, offset = log(oxy.full.Y$Total.Reads), family = "poisson", grouped = FALSE, penalty = pf)),
             out <- rbind(out, findBetas(i, grouped = FALSE, penalty = pf))
                                ))
} 

endT <- Sys.time()



write.csv(out,"Top10-unforcedOxalate.csv",row.names = FALSE) 

reruns<-setdiff(1:627,out$j) #list of failures
View(oxy.full.Y[,reruns]) #yeah, these are pretty much garbage

#Poisson.csv has the results from using non-normalized X counts
#Poisson2.csv has the newer results using a normalized X counts, see line 38 above
#Poisson2_standardized.csv has the newest results using a normalized X counts and standardized glmnet
#Top10.csv has the analysis of the OTUs above the 90th quantile for read counts (624 otus)
####### !!!!!!!!!!!!! I realized I ran all of the above analyses without including
####### potential offsets for different animals.
##### JUST USE TOP10.CSV!!!!
out <- read.csv("Top10-unforcedOxalate.csv")

out[out$i == 632,]$j #373
oxy.taxa[1,names(oxy.full.Y)[setdiff(out[out$i == 632,]$j,625:627)]] #who was affected by oxalate consumed
hist(out$dev.ratio,100,xlab="Pseudo-R^2", main="")
hist(table(out$i[out$i != 1]), 54, main = "Out Degree Distribution", xlab ="Out Degree", right = F)
hist(table(out$j[out$i != 1]), 31, main = "In Degree Distribution", xlab ="In Degree", right = F)

out.clean <- out[out$dev.ratio > 0.02,] #eliminates 38 models
hist(unique(out.clean$dev.ratio),30) #pretty even spread

meanReads <- colMeans(oxy.full.X)
meanReads[1:22] <- 1 #animal effects need to have a 1
meanReads <- c(1,meanReads) #add a 1 for the intercept

out.clean <- out.clean[out.clean$i != 1, ]
out.clean[out.clean$i == 632,]$j #forced
out.clean$x.scaled <- out.clean$x * meanReads[out.clean$i]

hist(table(out.clean$i[!out.clean$i %in% c(1,632)]), 20, main = "Out Degree Distribution", xlab ="Out Degree", right = F)
hist(table(out.clean$j), 20, main = "In Degree Distribution", xlab ="In Degree", right = F)

length(table(out.clean$j[!out.clean$i %in% c(1,632)])) #489 taxa
table(unique(out.clean[!out.clean$i %in% c(1,632),c("alpha","lambda")])$alpha)

alphaBetter <- unique(out.clean[,c("alpha","dev.ratio")])
plot(alphaBetter)
lines(supsmu(alphaBetter$alpha,alphaBetter$dev.ratio),col="red")

out.clean[out.clean$i %in% 2:7, ] #main effects of individual
table(out.clean[out.clean$i %in% 2:7, "i"])

oxyCoefs <- out.clean[out.clean$i == 632,c("j","x.scaled")]

hist(oxyCoefs$x.scaled, 50)
theTaxa <- names(oxy.full.Y)[oxyCoefs$j[c(-579,-580)]]
theOrder <- order(oxyCoefs$x.scaled[c(-579,-580)])
oxy.taxa[1,theTaxa[theOrder]] #throws an error
setdiff(theTaxa,names(oxy.taxa)) #somehow denovo804699 is not in oxy.taxa
orderedTaxa <- unname(apply(oxy.taxa[1,setdiff(theTaxa[theOrder],setdiff(theTaxa,names(oxy.taxa)))],2,as.character))
correctedX <- (oxyCoefs$x.scaled[theOrder])[oxyCoefs$j[c(-579,-580)] != which(theTaxa == setdiff(theTaxa,names(oxy.taxa)))]
oxyEffects <- data.frame(Taxa = orderedTaxa, Oxalate.Effect = correctedX)
oxyEffects$Taxa <- as.character(oxyEffects$Taxa)
unique(oxyEffects$Taxa) #hmm there are only 43 unique taxa names

hist(oxyEffects$Oxalate.Effect[oxyEffects$Taxa == "Ruminococcus flavefaciens"],20)

hist(aggregate(oxyEffects$Oxalate.Effect,by=list(oxyEffects$Taxa),mean)$x,20,main="Aggregated Oxalate Effects")

originalTaxa <- unique(apply(oxy.taxa[1,],2,as.character)) #96 of these
setdiff(originalTaxa,oxyEffects$Taxa) #things that got no love

colX <- c(0,colMeans(oxy.full.X))
out$effect <- out$x * colX[out$i]
hist(out$effect[out$effect != 0],1000, xlim = c(-0.5,0.5), xlab = "Effect Size", main = "")
out$meansize <- colX[out$i]

#test <- glmnet(as.matrix(oxy.full.X),oxy.full.Y[,1],offset= log(oxy.full.Y[,"Total.Reads"]), intercept=TRUE, alpha = 0.7833597, standardize=TRUE, family = "poisson")
#test.cv <- cv.glmnet(as.matrix(oxy.full.X),oxy.full.Y[,1],offset= log(oxy.full.Y[,"Total.Reads"]), intercept=TRUE, alpha = 0.7833597, standardize=TRUE, family = "poisson",nfolds = nrow(oxy.full.X), grouped = FALSE,pmax=nrow(oxy.full.X))


#preds <- predict(test.cv, as.matrix(oxy.full.X), offset = log(oxy.full.Y$Total.Reads), s = "lambda.min", type = "response") 
#plot(preds,type = "l", log = "y", ylim = c(10,12000), xlab = "Sample", ylab = "Read Count", main = "Predicted vs Observed Read Counts for OTU-1")
#points(10*oxy.full.Y[,out[1,"j"]], col = "red", type = "b")
####NOTE: IT WOULD BE BETTER TO PLOT THESE INDIVIDUALLY FOR ANIMALS



library(reshape2)
wide <- dcast(out, j ~ i, value.var = "x" )
names(wide) <- c("Y.OTU",c("Intercept",names(oxy.full.X))[as.numeric(names(wide)[-1])])
wide$Y.OTU <- names(oxy.full.Y)[wide$Y.OTU]
wide[is.na(wide)]<-0



library(cluster)
clusters <- clusGap(wide[,-1], FUN = kmeans, K.max = 10) #6 clusters
#fails: clusters <- clusGap(t(wide[,2:7]), FUN = fanny, K.max = 2)
plot(clusters)


#### make heatmap. works with "PoissonV5.R code"
#### most of the code was snagged from "process.data.R code"
cut.off <- 0.0
library("igraph")
library(stringr)
out.hist <- out.clean
out.hist$effect <- out.hist$x * colX[out.hist$i]  
out.hist$i.denovo <- as.factor(c("intercept", names(oxy.full.X))[out.hist$i])
out.hist$j.denovo <- as.factor(names(oxy.full.Y)[out.hist$j])
out.hist$j.taxa <- as.character(out.hist$j.denovo)
out.hist$i.taxa <- as.character(out.hist$i.denovo)
#fill in some names. will replace later with taxa look up
#### taxa lookup snagged from "process.data.R code"
for(i in 1:dim(out.hist)[1]){
  if(as.character(out.hist[i,]$j.denovo) %in% names(oxy.taxa)==TRUE){
    out.hist[i,]$j.taxa <- as.character(oxy.taxa[1,as.character(out.hist[i,]$j.denovo)])
  }}
for(i in 1:dim(out.hist)[1]){
  if(as.character(out.hist[i,]$i.denovo) %in% names(oxy.taxa)==TRUE){
    out.hist[i,]$i.taxa <- as.character(oxy.taxa[1,as.character(out.hist[i,]$i.denovo)])
  }}
out.hist$i.taxa <- as.character(str_trim(out.hist$i.taxa))
out.hist$j.taxa <- as.character(str_trim(out.hist$j.taxa))
tmp <- na.omit(out.hist[,c("i.taxa","j.taxa", "effect")])
names(tmp) <- c("from", "to", "Weight")
net <- graph.data.frame(d=tmp, directed=TRUE)
E(net)$weight <- abs(E(net)$Weight)
E(net)$color <- ifelse(E(net)$Weight > 0, "blue", "red")
E(net)$width <- E(net)$weight
net.ag <- simplify(net, remove.loops = F, remove.multiple=TRUE, edge.attr.comb = getIgraphOpt("edge.attr.comb"))
net.ag2 <- simplify(net, remove.loops=F, remove.multiple=TRUE, edge.attr.comb= list(Weight="mean"))
E(net.ag)$Weight <- E(net.ag2)$Weight
E(net.ag)$width<- log(abs(E(net.ag)$Weight)) + 3
E(net.ag)$color <- ifelse(E(net.ag)$Weight > 0, "blue", "red")
E(net.ag)$weight <- abs(E(net.ag)$Weight)+0.1
net.ag <- delete.edges(net.ag, E(net.ag)[weight<cut.off])
net.ag <- delete.vertices(net.ag,which(degree(net.ag)<1))
mat <- get.adjacency(net.ag, attr="Weight", sparse=F)
melt.mat <- melt(mat)
names(melt.mat) <- c("from", "to", "weight")
melt.mat$weight[melt.mat$weight==0]<-NA
melt.mat$trans.weight <- sign(melt.mat$weight)*abs(melt.mat$weight)^(1/3)
fromorder <- names(sort(table(melt.mat[abs(melt.mat$trans.weight) >0,]$from), decreasing=T))
#melt.mat <- subset(melt.mat, from!="Unassigned")
#melt.mat <- subset(melt.mat, to!="Unassigned")
melt.mat$from <- factor(melt.mat$from, level = fromorder)
melt.mat$to <- factor(melt.mat$to, level = fromorder)
names(melt.mat)[1:2] <- c("From", "To")
#lighten grey NAs, get grid lines, fix transformation
p <- ggplot(melt.mat, aes(x=From, y=To)) + 
  geom_tile(aes(fill = trans.weight),colour="white") + 
  theme(axis.text.x=element_text(angle = 90, hjust = 1)) + 
  scale_x_discrete(limits=melt.mat$From) + 
  scale_fill_gradient2(na.value="grey", name="Interaction Strength", low="blue", mid='white', high='red', breaks=c(min(na.omit(melt.mat$trans.weight)),max(na.omit(melt.mat$trans.weight))),labels=c("highly negative", "highly positive"))
p


#add phylotree
library(phytools)
library(phyloseq)
library(data.table)
#source("https://bioconductor.org/biocLite.R")
#biocLite("ggtree")
library(ggtree)
library(ggdendro)
library(grid)
library(gridExtra)
grntree <- read_tree_greengenes("gg_13_5_otus_99_annotated.tree")
grntax <- read.delim(file = "gg_13_5_taxonomy.txt", header=FALSE, stringsAsFactors=FALSE)
grntax[,2] <- gsub(";|\\[|\\]| ","",grntax[,2])
grntax[,2] <- gsub("Enterobacteriaceae","Enterobacteraceae",grntax[,2])
df <- matrix(,nrow=nrow(grntax),ncol=8)
for(i in 1:nrow(grntax)){
  df[i,1:length(df2)] <- df2 <- t(as.matrix(unlist(strsplit(grntax[i,2],".__"))))
  if( !is.na(df[i,8])){
    df[i,8] <- paste(df[i,7],df[i,8],sep=" ")
  }
}
df<-df[,-1]
grntax <- cbind(grntax[,1],grntax[,2],df)
colnames(grntax) <- c("idnum","string","k","p","c","o","f","g","s")
grntax[grntax==""] <- NA
grntax[,1] <- as.numeric(as.character(grntax[,1]))
tax <- as.matrix(unique(melt.mat$To))
tax[,1] <- str_trim(tax[,1])
colnames(tax) <- "name"
species <- merge(grntax,tax,by.x="s",by.y="name")
species <- subset(species, select = c(idnum,k,p,c,o,f,g,s))
species$idnum <- as.numeric(as.character(species$idnum))
species <- species[sample(nrow(species)),]
short_s <- aggregate(idnum ~ s, FUN = max, data=species)
genre <- merge(grntax,tax,by.x="g",by.y="name")
genre <- genre[which(is.na(genre$s)),]
genre <- subset(genre, select = c(idnum,k,p,c,o,f,g,s))
genre$idnum <- as.numeric(as.character(genre$idnum))
genre <- genre[sample(nrow(genre)),]
short_g <- aggregate(idnum ~ g, FUN = max, data=genre)
family <- merge(grntax,tax,by.x="f",by.y="name")            
family <- family[which(is.na(family$g)),]
family <- subset(family, select = c(idnum,k,p,c,o,f,g,s))
family$idnum <- as.numeric(as.character(family$idnum))
family <- family[sample(nrow(family)),]
short_f <- aggregate(idnum ~ f, FUN = max, data=family)
order <- merge(grntax,tax,by.x="o",by.y="name")
order <- order[which(is.na(order$f)),]
order <- subset(order, select = c(idnum,k,p,c,o,f,g,s))
order$idnum <- as.numeric(as.character(order$idnum))
order <- order[sample(nrow(order)),]
short_o <- aggregate(idnum ~ o, FUN = max, data=order)
class <- merge(grntax,tax,by.x="c",by.y="name")
class <- order[which(is.na(order$o)),]
class <- subset(class, select = c(idnum,k,p,c,o,f,g,s))  #######code is incomplete if you want to go past class (none in our case)
class <- class[sample(nrow(class)),]
tips_to_plot <- rbind(species,genre,family,order,class)
short_tips_to_plot <- rbindlist(list(short_s,short_g,short_f,short_o))
tips_to_plot$idnum <- as.character(tips_to_plot$idnum)

#write.table(tips_to_plot,file="tips_to_plot_all_forced.txt")
#tips_to_plot <- read.csv("tips_to_plot_all.txt", sep="")
#write.table(short_tips_to_plot,file="shor_tips_to_plot_all_forced.txt")
#short_tips_to_plot <- read.csv("short_tips_to_plot_all.txt", sep="")

#need to go list off taxa that need to be added to tree by hand
print("Some taxa that you supplied were not found in greengenes;")
tax.found <- list()
tax.unfound <- list()
for(i in 1:nrow(tax)){
  if(tax[i,1] %in% tips_to_plot$s){
    tax.found <- c(tax.found,tax[i,1])
  }else if(tax[i,1] %in% tips_to_plot$g){
    tax.found <- c(tax.found,tax[i,1])
  }else if(tax[i,1] %in% tips_to_plot$f){
    tax.found <- c(tax.found,tax[i,1])
  }else if(tax[i,1] %in% tips_to_plot$o){
    tax.found <- c(tax.found,tax[i,1])
  }else if(tax[i,1] %in% tips_to_plot$c){
    tax.found <- c(tax.found,tax[i,1])
  }else{
    tax.unfound <- c(tax.unfound,tax[i,1])
    print(tax[i,1])
  }
}
tax.found <- as.matrix(tax.found)
tax.unfound <- as.matrix(tax.unfound)

#plot all 
pruned.tree <- drop.tip(grntree,setdiff(grntree$tip.label,tips_to_plot$idnum))
#plot(pruned.tree)

#plot short list of tips
pruned.tree.short <- drop.tip(grntree,setdiff(grntree$tip.label,short_tips_to_plot$idnum))
plot(pruned.tree.short)

#well....there are way more taxa in the greengenes taxonomy file than in
#the greengenes tree file. Need to pick ids that exist in greengenes tree file
df <- tips_to_plot
new.tree <- pruned.tree
new.tree2 <- pruned.tree
rownames(df) <- df$idnum
df <- df[,-1]

tax.found2 <- tax.found
row.names(tax.found2) <- tax.found2[,1]

silva <- read.tree("LTPs123_SSU_tree.newick")
nodes <- matrix
for(i in 1:length(tax.found)){
  nodes <- c(nodes,grep(as.character(tax.found[i,1]),silva$node.label))
}

#loops through phylogeny and replaces tip numbers with names
for(i in 1:length(new.tree$tip.label)){
  count <- as.matrix(colSums(!is.na(df[pruned.tree$tip.label[i],])))
  name <- as.character(df[pruned.tree$tip.label[i],sum(count[,1])])
  new.tree$tip.label[i] <- as.character(tax.found2[name,1])
  tax.found2[name,1] <- paste(name,"1",sep=" ")
}
#gotta fix this damn clostridiales thing. my work around is not great
new.tree$tip.label[new.tree$tip.label=="Clostridiales"] <- "Clostridiales 1"
new.tree$tip.label[[7469]] <- "Clostridiales"

new.tree.pick <- drop.tip(new.tree,as.character(tax.found2[,1]))
plot.phylo(new.tree.pick,align.tip.label = TRUE)
plot(new.tree.pick,use.edge.length=FALSE,cex=0.75)

melt.mat2 <- melt.mat
fromorder2 <- rbind(tax.unfound,as.matrix(new.tree.pick$tip.label))
toorder2 <- rbind(tax.unfound,as.matrix(new.tree.pick$tip.label))
#check this out if you are running different data though here;
rm.from <- c("Oxalate.Excreted","Oxalate.Degraded","denovo804699")
rm.to <- c("Oxalate.Consumed","theAnimalNA66","theAnimalNA67","theAnimalNA58","theAnimalNA70","theAnimalNA71","theAnimalNA68","denovo804699")
for(i in 1:length(rm.from)){
  melt.mat2 <- subset(melt.mat2, From!=rm.from[i])
  fromorder2 <- subset(fromorder2, fromorder2[,1]!=rm.from[i])
}
for(i in 1:length(rm.to)){
  melt.mat2 <- subset(melt.mat2, To!=rm.to[i])
  toorder2 <- subset(toorder2, toorder2[,1]!=rm.to[i])
}
fromorder <- read.csv("C:/Users/JT/Dropbox/cmci_oxalate/fromorder.txt", header=FALSE)
fromorder <- as.matrix(apply(fromorder, 2, rev))
toorder <- read.csv("C:/Users/JT/Dropbox/cmci_oxalate/toorder.txt", header=FALSE)
toorder <- as.matrix(apply(toorder, 2, rev))
melt.mat2$From <- factor(melt.mat2$From, level = fromorder)
melt.mat2$To <- factor(melt.mat2$To, level = toorder)
melt.mat2 <- subset(melt.mat2, !is.na(To))
melt.mat2 <- subset(melt.mat2, !is.na(From))
#lighten grey NAs, get grid lines, fix transformation
p2 <- ggplot(melt.mat2, aes(x=From, y=To)) + 
  geom_tile(aes(fill = trans.weight),colour="white") + 
  theme(axis.text.y=element_blank(), axis.text.x=element_text(angle = 90, hjust = 1), axis.title.y=element_blank(), panel.grid.major=element_blank(), axis.ticks=element_blank()) + 
  scale_fill_gradient2(na.value="grey", low="blue", mid='white', high='red', breaks=c(min(na.omit(melt.mat$trans.weight)),max(na.omit(melt.mat$trans.weight,labels=element_blank()))))
p2


#alright, try again.
library(geiger)
library(phylobase)
library(ape)

#functions
find_best_mrca <- function(phy, group, threshold){
  
  group_matches <- phy$tip.label[grepl(group, phy$tip.label, ignore.case=TRUE)]
  group_mrca <- getMRCA(phy,phy$tip.label[grepl(group, phy$tip.label, ignore.case=TRUE)])
  group_leaves <- tips(phy, group_mrca)
  match_ratio <- length(group_matches)/length(group_leaves)
  
  if( match_ratio < threshold){
    
    #start searching for children nodes that have more than 95% of descendants matching to the search pattern
    mrca_children <- descendants(as(phy,"phylo4"), group_mrca, type="all")
    i <- 1
    new_ratios <- NULL
    nleaves <- NULL
    names(mrca_children) <- NULL
    
    for(new_mrca in mrca_children){
      child_leaves <- tips(tree.test, new_mrca)
      child_matches <- grep(group, child_leaves, ignore.case=TRUE)
      new_ratios[i] <- length(child_matches)/length(child_leaves)
      nleaves[i] <- length(tips(phy, new_mrca))
      i <- i+1
    }
    
    
    
    match_result <- data.frame(mrca_children, new_ratios, nleaves)
    
    
    match_result_sorted <- match_result[order(-match_result$nleaves,match_result$new_ratios),]
    found <- numeric(0);
    
    print(match_result_sorted)
    
    for(line in 1:nrow(match_result_sorted)){
      if(match_result_sorted$ new_ratios[line]>=threshold){
        return(match_result_sorted$mrca_children[line])
        found <- 1
      }
      
    }
    
    if(found==0){return(found)}
  }else{return(group_mrca)}
  
  
  
  
}

add_triangle <- function(phy, group,phylo_plot){
  
  group_node_labels <- phy$tip.label[grepl(group, phy$tip.label)]
  group_mrca <- getMRCA(phy,group_node_labels)
  group_nodes <- descendants(as(tree.test,"phylo4"), group_mrca, type="tips")
  names(group_nodes) <- NULL
  
  x<-phylo_plot$xx
  y<-phylo_plot$yy
  
  
  x1 <- max(x[group_nodes])
  x2 <-max(x[group_nodes])
  x3 <- x[group_mrca]
  
  y1 <- min(y[group_nodes])
  y2 <- max(y[group_nodes])
  y3 <-  y[group_mrca]
  
  xcoords <- c(x1,x2,x3)
  ycoords <- c(y1,y2,y3)
  
  polygon(xcoords, ycoords)
  
  return(c(x2,y3))
  
}



#main

#cat("((A_1:0.05,E_2:0.03,A_3:0.2,A_4:0.1,A_5:0.1,A_6:0.1,A_7:0.35,A_8:0.4,A_9:01,A_10:0.2):0.9,((((B_1:0.05,B_2:0.05):0.5,B_3:0.02,B_4:0.04):0.6,(C_1:0.6,C_2:0.08):0.7):0.5,(D_1:0.3,D_2:0.4,D_3:0.5,D_4:0.7,D_5:0.4):1.2):0.5);", file = "ex.tre", sep = "\n")
#tree.test <- read.tree("ex.tre")

tree.test <- new.tree
# Step 1: Find the best MRCA that matches to the keywords or search patten

groups <- as.character(tax.found)
group_labels <- groups

group_edges <- numeric(0)
edge.width <- rep(1, nrow(new.tree$edge))
count <- 1


for(group in groups){
  
  best_mrca <- find_best_mrca(tree.test, group, 0.90)
  
  group_leaves <- tips(tree.test, best_mrca)
  
  groups[count] <- paste(group_leaves, collapse="|")
  group_edges <- c(group_edges,best_mrca)
  
  #Step2: Remove the edges of the branches that will be collapsed, so they become invisible
  edge.width[tree.test$edge[,1] %in% c(group_edges[count],descendants(as(tree.test,"phylo4"), group_edges[count], type="all")) ] <- 0
  count = count +1
  
}


#Step 3: plot the tree hiding the branches that will be collapsed/grouped

last_plot.phylo <- plot(tree.test, show.tip.label = F, edge.width = edge.width)

#And save a copy of the plot so we can extract the xy coordinates of the nodes
#To get the x & y coordinates of a plotted tree created using plot.phylo
#or plotTree, we can steal from inside tiplabels:
last_phylo_plot<-get("last_plot.phylo",envir=.PlotPhyloEnv)

#Step 4: Add triangles and labels to the collapsed nodes
for(i in 1:length(groups)){
  
  text_coords <- add_triangle(tree.test, groups[i],last_phylo_plot)
  
  text(text_coords[1],text_coords[2],labels=group_labels[i], pos=4)
  
}
