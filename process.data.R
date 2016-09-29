out <- read.csv("Top10-unforcedOxalate.csv")
out.clean <- out[out$dev.ratio > 0.02,] #eliminates 38 models

meanReads <- colMeans(oxy.full.X)
meanReads[1:6] <- 1 #animal effects need to have a 1
meanReads <- c(1,meanReads) #add a 1 for the intercept

out.clean <- out.clean[out.clean$i != 1, ]
out.clean$x.scaled <- out.clean$x * meanReads[out.clean$i]

out.clean$meansize.i <- meanReads[out.clean$i]
out.clean$meansize.j <- meanReads[out.clean$j]

out.clean$j.denovo <- as.factor(names(oxy.full.Y)[out.clean$j])
out.clean$i.denovo <- as.factor(c("intercept", names(oxy.full.X))[out.clean$i])
out.clean$j.taxa <- NA
out.clean$i.taxa <- NA


tmp <- as.data.frame(t(oxy.taxa[,2:dim(oxy.taxa)[2]]))
names(tmp) <- c("taxa","level")
tmp$denovo <- rownames(tmp)
rownames(tmp) <- NULL
tmp$denovo <- as.factor(tmp$denovo)

#### Taking too long to figure out how to properly pull taxa because of some denovo not in the list, just do a damn loop and figure out slick way later
for(i in 1:dim(out.clean)[1]){
  if(out.clean[i,]$j.denovo %in% names(oxy.taxa)==T){
  out.clean[i,]$j.taxa <- as.character(tmp[tmp$denovo==as.character(out.clean[i,]$j.denovo),]$taxa)
  }}
out.clean$j.taxa <- as.factor(out.clean$j.taxa)

for(i in 1:dim(out.clean)[1]){
  if(out.clean[i,]$i.denovo %in% names(oxy.taxa)==T){
    out.clean[i,]$i.taxa <- as.character(tmp[tmp$denovo==as.character(out.clean[i,]$i.denovo),]$taxa)
  }}
out.clean$i.taxa <- as.factor(out.clean$i.taxa)


library("igraph")
tmp <- out.clean[,c("i.denovo","j.denovo", "x.scaled", "meansize.i", "meansize.j")]
names(tmp) <- c("from", "to", "Weight", "meansize.otu.from", "meansize.otu.to")
net <- graph.data.frame(d=tmp, directed=TRUE)
E(net)$weight <- abs(E(net)$Weight)

cl <- walktrap.community(net, weights = E(net)$weight)
modularity(cl)
membership(cl)
V(net)$color <- membership(cl)
V(net)$cluster <- membership(cl)

meansize.taxa <- aggregate(meansize.otu.from~from,data=tmp,FUN=sum)
### a few vertices are only in the "to" category.  Add these to the list.
missing.taxa <- setdiff(union(unique(tmp$from), unique(tmp$to)),unique(tmp$from))
missing.tmp <- subset(tmp, to %in% missing.taxa)
meansize.missing <- aggregate(meansize.otu.to~to, data=missing.tmp,FUN=sum)
names(meansize.taxa) <- c("name", "meansize.taxa")
names(meansize.missing) <- c("name", "meansize.taxa")
taxa.size <- rbind(meansize.taxa, meansize.missing)
taxa.size$name <- as.character(taxa.size$name)
taxa.size <- taxa.size[match(V(net)$name,taxa.size$name),]
V(net)$meansize <- taxa.size$meansize.taxa

net.ag <- simplify(net, remove.loops = F, remove.multiple=TRUE, edge.attr.comb = getIgraphOpt("edge.attr.comb"))
net.ag2 <- simplify(net, remove.loops=F, remove.multiple=TRUE, edge.attr.comb= list(Weight="sum"))
E(net.ag)$Weight <- E(net.ag2)$Weight









###################################################
topClusters <- sort(table(cl$membership), decreasing=TRUE)
topClusters[1:10]
plot(topClusters, main="Cluster size", ylab="Number of members", type="b", lwd=2)
l <-  layout.fruchterman.reingold(net, weights = E(net)$weight)
plot(net, edge.width=0.1,edge.arrow.size=0, edge.width=E(net)$weight, vertex.size=2, layout=l, vertex.label=NA)

cut.off <- 0.09
net.small <- delete_edges(net, E(net)[weight<cut.off])


plot(net.small, edge.arrow.size=0.1, edge.width=0.2,layout=layout.fruchterman.reingold, vertex.label=NA, edge.curved=T, vertex.size=2)
plot(net.small, edge.arrow.size=0.1, layout=layout.kamada.kawai, vertex.label.cex=0.7, vertex.label.dist=0.15,  edge.curved=0.2, vertex.label.color="black")

plot(net.small, edge.arrow.size=0.1, layout=layout.kamada.kawai, vertex.shape="none", vertex.label.cex=V(net)$size/3,  edge.curved=0.2, vertex.label.color="black")



tmp <- na.omit(out[,c("i.taxa","j.taxa", "effect", "meansize", "meansize.j")])
names(tmp) <- c("from", "to", "Weight", "meansize.otu", "meansize.otu.to")
net <- graph.data.frame(d=tmp, directed=TRUE)
meansize.taxa <- aggregate(meansize.otu~from,data=tmp,FUN=sum)
### a few vertices are only in the "to" category.  Add these to the list.
missing.taxa <- setdiff(union(unique(tmp$from), unique(tmp$to)),unique(tmp$from))
missing.tmp <- subset(tmp, to %in% missing.taxa)
meansize.missing <- aggregate(meansize.otu.to~to, data=missing.tmp,FUN=sum)
names(meansize.taxa) <- c("name", "meansize.taxa")
names(meansize.missing) <- c("name", "meansize.taxa")
taxa.size <- rbind(meansize.taxa, meansize.missing)
taxa.size$name <- as.character(taxa.size$name)
taxa.size <- taxa.size[match(V(net)$name,taxa.size$name),]
V(net)$meansize <- taxa.size$meansize.taxa
#tmp <- merge(tmp, meansize.taxa, by="from")
#tmp <- tmp[tmp$weight < -0.037 | tmp$weight > 0.037,]
#tmp$weight <- abs(tmp$weight)
V(net)$size <- 0.5*log(V(net)$meansize)+4
E(net)$weight <- abs(E(net)$Weight)
E(net)$color <- ifelse(E(net)$Weight > 0, "blue", "red")
E(net)$width <- E(net)$weight

plot(net, edge.arrow.size=0.1, layout=layout.fruchterman.reingold, vertex.label.cex=0.7, vertex.label.dist=0.2,  edge.curved=T)
plot(net, edge.arrow.size=0.1, layout=layout.circle, vertex.label.cex=0.7, vertex.label.dist=0.2,  edge.curved=T)


cut.off <- 0.2
net.small <- delete_edges(net, E(net)[weight<cut.off])


plot(net.small, edge.arrow.size=0.1, layout=layout.fruchterman.reingold, vertex.label.cex=0.7, vertex.label.dist=0.15,  edge.curved=T, vertex.label.color="black")
plot(net.small, edge.arrow.size=0.1, layout=layout.kamada.kawai, vertex.label.cex=0.7, vertex.label.dist=0.15,  edge.curved=0.2, vertex.label.color="black")

plot(net.small, edge.arrow.size=0.1, layout=layout.kamada.kawai, vertex.shape="none", vertex.label.cex=V(net)$size/3,  edge.curved=0.2, vertex.label.color="black")

###########################
#############################
#### make a network where edges are all aggregated if going between same vertices.
net <- 
net.ag <- simplify(net, remove.loops = F, remove.multiple=TRUE, edge.attr.comb = getIgraphOpt("edge.attr.comb"))
net.ag2 <- simplify(net, remove.loops=F, remove.multiple=TRUE, edge.attr.comb= list(Weight="sum"))
E(net.ag)$Weight <- E(net.ag2)$Weight
V(net.ag)$size <- 0.5*log(V(net.ag)$meansize)+4
E(net.ag)$width<- log(abs(E(net.ag)$Weight)) + 3
E(net.ag)$color <- ifelse(E(net.ag)$Weight > 0, "blue", "red")
E(net.ag)$weight <- abs(E(net.ag)$Weight)+0.1
cut.off <- 0.2
net.ag <- delete.edges(net.ag, E(net.ag)[weight<cut.off])
net.ag <- delete.vertices(net.ag,which(degree(net.ag)<1))

plot(net.ag, edge.arrow.size=0.2, layout=layout.fruchterman.reingold, vertex.label.cex=0.7, vertex.label.dist=0.15,  edge.curved=0.2, vertex.label.color="black")
plot(net.ag, edge.arrow.size=0.1, layout=layout.kamada.kawai, vertex.label.cex=0.7, vertex.label.dist=0.2,  edge.curved=0.2, vertex.label.color="black")
plot(net.ag, edge.arrow.size=0.2, layout=layout.circle, vertex.label.cex=0.9, vertex.label.dist=0.25,  edge.curved=0.2, vertex.label.color="black")

###mat2 <- get.edges(net.ag)
mat <- get.adjacency(net.ag, attr="Weight", sparse=F)


#oxy_heatmap <- heatmap(mat, Rowv=F, Colv=F)

#heatmap.2(mat,trace="none")
library(ggplot2)
library(reshape2)
melt.mat <- melt(mat)
names(melt.mat) <- c("from", "to", "weight")
melt.mat$trans.weight <- sign(melt.mat$weight)*abs(melt.mat$weight)^(1/3)
#melt.mat <- subset(melt.mat, abs(trans.weight) > 0)
fromorder <- names(sort(table(melt.mat[abs(melt.mat$trans.weight) >0,]$from), decreasing=T))
#toorder <- names(sort(table(melt.mat$to), decreasing=T))
melt.mat <- subset(melt.mat, from!="Unassigned")
melt.mat <- subset(melt.mat, to!="Unassigned")

melt.mat$from <- factor(melt.mat$from, level = fromorder)
melt.mat$to <- factor(melt.mat$to, level = fromorder)
names(melt.mat)[1:2] <- c("From", "To")

p <- ggplot(melt.mat, aes(x=From, y=To)) + geom_raster(aes(fill = trans.weight), hjust=0.5, vjust=0.5, interpolate=FALSE) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_gradient2(name="Interaction\n Strength", low="blue", mid='white', high='red', breaks=c(-0.9,3.1),labels=c("highly negative", "highly positive") )
p



edgelist <- get.edgelist(net,names=TRUE) 
edges <- as.data.frame(edgelist, edlist)
weighted.edges <- aggregate(edges["Weight"], by=list(V1=edges$V1,V2=edges$V2), FUN=sum )
collapsed.graph <- graph.data.frame( weighted.edges ) 
E(collapsed.graph)$Weight

plot(collapsed.graph, edge.width=log(E(collapsed.graph)$Weight+0.1)/5, edge.arrow.size=0.1, vertex.size=2, layout=layout.kamada.kawai, vertex.label.cex=0.7, vertex.label.dist=0.2)










cl <- walktrap.community(net, weights = E(net)$weight)
modularity(cl)
membership(cl)
V(net)$color <- membership(cl)
V(net)$cluster <- membership(cl)

topClusters <- sort(table(cl$membership), decreasing=TRUE)
topClusters[1:10]
plot(topClusters, main="Cluster size", ylab="Number of members", type="b", lwd=2)
l <-  layout.kamada.kawai(net, weights = E(net)$weight)
plot(net.ag, edge.width=0.1,edge.arrow.size=0, edge.width=E(net)$weight, vertex.size=2, layout=l)

#### check out which cluster includes oxalate consumed
membership(cl)[which(names(membership(cl))=="Oxalate.consumed")]
### looks like it's in cluster number 8, find out who else is in there
ox.graph <- induced_subgraph(net, which(V(net)$cluster %in% 8))
plot(ox.graph, layout=layout.fruchterman.reingold, edge.width=0.5, vertex.size=5, edge.width=E(net)$weight, edge.arrow.size=0.2)
