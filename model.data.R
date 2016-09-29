library(EnvStats)
library(glmnet)
library(ggplot2)
library(doParallel); packageVersion("doParallel")

#read data
setwd("C:/Users/JT/Dropbox/cmci_oxalate/microbiome_code")
load("X.MM.IF.rda")
load("Y.MM.IF.rda")
registerDoParallel(cores=20) #24 for big mac!
outfile <- "Top10-unforcedOxalate.csv"


#print out X's and Y's. Good to check.
names(oxy.full.X)
names(oxy.full.Y)

findBetas <- function (index, seq.alpha = seq(0.5,1,0.1), nCV = 500, pSuccess = 0.2, nlambda = 100, grouped = TRUE, ...){
  nAlpha <- length(seq.alpha)
  maxIter <- nAlpha*nCV #perform nCV runs per alpha level
  
  #loop across alpha values
  MSE <- foreach(i = 1:maxIter) %dopar% tryCatch(
    cv.glmnet(as.matrix(oxy.full.X),oxy.full.Y[,index], alpha=seq.alpha[(i-1) %/% nCV + 1], nlambda=nlambda, grouped = grouped, ...)$cvm, 
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
    cur.mod <- glmnet(as.matrix(oxy.full.X),oxy.full.Y[,index], alpha=seq.alpha[i], nlambda = nlambda, ...)
    AIC[i] <- findAIC(cur.mod,index,cur.mod$lambda[indices.min[[i]]])
  }
  AICg <<- AIC #write this out for diagnostic purposes
  
  if(min(AIC) == Inf) return(data.frame()) #return nothing if every alpha level failed
  
  mod.index <- which.min(AIC)
  final.model <- glmnet(as.matrix(oxy.full.X),oxy.full.Y[,index], alpha=seq.alpha[mod.index], nlambda = nlambda, ...)
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


registerDoParallel(cores=34) #24 for big mac!
startT <- Sys.time()

stopIndex <- (ncol(oxy.full.Y)-1) # last column is total reads, don't use it
#stopIndex <- 4
#pf <- c(rep(1,ncol(oxy.full.X)-1),0) #penalty factor vector that forces inclusion of oxalate.consumed
pf <- rep(1,ncol(oxy.full.X)) #do not force oxalate.consumed
out <- c()
usePoisson <- !names(oxy.full.Y) %in% c("Total.Reads") 
for(i in 1:stopIndex){
  print(i)
  try(ifelse(usePoisson[i],
             out <- rbind(out, findBetas(i, offset = log(oxy.full.Y$Total.Reads), family = "poisson", grouped = FALSE, penalty = pf)),
             out <- rbind(out, findBetas(i, grouped = FALSE, penalty = pf))
  ))
} 

endT <- Sys.time()

write.csv(out,outfile,row.names = FALSE) 

#reruns<-setdiff(1:627,out$j) #list of failures
#View(oxy.full.Y[,reruns]) #yeah, these are pretty much garbage



