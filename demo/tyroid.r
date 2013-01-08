library(sdnet)

data(tyr1)
data(tyr2)

n1 <- ncol(tyr1$cdata)
n2 <- ncol(tyr2$cdata)
N <- nrow(tyr1$cdata)
clslevs <- range(tyr1$cls)
n11 <- sum(tyr1$cls==clslevs[1])
n12 <- sum(tyr1$cls==clslevs[2])
ncats <- 3
genenames <- rownames(tyr1$cdata)
nodeCats <- lapply(1:N, function(i) 1:ncats)
names(nodeCats) <- genenames

## soft discretization
q <- sdnet::cnDiscretize(cbind(tyr1$cdata,tyr2$cdata), numcats=ncats, mode="soft", marginal="quantile", learnset=1:n1, cover=0.95)
ptest <- q$pdata[,(n1+1):(n1+n2)]
plearn <- q$pdata[,1:n1]
rm(q)

bnet <- sdnet::cnNew(genenames, cats=nodeCats, pars=vector("list",N), probs=NULL, dagonly=TRUE)
## sets P(X=k) \propto \sum_s q_k(y^s) = P(y^s|X=k)/(\sum_m P(y^s|X=m))
bnet <- sdnet::cnSetProb(bnet, plearn, nodeCats=nodeCats[bnet@nodes], softmode=TRUE)
net1 <- sdnet::cnSetProb(bnet, plearn[, tyr1$cls==clslevs[1]], nodeCats=nodeCats[bnet@nodes], softmode=TRUE)
net2 <- sdnet::cnSetProb(bnet, plearn[, tyr1$cls==clslevs[2]], nodeCats=nodeCats[bnet@nodes], softmode=TRUE)
rm(plearn)
  
pl1 <- -sapply(1:net1@numnodes, function(j) sum(bnet@probs[[j]]*log(net1@probs[[j]]/bnet@probs[[j]])))
pl1 <- 2*(n11+n12)*(n11/n12)*pl1
pvals <- 1 - pchisq(pl1, ncats-1)
ind <- order(pvals, decreasing=FALSE)
hc <- sapply(1:N, function(i) sqrt(N)*(i/N-pvals[ind[i]])/sqrt((i/N)*(1-i/N)))
k <- which(hc==max(hc[!is.nan(hc)]))[1]
if(is.na(k) || k<2) k <- 2
ind <- ind[1:k]

predict <- NULL
for(k in 1:ncol(ptest)) {
  dd <- t(matrix(ptest[,k], nrow=ncats))
  pl1 <- sapply(ind, function(j) sum(net1@probs[[j]]*log(dd[j,]/net1@probs[[j]])))
  pl2 <- sapply(ind, function(j) sum(net2@probs[[j]]*log(dd[j,]/net2@probs[[j]])))
  pdl <- pl2- pl1
  pdl <- sum(pdl[!is.nan(pdl)&abs(pdl)<Inf])
  sel <- clslevs[1]
  if(pdl >= 0)
    sel <- clslevs[2]
  predict <- c(predict, sel)
}
rm(ptest)

## compare the prediction to the true classes 
r <- cbind(tyr2$cls, predict)
colnames(r) <- c("tyr2$cls", "prediction")
r

acc <- 0.5*(sum(predict==clslevs[1]&tyr2$cls==clslevs[1])/sum(tyr2$cls==clslevs[1]) + sum(predict==clslevs[2]&tyr2$cls==clslevs[2])/sum(tyr2$cls==clslevs[2]))
cat("balanced accuracy of predicting (tyr1 -> tyr2) = ", acc, "\n")
