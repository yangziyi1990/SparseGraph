########## import library #############
library(R.matlab)
library(space)
library(igraph)
library(glasso)
library(poweRlaw)
library(reshape2)
library(corpcor)
library(GeneNet)
library(CDLasso)
WhereAmI <- "D:/Ziyi/School/PMO/13.Network/1-code/4-ziyi/"

source(paste0(WhereAmI,"lib/glasso_SF.R"))
source(paste0(WhereAmI,"lib/BayesianGLasso.R"))
source(paste0(WhereAmI,"lib/PCA_CMI.R"))
source(paste0(WhereAmI,"lib/CMI2NI.R"))
source(paste0(WhereAmI,"lib/ENA.R"))

########## import simulated data ##########
dat <- readMat(paste0(WhereAmI,"Sample.mat"))
theta <- readMat(paste0(WhereAmI,"Theta.mat"))
emp_cov <- readMat(paste0(WhereAmI,"emp_cov.mat"))
true.adj <- readMat(paste0(WhereAmI,"adjacent.mat"))

sample <- t(dat$X)
theta <- theta$Theta
emp_cov <- emp_cov$S
true.adj <- true.adj$G

n <- dim(sample)[1]
p <- dim(sample)[2]

########## view the true simulated network ##########
plot.adj <- true.adj
diag(plot.adj)=0
temp=graph.adjacency(adjmatrix=plot.adj, mode="undirected")
temp.degree=apply(plot.adj, 2, sum)
V(temp)$color=(temp.degree>7)+3
plot(temp, vertex.size=5, vertex.frame.color="white",layout=layout.fruchterman.reingold, vertex.label=NA, edge.connectivity = NA)


#################### estimate the parcial correlation matrix with various methods
expr.mat <- sample
est_edge <- list()
method = "inverse.sum"
n.perm = 0
sig.quant = 0.8

pcor_est <- pcor.shrink(expr.mat)
test_pcor <- network.test.edges(pcor_est)

########### the values of lam1 were selected to make the results of different methods comparable.
##-----------------------------
### 1. neighborhood selection
##-----------------------------
cat("Constructing using neighborhood selection ...")
est_res.ns <- matrix(0, p, p)
lambda <- 0.65
alpha <- (1-lambda)^2 * p^2
for (k in 1:p) {
  rsp <- expr.mat[, k]
  prd <- t(expr.mat[, -k])
  lam <- sqrt(sum(rsp ^ 2)) * qnorm(alpha / (2 * p ^ 2), lower.tail = F)
  result1 <- l2.reg(prd, rsp, lambda = lam)
  est_res.ns[k, -k] <- result1$estimate
}
est_edge[[1]] <- (abs(est_res.ns) + t(abs(est_res.ns))) / 2

## edges detected
fit.adj=abs(est_edge[[1]])>1e-6
sum(fit.adj==1)/2                     ##total number of edges detected
sum(fit.adj[true.adj==1]==1)/2        ##total number of true edges detected

## plot network
edge.result1 <- est_edge[[1]]
edge.result1[which(edge.result1!=0)] <- 1
temp <- graph.adjacency(adjmatrix = edge.result1, mode = "undirected")
temp.degree <- apply(edge.result1 , 2, sum)
V(temp)$color <- (temp.degree>12)+3
plot(temp, vertex.size=5, vertex.frame.color="white",
     layout=layout.fruchterman.reingold, vertex.label=NA, edge.connectivity = NA)

cat(" OK\n")

##------------------
### 2. glasso
##------------------
cat("Constructing using GLASSO ...")
S <- t(expr.mat) %*% expr.mat / n
result2 <- glasso(S, rho = 0.25)
est_edge[[2]] <- abs(result2$wi)


## edges detected
fit.adj=abs(est_edge[[2]])>1e-6
sum(fit.adj==1)/2                     ##total number of edges detected
sum(fit.adj[true.adj==1]==1)/2        ##total number of true edges detected

## plot network
edge.result2 <- est_edge[[2]]
diag(edge.result2) <- 0
edge.result2[which(edge.result2!=0)] <- 1
temp <- graph.adjacency(adjmatrix = edge.result2, mode = "undirected")
temp.degree <- apply(edge.result2 , 2, sum)
V(temp)$color <- (temp.degree>10)+3
plot(temp, vertex.size=5, vertex.frame.color="white",
     layout=layout.fruchterman.reingold, vertex.label=NA, edge.connectivity = NA)

cat(" OK\n")

##-------------------
### glasso-sf
##-------------------
cat("Constructing using GLASSO-SF ...")
result3 <- glasso_sf(expr.mat, alpha = 0.17)
est_edge[[3]] <- abs(result3$wi)
# store prob
#est_edge.prob <- cbind(est_edge.prob, GLasso.sf = mapply(FUN = function(x,y) est_edge[[4]][x,y], est_edge.prob$Node1, est_edge.prob$Node2))
cat(" OK\n") 

## edges detected
fit.adj=abs(est_edge[[3]])>1e-6
sum(fit.adj==1)/2                     ##total number of edges detected
sum(fit.adj[true.adj==1]==1)/2        ##total number of true edges detected

## plot network
edge.result3 <- est_edge[[3]]
diag(edge.result3) <- 0
edge.result3[which(edge.result3!=0)] <- 1
temp <- graph.adjacency(adjmatrix = edge.result3, mode = "undirected")
temp.degree <- apply(edge.result3 , 2, sum)
V(temp)$color <- (temp.degree>10)+3
plot(temp, vertex.size=5, vertex.frame.color="white",
     layout=layout.fruchterman.reingold, vertex.label=NA, edge.connectivity = NA)

cat(" OK\n")
##---------------
### pcacmi
##---------------
cat("Constructing using PCACMI ...")
invisible(capture.output(result4 <- pca_cmi(t(expr.mat), 0.03)))
est_edge[[4]] <- abs(result4$Gval)

## edges detected
fit.adj=abs(est_edge[[4]])>1e-6
sum(fit.adj==1)/2                     ##total number of edges detected
sum(fit.adj[true.adj==1]==1)/2        ##total number of true edges detected

## plot network
edge.result4 <- est_edge[[4]]
diag(edge.result4) <- 0
edge.result3[which(edge.result4!=0)] <- 1
temp <- graph.adjacency(adjmatrix = edge.result4, mode = "undirected")
temp.degree <- apply(edge.result4 , 2, sum)
V(temp)$color <- (temp.degree>4)+3
plot(temp, vertex.size=5, vertex.frame.color="white",
     layout=layout.fruchterman.reingold, vertex.label=NA, edge.connectivity = NA)

cat(" OK\n")

##------------------
### space
##------------------
cat("Constructing using SPACE ...")
invisible(capture.output(result5 <- space.joint(expr.mat, lam1 = 0.3 * n, iter = 5)))
est_edge[[5]] <- abs(result5$ParCor)

## edges detected
fit.adj=abs(est_edge[[5]])>1e-6
sum(fit.adj==1)/2                     ##total number of edges detected
sum(fit.adj[true.adj==1]==1)/2        ##total number of true edges detected

## plot network
edge.result5 <- est_edge[[5]]
diag(edge.result5) <- 0
edge.result5[which(edge.result5!=0)] <- 1
temp <- graph.adjacency(adjmatrix = edge.result5, mode = "undirected")
temp.degree <- apply(edge.result5 , 2, sum)
V(temp)$color <- (temp.degree>10)+3
plot(temp, vertex.size=5, vertex.frame.color="white",
     layout=layout.fruchterman.reingold, vertex.label=NA, edge.connectivity = NA)

cat(" OK\n")





