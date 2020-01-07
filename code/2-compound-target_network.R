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
dat <- read.table(paste0(WhereAmI,"3-combine-matrix.txt"), header = T)
name.target <- row.names(dat)
name.comp <- colnames(dat)

n <- dim(dat)[1]
p <- dim(dat)[2]


expr.mat <- dat
est_edge <- list()
method = "inverse.sum"
n.perm = 0
sig.quant = 0.8

pcor_est <- pcor.shrink(expr.mat)
test_pcor <- network.test.edges(pcor_est)


########### the values of lam1 were selected to make the results of different methods comparable.
##-----------------
### GeneNet
##-----------------
cat("Constructing using GeneNet ...")
invisible(capture.output(pcor_est <- pcor.shrink(expr.mat)))
invisible(capture.output(test_pcor <- network.test.edges(pcor_est)))
est_edge.GeneNet <- matrix(NA, p, p)

for (i in 1:nrow(test_pcor)) {
  est_edge.GeneNet[test_pcor$node1[i], test_pcor$node2[i]] <- test_pcor$prob[i]
}
est_edge[[1]] <- est_edge.GeneNet
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")

## plot network
edge.result0 <- est_edge[[1]]
edge.result0[which(edge.result0!=0)] <- 1
temp <- graph.adjacency(adjmatrix = edge.result0, mode = "undirected")
temp.degree <- apply(edge.result0 , 2, sum)
V(temp)$color <- (temp.degree>15)+3
plot(temp, vertex.size=5, vertex.frame.color="white",
     layout=layout.fruchterman.reingold, vertex.label=NA, edge.connectivity = NA)

# store prob
est_edge.prob <- test_pcor[test_pcor$node1 < test_pcor$node2, c(2,3,6)]
colnames(est_edge.prob) <- c("Node1", "Node2", "GeneNet")
est_edge.prob$GeneNet <- round(est_edge.prob$GeneNet, 5)
cat(" OK\n")

##-----------------------------
### 1. neighborhood selection
##-----------------------------
cat("Constructing using neighborhood selection ...")
est_res.ns <- matrix(0, p, p)
lambda <- 0.5
alpha <- (1-lambda)^2 * p^2
for (k in 1:p) {
  rsp <- expr.mat[, k]
  prd <- t(expr.mat[, -k])
  lam <- sqrt(sum(rsp ^ 2)) * qnorm(alpha / (2 * p ^ 2), lower.tail = F)
  result1 <- l2.reg(prd, rsp, lambda = lam)
  est_res.ns[k, -k] <- result1$estimate
}
est_edge[[2]] <- (abs(est_res.ns) + t(abs(est_res.ns))) / 2

## edges detected
fit.adj=abs(est_edge[[2]])>1e-6
sum(fit.adj==1)/2                     ##total number of edges detected

## plot network
edge.result1 <- est_edge[[2]]
edge.result1[which(edge.result1!=0)] <- 1
temp <- graph.adjacency(adjmatrix = edge.result1, mode = "undirected")
temp.degree <- apply(edge.result1 , 2, sum)
V(temp)$color <- (temp.degree>20)+3
plot(temp, vertex.size=5, vertex.frame.color="white",
     layout=layout.fruchterman.reingold, vertex.label=NA, edge.connectivity = NA)
# store prob
est_edge.prob <- cbind(est_edge.prob, NS = mapply(FUN = function(x,y) est_edge[[2]][x,y], 
                                                  est_edge.prob$Node1, est_edge.prob$Node2))

## select significant compound ##
hub.degree <- temp.degree[(which(temp.degree>20))]
hub.index <- which(temp.degree>20)
hub <- cbind(hub.degree,name.comp[hub.index])
hub <- hub[order(hub[,1],decreasing = TRUE),]
write.table(hub, file = paste0(WhereAmI,"result/2-neighborhood_selection.txt"), sep = "\t", row.names = FALSE,col.names = FALSE)
cat(" OK\n")


##------------------
### 2. glasso
##------------------
cat("Constructing using GLASSO ...")
expr.mat <- as.matrix(expr.mat)
S <- t(expr.mat) %*% expr.mat / n
result2 <- glasso(S, rho = 0.006)
est_edge[[3]] <- abs(result2$wi)


## edges detected
diag(est_edge[[3]]) <- 0 
fit.adj=abs(est_edge[[3]])>1e-6
sum(fit.adj==1)/2                     ##total number of edges detected

## plot network
edge.result2 <- est_edge[[3]]
diag(edge.result2) <- 0
edge.result2[which(edge.result2!=0)] <- 1
temp <- graph.adjacency(adjmatrix = edge.result2, mode = "undirected")
temp.degree <- apply(edge.result2 , 2, sum)
V(temp)$color <- (temp.degree>20)+3
plot(temp, vertex.size=5, vertex.frame.color="white",
     layout=layout.fruchterman.reingold, vertex.label=NA, edge.connectivity = NA)

# store prob
est_edge.prob <- cbind(est_edge.prob, NS = mapply(FUN = function(x,y) est_edge[[3]][x,y], 
                                                  est_edge.prob$Node1, est_edge.prob$Node2))

## select significant compound ##
hub.degree <- temp.degree[(which(temp.degree>20))]
hub.index <- which(temp.degree>20)
hub <- cbind(hub.degree,name.comp[hub.index])
hub <- hub[order(hub[,1],decreasing = TRUE),]
write.table(hub, file = paste0(WhereAmI,"result/3-glasso.txt"), sep = "\t", row.names = FALSE,col.names = FALSE)

cat(" OK\n")

##-------------------
### glasso-sf
##-------------------
cat("Constructing using GLASSO-SF ...")
result3 <- glasso_sf(expr.mat, alpha = 0.008)
est_edge[[4]] <- abs(result3$wi)

## edges detected
diag(est_edge[[4]]) <- 0 
fit.adj=abs(est_edge[[4]])>1e-6
sum(fit.adj==1)/2                     ##total number of edges detected


## plot network
edge.result3 <- est_edge[[4]]
diag(edge.result3) <- 0
edge.result3[which(edge.result3!=0)] <- 1
temp <- graph.adjacency(adjmatrix = edge.result3, mode = "undirected")
temp.degree <- apply(edge.result3 , 2, sum)
V(temp)$color <- (temp.degree>20)+3
plot(temp, vertex.size=5, vertex.frame.color="white",
     layout=layout.fruchterman.reingold, vertex.label=NA, edge.connectivity = NA)

# store prob
est_edge.prob <- cbind(est_edge.prob, NS = mapply(FUN = function(x,y) est_edge[[4]][x,y], 
                                                  est_edge.prob$Node1, est_edge.prob$Node2))

## select hub ##
hub.degree <- temp.degree[(which(temp.degree>20))]
hub.index <- which(temp.degree>20)
hub <- cbind(hub.degree,name.comp[hub.index])
hub <- hub[order(hub[,1],decreasing = TRUE),]

write.table(hub, file = paste0(WhereAmI,"result/4-glassosf.txt"), sep = "\t", row.names = FALSE,col.names = FALSE)

cat(" OK\n")
##---------------
### pcacmi
##---------------
cat("Constructing using PCACMI ...")
invisible(capture.output(result4 <- pca_cmi(t(expr.mat), 0.03)))
est_edge[[5]] <- abs(result4$Gval)

## edges detected
diag(est_edge[[5]]) <- 0 
fit.adj=abs(est_edge[[5]])>1e-6
sum(fit.adj==1)/2                     ##total number of edges detected
sum(fit.adj[true.adj==1]==1)/2        ##total number of true edges detected

## plot network
edge.result4 <- est_edge[[5]]
diag(edge.result4) <- 0
edge.result3[which(edge.result4!=0)] <- 1
temp <- graph.adjacency(adjmatrix = edge.result4, mode = "undirected")
temp.degree <- apply(edge.result4 , 2, sum)
V(temp)$color <- (temp.degree>4)+3
plot(temp, vertex.size=5, vertex.frame.color="white",
     layout=layout.fruchterman.reingold, vertex.label=NA, edge.connectivity = NA)

hub.degree <- temp.degree[(which(temp.degree>20))]
hub.index <- which(temp.degree>20)
hub <- cbind(hub.degree,name.comp[hub.index])
hub <- hub[order(hub[,1],decreasing = TRUE),]
write.table(hub, file = paste0(WhereAmI,"result/5-pcacmi.txt"), sep = "\t", row.names = FALSE,col.names = FALSE)


cat(" OK\n")


##------------------
### space
##------------------
cat("Constructing using SPACE ...")
invisible(capture.output(result5 <- space.joint(expr.mat, lam1 = 0.009 * n, iter = 5)))
est_edge[[6]] <- abs(result5$ParCor)

## edges detected
diag(est_edge[[6]]) <- 0 
fit.adj=abs(est_edge[[6]])>1e-6
sum(fit.adj==1)/2                     ##total number of edges detected

## plot network
edge.result5 <- est_edge[[6]]
diag(edge.result5) <- 0
edge.result5[which(edge.result5!=0)] <- 1
temp <- graph.adjacency(adjmatrix = edge.result5, mode = "undirected")
temp.degree <- apply(edge.result5 , 2, sum)
V(temp)$color <- (temp.degree>20)+3
plot(temp, vertex.size=5, vertex.frame.color="white",
     layout=layout.fruchterman.reingold, vertex.label=NA, edge.connectivity = NA)

## select significant compound ##
hub.degree <- temp.degree[(which(temp.degree>20))]
hub.index <- which(temp.degree>20)
hub <- cbind(hub.degree,name.comp[hub.index])
hub <- hub[order(hub[,1],decreasing = TRUE),]
write.table(hub, file = paste0(WhereAmI,"result/6-space.txt"), sep = "\t", row.names = FALSE,col.names = FALSE)

cat(" OK\n")
