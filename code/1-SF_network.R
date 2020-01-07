########## import library #############
library(R.matlab)
library(space)
library(igraph)
library(glasso)
library(poweRlaw)
library(pheatmap)
WhereAmI <- "D:/Ziyi/School/PMO/13.Network/1-code/4-ziyi/"

source(paste0(WhereAmI,"lib/glasso_SF.R"))
source(paste0(WhereAmI,"lib/BayesianGLasso.R"))
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
V(temp)$color=(temp.degree>6)+3
plot(temp, vertex.size=5, vertex.frame.color="white",layout=layout.fruchterman.reingold, vertex.label=NA, edge.connectivity = NA)


#################### estimate the parcial correlation matrix with various methods
alpha=1
l1=1/sqrt(n)*qnorm(1-alpha/(2*p^2))
iter=6

########### the values of lam1 were selected to make the results of different methods comparable.
#### 1. MB method
result1=space.neighbor(sample, lam1=l1*0.3, lam2=0)
fit.adj = abs(result1$ParCor)>1e-6
sum(fit.adj==1)/2                      ##total number of edges detected
sum(fit.adj[true.adj==1]==1)/2         ##total number of true edges detected

## plot network
edge.result1 <- result1$ParCor
diag(edge.result1) <- 0
edge.result1[which(edge.result1!=0)] <- 1
temp <- graph.adjacency(adjmatrix = edge.result1, mode = "undirected")
temp.degree <- apply(edge.result1 , 2, sum)
V(temp)$color <- (temp.degree>10)+3
plot(temp, vertex.size=5, vertex.frame.color="white",
     layout=layout.fruchterman.reingold, vertex.label=NA, edge.connectivity = NA)


#### 2. JSRM
result2=space.joint(sample, lam1=l1*n*0.9, lam2=0, iter=iter)

fit.adj=abs(result2$ParCor)>1e-6
sum(fit.adj==1)/2                      ##total number of edges detected
sum(fit.adj[true.adj==1]==1)/2         ##total number of true edges detected

## plot network
edge.result2 <- result2$ParCor
diag(edge.result2) <- 0
edge.result2[which(edge.result2!=0)] <- 1
temp <- graph.adjacency(adjmatrix = edge.result2, mode = "undirected")
temp.degree <- apply(edge.result2 , 2, sum)
V(temp)$color <- (temp.degree>8)+3
plot(temp, vertex.size=5, vertex.frame.color="white",
     layout=layout.fruchterman.reingold, vertex.label=NA, edge.connectivity = NA)


#### 3. JSRM.dew
result3=space.joint(sample, lam1=l1*n*0.9, lam2=0, weight=2, iter=iter)

fit.adj=abs(result3$ParCor)>1e-6
sum(fit.adj==1)/2                      ##total number of edges detected
sum(fit.adj[true.adj==1]==1)/2         ##total number of true edges detected


## plot network
edge.result3 <- result3$ParCor
diag(edge.result3) <- 0
edge.result3[which(edge.result3!=0)] <- 1
temp <- graph.adjacency(adjmatrix = edge.result3, mode = "undirected")
temp.degree <- apply(edge.result3 , 2, sum)
V(temp)$color <- (temp.degree>10)+3
plot(temp, vertex.size=5, vertex.frame.color="white",
     layout=layout.fruchterman.reingold, vertex.label=NA, edge.connectivity = NA)


#### 4. glasso
result4 <- glasso(emp_cov, rho = 0.2)

fit.adj=abs(result4$wi)>1e-6
sum(fit.adj==1)/2                      ##total number of edges detected
sum(fit.adj[true.adj==1]==1)/2         ##total number of true edges detected

## plot network
edge.result4 <- result4$wi
diag(edge.result4) <- 0
edge.result4[which(edge.result4!=0)] <- 1
temp <- graph.adjacency(adjmatrix = edge.result4, mode = "undirected")
temp.degree <- apply(edge.result4 , 2, sum)
V(temp)$color <- (temp.degree>10)+3
plot(temp, vertex.size=5, vertex.frame.color="white",
     layout=layout.fruchterman.reingold, vertex.label=NA, edge.connectivity = NA)



#### 5. glasso-SF
result5 <- glasso_sf(sample, alpha = 0.16)

fit.adj=abs(result5$wi)>1e-6
sum(fit.adj==1)/2                     ##total number of edges detected
sum(fit.adj[true.adj==1]==1)/2        ##total number of true edges detected


## plot network
edge.result5 <- result5$wi
diag(edge.result5) <- 0
edge.result5[which(edge.result5!=0)] <- 1
temp <- graph.adjacency(adjmatrix = edge.result5, mode = "undirected")
temp.degree <- apply(edge.result5 , 2, sum)
V(temp)$color <- (temp.degree>10)+3
plot(temp, vertex.size=5, vertex.frame.color="white",
     layout=layout.fruchterman.reingold, vertex.label=NA, edge.connectivity = NA)

## heatmap
pheatmap(theta, cluster_rows = FALSE, cluster_cols = FALSE)
data1 <- result5$wi
pheatmap(data1, cluster_rows = FALSE, cluster_cols = FALSE)