
#Generate cluster data
#https://search.r-project.org/CRAN/refmans/fungible/html/monte.html

#The following script is designed for large dimensional data with nodes of 200
#Irene Zeng April 2024

####larger cluster number of cluster 200, N=49530-subjects with 200 cluster, for Node=50 
#n=300, p=50, number of cluster 200, N=39800

library(fungible)
corarray=array(dim=c(50,50,200)) 


randcorr=function(d,k)
 { 
   W = matrix(rnorm(k*d,0,1),k)
   S = t(W)%*%W + diag(runif(d))
   SS = diag(1/sqrt(diag(S))) %*% S %*% diag(1/sqrt(diag(S)))
return(SS)
 }

set.seed(100)
for (i in 1:200)
 {
 corarray[,,i]=randcorr(d=50,k=200)
 }

list1 <- list()
# fill the list with arrays
for (i in 1:200) {
  list1[[i]] <- corarray[,,i]
}
cormatlist=list1

cluster=c(1:200)
set.seed(10)
clusize=round(rnorm(n=200,mean=5,sd=1),0)
clusize
sum(clusize)
#[1] 987


clusize_pop=clusize*50
sum(clusize_pop)
#[1] 49350 

popclus<-monte(seed = 123, nvar = 50, nclus = 200, 
               cor.list = cormatlist, clus.size = clusize_pop,
               eta2 = runif(200))

popclusdata=data.frame(popclus$data)
dim(popclusdata)
[1] 49350    51


popdata_clus=as.matrix(popclusdata[,2:51])
pop_graph=huge(popdata_clus)
pop_graph_adj_5=as.matrix(pop_graph$path[[5]])
pop_graph$lambda
[1] 0.97422111 0.75430402 0.58403021 0.45219338 0.35011691 0.27108281
 [7] 0.20988957 0.16250988 0.12582549 0.09742211



#####################################################node p=200
#first set of cluster 20
library(huge)
cluster=c(1:20)
set.seed(10)
clusize=round(rnorm(n=20,mean=10,sd=1),0)
clusize
sum(clusize)
#[1] 199

clusize_pop=clusize*200
sum(clusize_pop)
#[1] 39800


library(fungible)
randcorr=function(d,k)
 { 
   W = matrix(rnorm(k*d,0,1),k)
   S = t(W)%*%W + diag(runif(d))
   SS = diag(1/sqrt(diag(S))) %*% S %*% diag(1/sqrt(diag(S)))
return(SS)
 }

corarray=array(dim=c(200,200,20)) 
set.seed(100)
for (i in 1:20)
 {
 corarray[,,i]=randcorr(d=200,k=20)
 }

list1 <- list()
# fill the list with arrays
for (i in 1:20) {
  list1[[i]] <- corarray[,,i]
}
cormatlist=list1

popclus<-monte(seed = 123, nvar = 200, nclus = 20, 
               cor.list = cormatlist, clus.size = clusize_pop,
               eta2 = runif(200))

popclusdata=data.frame(popclus$data)
dim(popclusdata)
[1] 39800   201

#data use the monte function to generate cluste data popclusdata 
popdata_clus=as.matrix(popclusdata[,2:201])
pop_graph=huge(popdata_clus)
pop_graph_adj_5=as.matrix(pop_graph$path[[5]])
pop_graph$lambda
 [1] 0.98253667 0.76074246 0.58901526 0.45605312 0.35310537 0.27339666
 [7] 0.21168111 0.16389699 0.12689949 0.09825367

####larger cluster number of cluster 200, N=49530-subjects with 200 cluster, for Node=200 

library(fungible)
corarray=array(dim=c(200,200,200)) 


randcorr=function(d,k)
 { 
   W = matrix(rnorm(k*d,0,1),k)
   S = t(W)%*%W + diag(runif(d))
   SS = diag(1/sqrt(diag(S))) %*% S %*% diag(1/sqrt(diag(S)))
return(SS)
 }

set.seed(100)
for (i in 1:200)
 {
 corarray[,,i]=randcorr(d=200,k=200)
 }

list1 <- list()
# fill the list with arrays
for (i in 1:200) {
  list1[[i]] <- corarray[,,i]
}
cormatlist=list1

cluster=c(1:200)
set.seed(10)
clusize=round(rnorm(n=200,mean=5,sd=1),0)
clusize
sum(clusize)
#[1] 987


clusize_pop=clusize*50
sum(clusize_pop)
#[1] 49350 

popclus<-monte(seed = 123, nvar = 200, nclus = 200, 
               cor.list = cormatlist, clus.size = clusize_pop,
               eta2 = runif(200))

popclusdata=data.frame(popclus$data)
dim(popclusdata)
[1] 49350   201

popdata_clus=as.matrix(popclusdata[,2:201])
pop_graph=huge(popdata_clus)
pop_graph_adj_5=as.matrix(pop_graph$path[[5]])
pop_graph$lambda

 [1] 0.98090810 0.75948152 0.58803896 0.45529721 0.35252009 0.27294351
 [7] 0.21133024 0.16362533 0.12668915 0.09809081


