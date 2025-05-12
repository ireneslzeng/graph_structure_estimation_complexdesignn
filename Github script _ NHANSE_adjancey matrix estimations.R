#The following script is for network analysis in NHANSE
#Irene Zeng Feb 2024
library(huge)
library(igraph)
Umarkers=read.csv(file.choose(),header=T,sep=",")
attach(Umarkers)

Umatrix=as.matrix(cbind(URXUBA,URXUCD,URXUCO,URXUCS,URXUMN,URXUMO,URXUPB,URXUSB,URXUSN,URXUTL,URXUTU,
                        URXUNI,URXUHG,URXUIO,URXUAB,URXUAC,URXUAS3,URXUAS5,URXUDMA,URXUMMA,URXUAS))
nm=c("BA","CD","CO","CS","MN","MO","PB","SB","SN","TL","TU","NI","HG","IO","AB","AC","AS3","AS5","DMA","MMA","AS")
############################################################Construct different age cohort on the original scale########################
Umarkers_cohort1=subset(Umarkers2,((Umarkers2$RIDAGEYR>=0)&(Umarkers2$RIDAGEYR < 10)))

Umarkers_cohort2=subset(Umarkers2,((Umarkers2$RIDAGEYR>=10)&(Umarkers2$RIDAGEYR < 20)))

Umarkers_cohort3=subset(Umarkers2,((Umarkers2$RIDAGEYR>=20) & (Umarkers2$RIDAGEYR < 65)))

Umarkers_cohort4=subset(Umarkers2,((Umarkers2$RIDAGEYR>=65)))
Umatrix_cohort1=na.omit(Umatrix[((RIDAGEYR>=0) & (RIDAGEYR < 10)),])
dim(Umatrix_cohort1)
#[1] 376  21

Umatrix_cohort2=na.omit(Umatrix[((RIDAGEYR>=10) & (RIDAGEYR < 20)),])
dim(Umatrix_cohort2)
#[1] 788  21

Umatrix_cohort3=na.omit(Umatrix[((RIDAGEYR>=20) & (RIDAGEYR < 65)),])
dim(Umatrix_cohort3) 
#[1] 2052   21

Umatrix_cohort4=na.omit(Umatrix[((RIDAGEYR>=65)),])
dim(Umatrix_cohort4)
#[1] 693  21

#cohort 1- on original scale
#adjusted by survey sampling weights
library(Matrix)

###Weighted analysis
cohort1w=subset(Umarkers_cohort1,(!(is.na(Umarkers_cohort1$WTSAPRP))&(Umarkers_cohort1$WTSAPRP>0)))
cohort1_matrix=with(as.matrix(cbind(URXUBA,URXUCD,URXUCO,URXUCS,URXUMN,URXUMO,URXUPB,URXUSB,URXUSN,URXUTL,URXUTU,
                                    URXUNI,URXUHG,URXUIO,URXUAB,URXUAC,URXUAS3,URXUAS5,URXUDMA,URXUMMA,URXUAS)),data=cohort1w)

cohort1w_2=na.omit(subset(cohort1w,select=c(SDMVSTRA,WTSAPRP,URXUBA,URXUCD,URXUCO,URXUCS,URXUMN,URXUMO,URXUPB,URXUSB,URXUSN,URXUTL,URXUTU,
                                            URXUNI,URXUHG,URXUIO,URXUAB,URXUAC,URXUAS3,URXUAS5,URXUDMA,URXUMMA,URXUAS)))


strata1=round(cohort1w$WTSAPRP,digit=0)
nstrata1=length(table(strata1))
prob_select_1=strata1/sum(na.omit(strata1))
npop_1=sum(na.omit(strata1))
#> npop_1
#[1] 28536154

population_1=vector(mode="numeric",length=npop_1)
id_1=c(1:dim(cohort1w)[1])


#selection functions for adjencey function to be presented#
select=function (x) {ifelse (x>0.5,1,0)}
select2=function (x) {ifelse (x==1,1,0)}

##psudo-population estimate 
corr_sample_1_adj=array(dim=c(21,21,10,1000))
eigenvalueM_p1=array(dim=c(1000,10))

for (rep in 1:1000)
{
  corrsample=sample(id_1,size=npop_1,replace=T,prob=prob_select_1)
  corr_samples=na.omit(cohort1_matrix[corrsample,])
  Q3=huge.npn(corr_samples, npn.func = "skeptic")
  out.mb_trans3=huge(Q3, method = "mb")
  for (j in 1:10) {corr_sample_1_adj[,,j,rep]=as.matrix(out.mb_trans3$path[[j]])
  eigenvalueM_p1[rep,j]=spectrum(graph_from_adjacency_matrix(out.mb_trans$path[[j]],weighted=TRUE,mode="undirected"),which = list(pos="LA",howmany=1))$value
  }
}

weigh_adj=apply(corr_sample_1_adj[,,5,1:100],1:2,mean)  #fix the path with penalized parameter of [5] solution
select=function (x) {ifelse (x>0.5,1,0)}
weigh_adj2=apply(weigh_adj,1:2,select)  			  #selecte only replicate average p>0.5

select2=function (x) {ifelse (x==1,1,0)}
weigh_adj3=apply(weigh_adj,1:2,select2) 			  #selecte only replicate average p=1

write.csv(weigh_adj2,file="pathname/weigh_adj2.csv")
write.csv(weigh_adj3,file="pathname/weigh_adj3.csv")
write.csv(eigenvalueM_p1,file="pathname/eigenvaluematrix_pop_cohort1.csv")


##Use BRR estimate of the adjacency matrix##


#cohort 1#
nsize1=dim(cohort1w)[1]
id_1=c(1:dim(cohort1w)[1])

half_sample_1_adj=array(dim=c(21,21,10,1000))
eigenvalueM=array(dim=c(1000,10))
for (rep in 1:1000)
{
  halfsample=sample(id_1,size=nsize1/2,replace=T)
  corr_samples=na.omit(cohort1_matrix[halfsample,])
  Q3=huge.npn(corr_samples, npn.func = "skeptic")
  out.mb_trans=huge(Q3, method = "mb")
  for (j in 1:10) {half_sample_1_adj[,,j,rep]=as.matrix(out.mb_trans$path[[j]])
  eigenvalueM[rep,j]=spectrum(graph_from_adjacency_matrix(out.mb_trans$path[[j]],weighted=TRUE,mode="undirected"),which = list(pos="LA",howmany=1))$value
  }
}

half_adj_co1=apply(half_sample_1_adj[,,5,1:1000],1:2,mean)      #fix the path with penalized parameter of [5] solution
half_adj2_co1=apply(half_adj_co1,1:2,select)  			   #selecte only replicate average p>0.5
half_adj3_co1=apply(half_adj_co1,1:2,select2) 			   #selecte only replicate average p=1


write.csv(half_adj2_co1,file="pathname/half_adj2_cohort1.csv")
write.csv(half_adj3_co1,file="pathname/half_adj3_cohort1.csv")
write.csv(eigenvalueM,file="pathname/eigenvaluematrix_cohort1.csv")



#Including vertix names
colnames(half_adj2_co1)=nm
colnames(half_adj3_co1)=nm

g3=graph_from_adjacency_matrix(half_adj2_co1, weighted=TRUE,mode="undirected",add.colnames="name")
windows(7,7)
plot(g3, main="Survey sample using BRR estimated network -cohort 1 (>0.5)")#the middle solution

g4=graph_from_adjacency_matrix(half_adj3_co1, weighted=TRUE,mode="undirected",add.colnames="name")
windows(7,7)
plot(g4, main="Survey sample usin BRR estimated network- cohort 1 (=1)")#the middle solution



##Use JKR estimate of the adjacency matrix##
jk_sample_1_adj=array(dim=c(21,21,10,1000))
eigenvalueM=array(dim=c(1000,10))

/*leave 10% out, resampling without replacment*/
  
  for (rep in 1:1000)
  {
    jksample=sample(id_1,size=ceiling(nsize1*0.9),replace=F)
    corr_samples=na.omit(cohort1_matrix[jksample,])
    Q3=huge.npn(corr_samples, npn.func = "skeptic")
    out.mb_trans=huge(Q3, method = "mb")
    for (j in 1:10) {jk_sample_1_adj[,,j,rep]=as.matrix(out.mb_trans$path[[j]])
    eigenvalueM[rep,j]=spectrum(graph_from_adjacency_matrix(out.mb_trans$path[[j]],weighted=TRUE,mode="undirected"),which = list(pos="LA",howmany=1))$value
    }
  }

jk_adj_co1=apply(jk_sample_1_adj[,,5,1:1000],1:2,mean)      #fix the path with penalized parameter of [5] solution
jk_adj2_co1=apply(jk_adj_co1,1:2,select)  			   		      #selecte only replicate average p>0.5
jk_adj3_co1=apply(jk_adj_co1,1:2,select2) 			   		      #selecte only replicate average p=1


##Bootstrapping approach##
nsample_1=dim(cohort1_matrix)[1]
corr_sample_1_adj=array(dim=c(21,21,10,1000))
for (rep in 1:1000)
{
  corrsample=sample(id_1,nsample_1,replace=T,prob=prob_select_1)
  corr_samples=na.omit(cohort1_matrix[corrsample,])
  Q3=huge.npn(corr_samples, npn.func = "skeptic")
  out.mb_trans3=huge(Q3, method = "mb")
  for (j in 1:10) corr_sample_1_adj[,,j,rep]=as.matrix(out.mb_trans3$path[[j]])
}

weigh_adj=apply(corr_sample_1_adj[,,5,1:1000],1:2,mean)  #fix the path with penalized parameter of [5] solution

select=function (x) {ifelse (x>0.5,1,0)}
weigh_adj2=apply(weigh_adj,1:2,select)  #selecte only replicate average >0.5

select2=function (x) {ifelse (x==1,1,0)}
weigh_adj3=apply(weigh_adj,1:2,select2) #selecte only replicate average =1

#Including vertix names
colnames(weigh_adj2)=nm
colnames(weigh_adj3)=nm

g3=graph_from_adjacency_matrix(weigh_adj2, weighted=TRUE,mode="undirected",add.colnames="name")
windows(7,7)
plot(g3, main="Survey sample weighted network -cohort 1 (>0.5)")#the middle solution

g4=graph_from_adjacency_matrix(weigh_adj3, weighted=TRUE,mode="undirected",add.colnames="name")
windows(7,7)
plot(g4, main="Survey sample weighted network- cohort 1 (=1)")#the middle solution



