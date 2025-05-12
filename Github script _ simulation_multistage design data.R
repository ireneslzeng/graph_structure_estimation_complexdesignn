
#Generate multistage design data 
#Irene Zeng Feb, 2024

#experiment 1-node=50
library(huge)
pop1= huge.generator(n=5248928,d=50,graph="random",vis=TRUE)
data_pop1=pop1$data
pop1_graph=huge(data_pop1)
pop1_graph_adj_5=as.matrix(pop1_graph$path[[5]])
pop1_graph$lambda
[1] 0.48436137 0.37502342 0.29036701 0.22482063 0.17407045 0.13477643
[7] 0.10435249 0.08079635 0.06255768 0.04843614

set.seed(100)
mbsize=round(rnorm(n=100,mean=700,sd=5),0)
mbsize
sum(mbsize)
#[1] 70002

set.seed(99)
psusize=round(rnorm(n=70002,mean=25,sd=1),0)
sum(psusize)
#[1] 1749480

set.seed(80)
housesize=round(rnorm(n=1749480,mean=3,sd=0.5),0)
sum(housesize)
#[1] 5248928

personid=c(1:5248928)
person_housenum=vector(mode="numeric",length=5248928)
person_mbnum=vector(mode="numeric",length=5248928)
person_aunum=vector(mode="numeric",length=5248928)

aunum=c(1:100)
mbnum=c(1:70002)
housenum=c(1:1749480)

m=1
for (k in 1:1749480)
{t=1 
while (t<=housesize[k]) 
{person_housenum[m]=k
t=t+1; m=m+1}
}

m=1
for (k in 1:70002)
{t=1 
while (t<=psusize[k]) 
{person_mbnum[m]=k
t=t+1; m=m+1}
}
m=1
for (k in 1:100)
{t=1 
while (t<=mbsize[k]) 
{person_aunum[m]=k
t=t+1; m=m+1}
}

sampleindex_1=data.frame(personid,person_housenum,person_mbnum,person_aunum)

#Sampling and weight generation 
#Adding nonresponses rate varied by GEO units?
#unequal probablity of sampling across AU
#nonresponse & poststratified callibration
p1=c(0.2,0.1,0.1,0.05,0.075,rep(0.005,95))
#within each AU, 4 stratas were formed according to the population subgroups (rural/urban, ethnicity) and each AU 
#has different selection probability  
set.seed(50)
#p2=rgeom(70002,prob=c(0.1,0.4,0.1,0.4))/70002 				#prob of selection for each strata within the AU-now each MB also had a selection prob
p2=rbinom(70002,100,prob=c(0.1,0.4,0.1,0.4))/100 			#prob of MB selected

f1=1/p1   										#adjusting factor for AU in sampling
prob_adjsele=vector(mode="numeric",length=70002)
sum_f=sum(f1)

m=1
for (k in 1:100)
{t=1
while (t<=mbsize[k]) {prob_adjsele[m]=(p2[m]*p1[k]);t=t+1;m=m+1} #PPS method 
}

#Scale to between 0-100;
prob_adjsele_scaled=prob_adjsele/sum(prob_adjsele)*100
weight=1/prob_adjsele_scaled
