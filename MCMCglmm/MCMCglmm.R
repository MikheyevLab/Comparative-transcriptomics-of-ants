
```{r}
library(ape)
library(caper)
library(geiger)
library(phytools)
library(MCMCglmm)
library(coda)
library(batchmeans)
```


```{r}
#Data
phylo1<- read.tree(file="ants.tre")

data = read.table("data.txt", header =T)
attach(data)

#Cross validation¨
unrooted_tr <- unroot(phylo1)
l <- 10^(-1:6)
cv <- numeric(length(l))

cv <- sapply(l, function(x) sum(attr(chronopl(unrooted_tr, lambda=1, CV=TRUE), "D2")))
plot(l,cv)

#Ultrametric tree
chronopl(phylo1, lambda=0.1) -> phylo

#Inverse tree to include phylogeny in the analysis
inv.phylo<-inverseA(phylo,nodes="TIPS")

#priors 3 different ones
prior<-list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))
```

```{r}
#Module 1
model_1<-MCMCglmm(MEmodule1~ Caste + Invasiveness + Sterility + Queen_number,random=~phylo,
                 ginverse=list(phylo=inv.phylo$Ainv),prior=prior,data=data,nitt=500000,burnin=150000,thin=500,verbose = FALSE)
summary(model_1)

model_2<-MCMCglmm(MEmodule1~ Caste + Invasiveness + Sterility + Queen_number,random=~phylo,
                 ginverse=list(phylo=inv.phylo$Ainv),prior=prior,data=data,nitt=500000,burnin=150000,thin=500,verbose = FALSE)
summary(model_2)

model_3<-MCMCglmm(MEmodule1~ Caste + Invasiveness + Sterility + Queen_number,random=~phylo,
                 ginverse=list(phylo=inv.phylo$Ainv),prior=prior,data=data,nitt=500000,burnin=150000,thin=500,verbose = FALSE)
summary(model_3)

chainList_DG <- mcmc.list(model_1$Sol,model_2$Sol,model_3$Sol)
gelman.diag(chainList_DG)
effectiveSize(chainList_DG)
plot(chainList_DG)
summary(chainList_DG)
```
