library(PPMR)

#load the scaled genenotype matrix in eQTL data (e.g. cis-SNPs of BACE1 gene from GEUVADIS data)
zx<-read.table("zx.txt")
zx<-as.matrix(zx)
n1 = dim(zx)[1]
q = dim(zx)[2]

#load the scaled genenotype matrix in GWAS data (e.g. the same cis-SNPs from GERA data)
zy<-read.table("zy.txt")
zy<-as.matrix(zy)
n2<-dim(zy)[1]

#set PVE_zx to be 10%
squaresigma_beta<-0.1/q
squaresigma_x<-0.9

#set the common pleiotropy effect to be 0.001
gamma=0.001

#set the causal effect to be 0
alpha=0

#get the simulated gene expression data
beta <- matrix(rnorm(q, 0, sd = sqrt(squaresigma_beta)),q,1)
epison_x <- matrix(rnorm(n1, 0, sd = sqrt(squaresigma_x)), n1, 1)
x <- zx %*% beta + epison_x
x<-as.matrix(x)

#get the simulated phenotype
y_mean=as.vector(zy%*%rep(betaa,q))
squaresigma_y<-1-alpha^2
epison_y<-matrix(rnorm(n2, 0, sd = sqrt(squaresigma_y)), n2, 1)
y<-y_mean+epison_y

#run the PPMR model using PMR_individual function
fmH1 = PMR_individual(x, y, zx, zy,gammain=0,alphain = 0,max_iterin =1000,epsin=1e-5)
fmH0gamma = PMR_individual(x, y, zx, zy,gammain=1, alphain = 0,max_iterin =1000,epsin=1e-5)
fmH0alpha = PMR_individual(x, y, zx, zy,gammain=0,alphain = 1,max_iterin =1000, epsin=1e-5)
loglikH1=max(fmH1$loglik,na.rm=T)
loglikH0gamma=max(fmH0gamma$loglik,na.rm=T)
loglikH0alpha=max(fmH0alpha$loglik,na.rm=T)
stat_alpha = 2 * (loglikH1 - loglikH0alpha)
stat_gamma = 2 * (loglikH1 - loglikH0gamma)

#get the estimate of the causal effect
alpha<-fmH1$alpha

#get the estimate of the pleiotropy effect
gamma<-fmH1$gamma

#get the pvalue for the causal test
pvalue_alpha = pchisq(stat_alpha,1,lower.tail=F)

#get the pvalue for the pleiotropy effect
pvalue_gamma = pchisq(stat_gamma,1,lower.tail=F)
#################################
###################################
