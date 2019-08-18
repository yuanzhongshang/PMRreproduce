library(CoMM)
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


#run CoMM model###
w2 = matrix(rep(1,n2),ncol=1)
w1 = matrix(rep(1,n1),ncol=1)
H1 = CoMM_covar_pxem(x, y, zx, zy, w1, w2,constr = 0)
H0 = CoMM_covar_pxem(x, y, zx, zy, w1, w2,constr = 1)
loglikH1 = max(H1$loglik,na.rm=T)
loglikH0 = max(H0$loglik,na.rm=T)
stat = 2 * (loglikH1 - loglikH0)
pvalue = pchisq(stat,1,lower.tail=F)

#################################
###################################