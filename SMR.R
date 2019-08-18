##suppose you have got the following summary statistics#####
######betax for eQTL effects vector,sigmax for the standard error of betax
####betay for GWAS effects vector,sigmay for the standard error of betay

Z1<-(betax/sigmax)^2
newindex<-which(Z1==max(Z1))
index<-newindex[1]
Z2<-(betay/sigmay)^2
TSMR<-(Z1[index]*Z2[index])/(Z1[index]+Z2[index])
pvalue<-1-pchisq(TSMR,df=1)
