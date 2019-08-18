##load the following R fucntions
##Thanks the authors of LDA MR-Egger paper to share the main LDA.MREgger function. 
inversesigma<-function(sigmaG){
  sigmaG<-as.matrix(sigmaG)
  tmp=sigma*(sigmaG%*%t(sigmaG))
  tmp<-0.99*tmp+0.01*diag(nrow(tmp))
  solvetmp<-solve(tmp)
  re<-sigma%*%solvetmp%*%sigma
  return(re)
}


LDA.MREgger<-function(X,Y,W){
  bX<-cbind(1,X)
  bread<-solve(crossprod(bX,W)%*%bX)
  theEsts<-bread%*%crossprod(bX,W%*%Y)
  theresid<-c(Y-theEsts[1]-X*theEsts[2])
  Sig.Est<-c(crossprod(theresid,W%*%theresid))/(length(X)-2)
  finresults<- cbind(theEsts,diag(bread)*Sig.Est)
  TestStat<-theEsts/sqrt(finresults[,2])
  Pvals<-2*pt(abs(TestStat),df = nrow(bX)-2,lower.tail = F)
  return(cbind(finresults,TestStat,Pvals))
}
##suppose you have got the following summary statistics#####
######betax for eQTL effects vector,sigmax for the standard error of betax
####betay for GWAS effects vector,sigmay for the standard error of betay

#load the scaled genenotype matrix in eQTL data (e.g. cis-SNPs of BACE1 gene from GEUVADIS data)
zx<-read.table("zx.txt")
zx<-as.matrix(zx)
n1 = dim(zx)[1]
q = dim(zx)[2]
###################################

aa<-t(zx)%*%(zx)/(n1-1)
aa<-0.99*aa+0.01*diag(nrow(aa))
LDsolvesigma<-solve(aa)
sigma<-aa

###get the conditional estimate for eQTL effect 
betaE<-as.vector(LDsolvesigma%*%betax)

##get the conditional estimate for GWAS effect
betaG<-as.vector(LDsolvesigma%*%betay)

solvesigma<-inversesigma(sigmay)
Z<-LDA.MREgger(betaE,betaG,solvesigma)

pvalue_gamma<-Z[1,4]
pvalue_alpha<-Z[2,4]
#################################
########################################

