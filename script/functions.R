
#FUNCTION TO SIMULATE ORDINAL DATA
#n= number of observations
#probs = list of marginal probabilities
# Cor=correlation matrix
#showCor_norm= show the value of the correlation matrix
#comes from the orddata package install.packages("orddata", repos="http://R-Forge.R-project.org")

#c<-rmvord(n=2,probs=probs1,Cor=Cor1)

rmvord<-function(n=1,probs,Cor,showCor_norm=FALSE)
{
  q=length(probs)
  categ_probs=0
  cumul_probs=list(0)
  quant_probs=list(0)
  means=0
  vars=0
  var.wt=function(x,w)
  {
    m=weighted.mean(x=x,w=w)
    sum((x[1:length(x)]-m)^2*w[1:length(x)])
  }
  for (i in 1:q)
  {
    categ_probs[i]=length(probs[[i]])
    cumul_probs[[i]]=cumsum(1:categ_probs[i]/10^12+probs[[i]])
    cumul_probs[[i]][categ_probs[i]]=1
    quant_probs[[i]]=qnorm(p=cumul_probs[[i]],mean=0,sd=1)
    means[i]=weighted.mean(x=1:categ_probs[i],w=probs[[i]])
    vars[i]=var.wt(x=1:categ_probs[i],w=probs[[i]])
  }
  Cor_norm=Cor
  for (i in 1:(q-1))
  {
    for (j in (i+1):q)
    {
      gridd=rep(0,times=201)
      for (ci in 1:(categ_probs[i]-1))
      {
        for (cj in 1:(categ_probs[j]-1))
        {
          for (steps in -100:100)
          {
            gridd[101+steps]=gridd[101+steps]+pmvnorm(upper=c(quant_probs[[i]][ci],quant_probs[[j]][cj]),corr=matrix(2,2,data=c(1,steps/100,steps/100,1)))[1]
          }
        }
      }
      f=suppressWarnings(approxfun(y=-100:100/100,x=gridd))
      Cor_norm[i,j]=Cor_norm[j,i]=f(Cor[i,j]*sqrt(vars[i]*vars[j])+means[i]*means[j]-categ_probs[i]*categ_probs[j]+categ_probs[j]*sum(cumul_probs[[i]][1:(categ_probs[i]-1)])+categ_probs[i]*sum(cumul_probs[[j]][1:(categ_probs[j]-1)]))
    }
  }
  if(showCor_norm) print(Cor_norm)
  retval=rmvnorm(n=n,sigma=Cor_norm)
  for (i in 1:q)
  {
    retval[,i]=cut(x=retval[,i],breaks=c(-1/0,quant_probs[[i]]),right=FALSE)
  }
  retval
}



#' Computation of the weighted kappa coefficient for more than two raters
#'
#' This function computes a weighted kappa coefficient in the presence of two or more observers and its standard error using the delta method. Note that this function requires a dataset with one rating per patient/object and observer.
#'
#' @param data: NxR matrix. dataset with the R ratings for the N participants coded coded from 1 to K, where K is the number of categories of the scale.
#' @param weight: weighting scheme to be used, "binary" for equal disagreement weights, "linear" for linear disagreement weights and "quadratic" for quadratic disagreement weights
#' @param chance: chance definition to be used, 1: each observer is assumed to classify the participants by rolling a fair dice, 2: the probability to classify participants in the K categories is assumed to be different for each observer, 3: the probabity to classify participants in the K categories is assumed to be the same for all observers (estimated as the average over all observers)
#' @return a vector where the first element is the weighted kappa coefficient and the second its standard error computed by the delta method
#' @author Sophie Vanbelle \email{sophie.vanbelle@@maastrichtuniversity.nl}
#' @details This function computes weighted kappa coefficient and its standard error using the delta method.    
#' @references H. J. A. Schouten. Measuring pairwise agreement among many observers. Biometrical Journal, 22(6):497–504, 1980.
#' @references H. J. A. Schouten. Measuring pairwise interobserver agreement when all subjects are judged by the same observers. Statistica Neerlandica, 36:45–61, 1982
#' @references Vanbelle S. et al. (submitted). Agreement between several observers on a ordinal scale
#' @export
#' @examples
#'
#' #dataset from Ayres de Campos
#'
#' data(CTG)
#' kappaw(CTG,"quadratic",2)
#'

kappaw<-function(data,weight,chance)
{
  
  n<-nrow(data)
  nrat<-ncol(data)
  npair<-nrat*(nrat-1)/2
  if (min(data==0)){data<-data+1}
  all<-c(data.matrix(data))
  all.f<-factor(all,levels=names(table(all)))
  ncat<-nlevels(all.f)    		#number of categories of the scale
  score <- matrix(1:ncat, nrow = ncat, ncol = ncat, byrow = TRUE)
  
  if (weight=="binary"){
    w<-diag(ncat)
  }
  if (weight=="linear"){
    w<-1-abs(score - t(score))/(ncat-1)
  }
  if (weight=="quadratic"){
    w<-1-(score - t(score))^2/(ncat-1)^2
  }
  
  freq<-t(sapply(1:n,fun2<-function(i){ table(factor(data[i,],levels=names(table(all))))}))
  
  p<-matrix(NA,nrow=ncat,ncol=ncat)
  oprim<-rep(NA,n)
  for (h in 1:n){
    for (i in 1:ncat){
      for (j in 1:ncat){
        if (i==j){
          p[i,j]<-w[i,j]*freq[h,i]*(freq[h,i]-1)/(nrat*(nrat-1))
        }
        if (i!=j){
          p[i,j]<-w[i,j]*freq[h,i]*freq[h,j]/(nrat*(nrat-1))
        }
      }
      
    }
    oprim[h]<-sum(p)
  }
  
  P_o<-mean(oprim)
  
  
  if (chance==1){
    #chance=1
    eprim<-rep(NA,n)
    p_j<-rep(1/ncat,ncat)
    eprim_tmp<-matrix(NA,nrow=n,ncol=nrat)
    eprim<-rep(NA,n)
    
    for (h in 1:n){
      for (r in 1:nrat){
        temp<-rep(NA,nrat-1)
        l<-1
        for (s in setdiff(1:nrat,r)){
          temp[l]<-sum(p_j*w[,data[h,s]])
          l<-l+1
        }
        eprim_tmp[h,r]<-sum(temp)
      }
      eprim[h]<-sum(eprim_tmp[h,])/(nrat*(nrat-1))
    }
    
    P_e<-mean(eprim)
    
    d1<-sum(((1-P_e)*oprim-2*(1-P_o)*eprim)^2)/n
    d2<-(P_o*P_e-2*P_e+P_o)^2
    var<-(d1-d2)/(n*(1-P_e)^4)
    kappa<-(P_o-P_e)/(1-P_e)
  }
  
  if (chance==3){
    #chance=3
    eprim<-rep(NA,n)
    p_j<-colMeans(sweep(freq,1,nrat,"/"))
    eprim_tmp<-matrix(NA,nrow=n,ncol=nrat)
    eprim<-rep(NA,n)
    
    for (h in 1:n){
      for (r in 1:nrat){
        temp<-rep(NA,nrat-1)
        l<-1
        for (s in seq(1,nrat)[-r]){
          temp[l]<-sum(p_j*w[,data[h,s]])
          l<-l+1
        }
        eprim_tmp[h,r]<-sum(temp)
      }
      eprim[h]<-sum(eprim_tmp[h,])/(nrat*(nrat-1))
    }
    
    P_e<-mean(eprim)
    
    d1<-sum(((1-P_e)*oprim-2*(1-P_o)*eprim)^2)/n
    d2<-(P_o*P_e-2*P_e+P_o)^2
    var<-(d1-d2)/(n*(1-P_e)^4)
    kappa<-(P_o-P_e)/(1-P_e)
  }
  
  
  
  if (chance==2){
    #chance=2
    eprim<-rep(NA,n)
    M<-replicate(nrat,matrix(0,ncol=ncat,nrow=1), simplify=FALSE)
    M<-lapply(1:nrat,f1<-function(idd){prop.table(table(factor(data[,idd],levels=names(table(all)))))})
    eprim_tmp<-matrix(NA,nrow=n,ncol=nrat)
    eprim<-rep(NA,n)
    
    for (h in 1:n){
      for (r in 1:nrat){
        temp<-rep(NA,nrat-1)
        l<-1
        for (s in seq(1,nrat)[-r]){
          temp[l]<-sum(M[[r]]*w[,data[h,s]])
          l<-l+1
        }
        eprim_tmp[h,r]<-sum(temp)
      }
      eprim[h]<-sum(eprim_tmp[h,])/(nrat*(nrat-1))
    }
    
    P_e<-mean(eprim)
    
    d1<-sum(((1-P_e)*oprim-2*(1-P_o)*eprim)^2)/n
    d2<-(P_o*P_e-2*P_e+P_o)^2
    var<-(d1-d2)/(n*(1-P_e)^4)
    kappa<-(P_o-P_e)/(1-P_e)
  }
  
  return(c(kappa,sqrt(var)))
}



#' Computation of the proportion of (dis)agreement, MAE and MSE for more than two raters
#'
#' This function computes a weighted (dis)agreement in the presence of two or more observers and its standard error using the delta method. Note that this function requires a dataset with one rating per patient/object and observer.
#'
#' @param data: NxR matrix. dataset with the R ratings for the N participants coded coded from 1 to K, where K is the number of categories of the scale.
#' @param weight: weighting scheme to be used, "binary_disagree" for equal disagreement weights leading to the proportion of disagreement, "binary_agree" for equal disagreement weights leading to the proportion of agreement, "linear" for linear disagreement weights (leading to the mean absolute error) and "quadratic" for quadratic disagreement weights (leading to the mean squared error)
#' @return a vector where the first element is the weighted kappa coefficient and the second its standard error computed by the delta method
#' @author Sophie Vanbelle \email{sophie.vanbelle@@maastrichtuniversity.nl}
#' @details This function computes weighted disagreement and its standard error using the delta method.    
#' @references H. J. A. Schouten. Measuring pairwise agreement among many observers. Biometrical Journal, 22(6):497–504, 1980.
#' @references H. J. A. Schouten. Measuring pairwise interobserver agreement when all subjects are judged by the same observers. Statistica Neerlandica, 36:45–61, 1982
#' @references Vanbelle S. et al. (submitted). Agreement between several observers on a ordinal scale
#' @export
#' @examples
#'
#' #dataset from Ayres de Campos
#'
#' data(CTG)
#' qow(CTG,"binary_agree")
#'



qow<-function(data,weight){
  
  #CREATE THE COMPLETE CASE DATASET
  data_<-na.omit(cbind(data,seq(1,length(data[,1]))))
  data2<-data.matrix(data_[,1:ncol(data)])
  idd2<-seq(1,nrow(data))
  N<-nrow(data)
  
  #CREATE FACTORS WITH THE CLASSIFICATION OF THE OBSERVERS (not useful)
  all<-c(data2)
  all.f<-factor(all,levels=names(table(all)))
  ncat<-nlevels(all.f)        #number of categories of the scale
  
  K<-max(idd2)						    #number of items
  nrater<-ncol(data)					#number of raters
  npair<-nrater*(nrater-1)/2  #number of distinct pairs 
  
  
  score <- matrix(1:ncat, nrow = ncat, ncol = ncat, byrow = TRUE)
  
  #CREATE THE WEIGHTED MATRX
  
  if (weight=="binary_disagree"){
    w<-1-diag(ncat)
  }
  if (weight=="binary_agree"){
    w<-diag(ncat)
  }
  if (weight=="linear"){
    w<-abs(score - t(score))
  }
  if (weight=="quadratic"){
    w<-(score - t(score))^2
  }
  
  
  #CREATE THE OBSERVED DISAGREEMENT
  Po<-rep(NA,npair)
  ps<-matrix(NA,ncat,ncat)
  
  count<-1
  for (t in 1:(nrater-1)){
    for (u in (t+1):nrater){
      ps<-table(factor(data2[,t],levels=names(table(all))),factor(data2[,u],levels=names(table(all))))/sum(table(data2[,t],data2[,u]))
      Po[count]<-sum(w*ps)
      
      count<-count+1
    }
    
  }
  
  PO<-mean(Po)
  
  #MATRIX WITH THE WEIGHTS w(ia,ib)
  Wo<-replicate(npair,matrix(0,ncol=K,nrow=1), simplify=FALSE)
  
  count<-1
  for (t in 1:(nrater-1)){
    for (u in (t+1):nrater){
      for (s in 1:K){ 
        Wo[[count]][s]<-w[data[s,t],data[s,u]]^2
      }
      count<-count+1
    }
  }
  
  wia_ib<-2*Reduce("+",Wo)/(nrater*(nrater-1))
  
  
  var<-(sum(wia_ib)-N*PO^2)/N^2
  
  return(c(PO,sqrt(var)))
  
}


#' Computation of the weighted kappa coefficient for more than two raters when only a summary of the assessments made by the different observers is available per participant
#'
#' This function computes a weighted kappa coefficient under chance definition 3 (all observers are assumed to have the same probability to classify the participants in the K categories) and its standard error using the delta method. Note that this function requires a dataset with the distribution of the ratings in the K categories for each participant, not the individual ratings.
#'
#' @param data: NxK matrix. dataset with the summary of the ratings per participant for the N participants, where K is the number of categories of the scale.
#' @param weight: weighting scheme to be used, "binary" for equal disagreement weights, "linear" for linear disagreement weights and "quadratic" for quadratic disagreement weights
#' @return a vector where the first element is the weighted kappa coefficient and the second its standard error computed by the delta method
#' @author Sophie Vanbelle \email{sophie.vanbelle@@maastrichtuniversity.nl}
#' @details This function computes weighted kappa coefficient and its standard error using the delta method.    
#' @references H. J. A. Schouten. Measuring pairwise agreement among many observers. Biometrical Journal, 22(6):497–504, 1980.
#' @references Vanbelle S. et al. (submitted). Agreement between several observers on a ordinal scale
#' @export
#' @examples
#'
#' #dataset from Ayres de Campos 
#'
#' data(CTG)
#' kappaw_sum(CTG,"quadratic",2)
#'
# R1<-c(0,0,2,0,1,0,1,1,0,0,0,1,1,2,2,2,1,1,0,0,1,0,0,1,0,1,1,1,1,2,2,0,0)
# R2<-c(0,0,2,0,2,0,1,1,0,0,0,1,2,2,2,2,1,1,1,0,2,0,1,1,1,2,2,1,2,2,2,0,1)
# R3<-c(1,1,2,0,1,0,1,1,0,1,0,1,1,2,1,2,1,1,1,0,2,0,0,0,0,1,1,1,2,1,1,0,1)
# 
# data<-cbind(R1,R2,R3)+1
# freq<-t(sapply(1:33,fun2<-function(i){table(factor(data[i,],levels=c(1,2,3)))}))


kappaw_sum<-function(data,weight){
  
n<-nrow(data)       #number of subjects
nrat<-rowSums(data)[1] #number of raters
npair<-nrat*(nrat-1)/2 #number of pairs of raters
ncat<-ncol(data)    		#number of categories of the scale


#CREATE THE WEIGHTED MATRX

score <- matrix(1:ncat, nrow = ncat, ncol = ncat, byrow = TRUE)

if (weight=="binary"){
  w<-diag(ncat)
}
if (weight=="linear"){
  w<-1-abs(score - t(score))/(ncat-1)
}
if (weight=="quadratic"){
  w<-1-(score - t(score))^2/(ncat-1)^2
}


#COMPUTATION OF WEIGHTED DISAGREEMENT

p<-matrix(NA,nrow=ncat,ncol=ncat)
oprim<-rep(NA,n)
for (h in 1:n){
  for (i in 1:ncat){
    for (j in 1:ncat){
      if (i==j){
        p[i,j]<-w[i,j]*data[h,i]*(data[h,i]-1)/(nrat*(nrat-1))
      }
      if (i!=j){
        p[i,j]<-w[i,j]*data[h,i]*data[h,j]/(nrat*(nrat-1))
      }
    }
    
  }
  oprim[h]<-sum(p)
}

P_o<-mean(oprim)


#COMPUTATION OF WEIGHTED DISAGREEMENT

eprim<-rep(NA,n)
p_j<-colMeans(sweep(freq,1,nrat,"/"))
p<-matrix(NA,nrow=ncat,ncol=ncat)
for (h in 1:n){
  for (i in 1:ncat){
    for (j in 1:ncat){
      p[i,j]<-w[i,j]*data[h,i]*p_j[j]/nrat
    }
    
  }
  eprim[h]<-sum(p)
}

P_e<-mean(eprim)

#COMPUTATION OF WEIGHTED KAPPA AND STANDARD ERROR

d1<-sum(((1-P_e)*oprim-2*(1-P_o)*eprim)^2)/n
d2<-(P_o*P_e-2*P_e+P_o)^2
var<-(d1-d2)/(n*(1-P_e)^4)
kappa<-(P_o-P_e)/(1-P_e)

return(c(kappa,sqrt(var)))

}



#' Computation of the weighted (dis)agreement for more than two raters when only a summary of the assessments made by the different observers is available per participant
#'
#' This function computes a weighted (dis)agreement coefficient and its standard error using the delta method. Note that this function requires a dataset with the distribution of the ratings in the K categories for each participant, not the individual ratings.
#'
#' @param data: NxK matrix. dataset with the summary of the ratings per participant for the N participants, where K is the number of categories of the scale.
#' @param weight: weighting scheme to be used, "binary_disagree" for equal disagreement weights leading to the proportion of disagreement, "binary_agree" for equal disagreement weights leading to the proportion of agreement, "linear" for linear disagreement weights (leading to the mean absolute error) and "quadratic" for quadratic disagreement weights (leading to the mean squared error)
#' @return a vector where the first element is the weighted (dis)agreement coefficient and the second its standard error computed by the delta method
#' @author Sophie Vanbelle \email{sophie.vanbelle@@maastrichtuniversity.nl}
#' @details This function computes weighted (dis)agreement coefficient and its standard error using the delta method.    
#' @references H. J. A. Schouten. Measuring pairwise agreement among many observers. Biometrical Journal, 22(6):497–504, 1980.
#' @references Vanbelle S. et al. (submitted). Agreement between several observers on a ordinal scale
#' @export
#' @examples
#'
#' #dataset from Ayres de Campos 
#'
#' data(CTG)
#' kappaw_sum(CTG,"quadratic",2)
#'
# R1<-c(0,0,2,0,1,0,1,1,0,0,0,1,1,2,2,2,1,1,0,0,1,0,0,1,0,1,1,1,1,2,2,0,0)
# R2<-c(0,0,2,0,2,0,1,1,0,0,0,1,2,2,2,2,1,1,1,0,2,0,1,1,1,2,2,1,2,2,2,0,1)
# R3<-c(1,1,2,0,1,0,1,1,0,1,0,1,1,2,1,2,1,1,1,0,2,0,0,0,0,1,1,1,2,1,1,0,1)
# 
# data<-cbind(R1,R2,R3)+1
# freq<-t(sapply(1:33,fun2<-function(i){table(factor(data[i,],levels=c(1,2,3)))}))


qow_sum<-function(data,weight){
  
  n<-nrow(data)          #number of subjects
  nrat<-rowSums(data)[1] #number of raters
  npair<-nrat*(nrat-1)/2 #number of pairs of raters
  ncat<-ncol(data)    	 #number of categories of the scale
  
  
  #CREATE THE WEIGHTED MATRX
  
  score <- matrix(1:ncat, nrow = ncat, ncol = ncat, byrow = TRUE)
  
  if (weight=="binary_agree"){
    w<-diag(ncat)
  }
  
  if (weight=="binary_disagree"){
    w<-1-diag(ncat)
  }
  if (weight=="linear"){
    w<-abs(score - t(score))
  }
  if (weight=="quadratic"){
    w<-(score - t(score))^2
  }
  
  
  #COMPUTATION OF WEIGHTED DISAGREEMENT
  
  p<-matrix(NA,nrow=ncat,ncol=ncat)
  oprim<-rep(NA,n)
  wirit2<-rep(NA,n)
  for (h in 1:n){
    for (i in 1:ncat){
      for (j in 1:ncat){
        if (i==j){
          p[i,j]<-w[i,j]*data[h,i]*(data[h,i]-1)/(nrat*(nrat-1))
        }
        if (i!=j){
          p[i,j]<-w[i,j]*data[h,i]*data[h,j]/(nrat*(nrat-1))
        }
      }
      
    }
    oprim[h]<-sum(p)
  }
  
  for (h in 1:n){
    for (i in 1:ncat){
      for (j in 1:ncat){
        if (i==j){
          p[i,j]<-w[i,j]^2*data[h,i]*(data[h,i]-1)/(nrat*(nrat-1))
        }
        if (i!=j){
          p[i,j]<-w[i,j]^2*data[h,i]*data[h,j]/(nrat*(nrat-1))
        }
      }
      
    }
    wirit2[h]<-sum(p)
  }
  
  P_o<-mean(oprim)
  
  #COMPUTATION OF STANDARD ERROR
  
  d1<-sum(wirit2)
  var<-(d1-n*P_o^2)/(n^2)
 
  
  return(c(P_o,sqrt(var)))
  
}


#' Determination of the length of the confidence interval around a kappa coefficient
#' 
#' This function simulates data to find the distribution of the length of the confidence interval around a weigted kappa for a given number of observers, a given expected agreement level, given marginal distribution, number of participants and given level of uncertainty
#'
#'
#' @param wkappa: value of the expected weighted kappa coefficient
#' @param weight: weighting scheme to be used, "binary" for equal disagreement weights, "linear" for linear disagreement weights and "quadratic" for quadratic disagreement weights
#' @param nrat: number of raters
#' @param npar: number of participants rated by the the raters
#' @param prob: a vector of length K with the proportion of participants expected in each category of the scale, where K is the number of categories on the scale. 
#' @param alpha: level of uncertainty
#' @param nsim: number of simulated datasets
#' @param seed: seed for the generation of the datasets
#' @return a vector of length nsim with the length of the confidence interval for each simulation
#' @author Sophie Vanbelle \email{sophie.vanbelle@@maastrichtuniversity.nl}
#' @details This function computes weighted (dis)agreement coefficient and its standard error using the delta method.    
#' @references H. J. A. Schouten. Measuring pairwise agreement among many observers. Biometrical Journal, 22(6):497–504, 1980.
#' @references Vanbelle S. et al. (submitted). Agreement between several observers on a ordinal scale
#' @importFrom orddata rmvord
#' @export
#' @examples
#'

sizeCI_kappaw<-function(wkappa,weight,nrat,npar,prob,alpha,nsim,seed){
 
  probs1<-rep(list(prob),nrat) #vector with the marginal probability distribution
  Cor1<-matrix(wkappa,ncol=nrat,nrow=nrat) #correlation matrix
  diag(Cor1)<-1
  
  set.seed(seed)
  
  width<-rep(NA,nsim)
  for (i in 1:nsim){
    data.sim<-rmvord(n=npar,probs=probs1,Cor=Cor1,showCor_norm = FALSE)
    a<-kappaw(data.sim,weight,3)
    width[i]<-2*qnorm(1-alpha/2)*a[2]
  }
  return(width)
}

# library(orddata)
# w<-sizeCI_kappaw(0.75,"linear",3,30,c(0.4,0.4,0.2),0.05,100,1234)
# 
# plot(density(w))
# percentage<-ifelse(w<0.20,1,0)
# mean(percentage)
# 
# summary(w)





#' Determination of the length of the confidence interval around a (dis)agreement coefficient
#' 
#' This function simulates data to find the distribution of the length of the confidence interval around a (dis)agreement coefficient for a given number of observers, a given expected agreement level, given marginal distribution, number of participants and given level of uncertainty
#'
#'
#' @param wqo: value of the expected (dis)agreement coefficient
#' @param weight: weighting scheme to be used, "binary_disagree" for equal disagreement weights leading to the proportion of disagreement, "binary_agree" for equal disagreement weights leading to the proportion of agreement, "linear" for linear disagreement weights (leading to the mean absolute error) and "quadratic" for quadratic disagreement weights (leading to the mean squared error)
#' @param nrat: number of raters
#' @param npar: number of participants rated by the the raters
#' @param prob: a vector of length K with the proportion of participants expected in each category of the scale, where K is the number of categories on the scale. 
#' @param alpha: level of uncertainty
#' @param nsim: number of simulated datasets
#' @param seed: seed for the generation of the datasets
#' @return a vector of length nsim with the length of the confidence interval for each simulation
#' @author Sophie Vanbelle \email{sophie.vanbelle@@maastrichtuniversity.nl}
#' @details This function computes weighted (dis)agreement coefficient and its standard error using the delta method.    
#' @references H. J. A. Schouten. Measuring pairwise agreement among many observers. Biometrical Journal, 22(6):497–504, 1980.
#' @references Vanbelle S. et al. (submitted). Agreement between several observers on a ordinal scale
#' @importFrom orddata rmvord
#' @export
#' @examples
#'


sizeCI_qow<-function(wqo,weight,nrat,npar,prob,alpha,nsim,seed){
  
  probs1<-rep(list(prob),nrat)#vector with the marginal probability distribution
  Cor1<-matrix(wqo,ncol=nrat,nrow=nrat) #correlation matrix
  diag(Cor1)<-1
  
  set.seed(seed)
  
  width<-rep(NA,nsim)
  for (i in 1:nsim){
    data.sim<-rmvord(n=npar,probs=probs1,Cor=Cor1,showCor_norm = FALSE)
    a<-qow(data.sim,weight)
    width[i]<-2*qnorm(1-alpha/2)*a[2]
  }
  return(width)
}

# library(orddata)
# w<-sizeCI_qow(0.80,"binary_agree",3,30,c(0.4,0.4,0.2),0.05,100,1234)
# 
# plot(density(w))
# percentage<-ifelse(w<0.20,1,0)
# mean(percentage)
# 
# summary(w)


#' Determination of the empirical power in hypothesis testing about a kappa coefficient
#' 
#' This function simulates data to find the distribution of the empirical power when testing statistical hypotheses about weigthed kappa for a given number of observers, a given expected agreement level, given marginal distribution, number of participants and given type 1 error
#'
#'
#' @param wkappa0: value of the expected weighted kappa coefficient under the null hypothesis
#' @param wkappa1: value of the expected weighted kappa coefficient under the alternative hypothesis
#' @param weight: weighting scheme to be used, "binary" for equal disagreement weights, "linear" for linear disagreement weights and "quadratic" for quadratic disagreement weights
#' @param nrat: number of raters
#' @param npar: number of participants rated by the the raters
#' @param prob: a vector of length K with the proportion of participants expected in each category of the scale, where K is the number of categories on the scale. 
#' @param alpha: type I error
#' @param nsim: number of simulated datasets
#' @param seed: seed for the generation of the datasets
#' @return a vector of length nsim with the length of the confidence interval for each simulation
#' @author Sophie Vanbelle \email{sophie.vanbelle@@maastrichtuniversity.nl}
#' @details This function computes weighted (dis)agreement coefficient and its standard error using the delta method.    
#' @references H. J. A. Schouten. Measuring pairwise agreement among many observers. Biometrical Journal, 22(6):497–504, 1980.
#' @references Vanbelle S. et al. (submitted). Agreement between several observers on a ordinal scale
#' @importFrom orddata rmvord
#' @export
#' @examples
#'

test_kappaw<-function(wkappa0,wkappa1,weight,nrat,npar,prob,alpha,nsim,seed){
  
  probs1<-rep(list(prob),nrat) #vector with the marginal probability distribution
  Cor1<-matrix(wkappa1,ncol=nrat,nrow=nrat) #correlation matrix
  diag(Cor1)<-1
  
  set.seed(seed)
  
  power.kw<-rep(NA,nsim)
  for (i in 1:nsim){
    data.sim<-rmvord(n=npar,probs=probs1,Cor=Cor1,showCor_norm = FALSE)
    se_wkappa1<-kappaw(data.sim,weight,3)[2]
    power.kw[i]<-pnorm((wkappa1-wkappa0)/se_wkappa1-qnorm(1-alpha))
  }
  return(power.kw)
}
# 
# library(orddata)
# result<-test_kappaw(0.70,0.80,"linear",3,30,c(0.4,0.4,0.2),0.05,100,1234)
# 
# plot(density(result))
# percentage<-ifelse(result>=(1-0.2),1,0)
# mean(percentage)
# 
# summary(result)






#' Determination of the empirical power in hypothesis testing about a (dis)agreement coefficient
#' 
#' This function simulates data to find the distribution of the empirical power when testing statistical hypotheses about a (dis)agreement coefficient for a given number of observers, a given expected agreement level, given marginal distribution, number of participants and given type 1 error
#'
#'
#' @param wqo0: value of the expected (dis)agreement coefficient under the null hypothesis
#' @param wqo1: value of the expected (dis)agreement coefficient under the alternative hypothesis
#' @param weight: weighting scheme to be used, "binary" for equal disagreement weights, "linear" for linear disagreement weights and "quadratic" for quadratic disagreement weights
#' @param nrat: number of raters
#' @param npar: number of participants rated by the the raters
#' @param prob: a vector of length K with the proportion of participants expected in each category of the scale, where K is the number of categories on the scale. 
#' @param alpha: type I error
#' @param nsim: number of simulated datasets
#' @param seed: seed for the generation of the datasets
#' @return a vector of length nsim with the length of the confidence interval for each simulation
#' @author Sophie Vanbelle \email{sophie.vanbelle@@maastrichtuniversity.nl}
#' @details This function computes weighted (dis)agreement coefficient and its standard error using the delta method.    
#' @references H. J. A. Schouten. Measuring pairwise agreement among many observers. Biometrical Journal, 22(6):497–504, 1980.
#' @references Vanbelle S. et al. (submitted). Agreement between several observers on a ordinal scale
#' @importFrom orddata rmvord
#' @export
#' @examples
#'

test_qow<-function(wqo0,wqo1,weight,nrat,npar,prob,alpha,nsim,seed){
  
  probs1<-rep(list(prob),nrat) #vector with the marginal probability distribution
  Cor1<-matrix(wqo1,ncol=nrat,nrow=nrat) #correlation matrix
  diag(Cor1)<-1
  
  set.seed(seed)
  
  power.kw<-rep(NA,nsim)
  for (i in 1:nsim){
    data.sim<-rmvord(n=npar,probs=probs1,Cor=Cor1,showCor_norm = FALSE)
    se_wqo1<-qow(data.sim,weight)[2]
    power.kw[i]<-pnorm((wqo1-wqo0)/se_wqo1-qnorm(1-alpha))
  }
  return(power.kw)
}
# 
# 
# result<-test_qow(0.80,0.90,"binary_agree",3,30,c(0.4,0.4,0.2),0.05,100,1234)
# 
# plot(density(result))
# percentage<-ifelse(result>=(1-0.2),1,0)
# mean(percentage)
# 
# summary(result)




correlo1<-function(data){
  
  nrat<-ncol(data)
  
  result<-matrix(NA,ncol=nrat,nrow=nrat)
  
  for (i in 1:nrat){
    for (j in 1:nrat){
      
      data_temp<-as.data.frame(data[,c(i,j)])
      
      result[i,j]<-qow(data_temp,"binary_agree")[1]
    }
  }
  
  colnames(result)<-colnames(data)
  rownames(result)<-colnames(data)
  
  corrplot(result, method="shade",addCoef.col = rgb(0,0,0, alpha = 1),cl.lim=c(0,1),number.cex=1.5,type="upper",diag=F)
  
}


correlo2<-function(data){
  
  nrat<-ncol(data)
  
  result<-matrix(NA,ncol=nrat,nrow=nrat)
  
  for (i in 1:nrat){
    for (j in 1:nrat){
      
      data_temp<-as.data.frame(data[,c(i,j)])
      
      result[i,j]<-qow(data_temp,"binary_disagree")[1]
    }
  }
  
  colnames(result)<-colnames(data)
  rownames(result)<-colnames(data)
  
  corrplot(result, method="shade",addCoef.col = rgb(0,0,0, alpha = 1),cl.lim=c(0,1),number.cex=1.5,type="upper",diag=F)
  
}


correlo3<-function(data){
  
  nrat<-ncol(data)
  
  result<-matrix(NA,ncol=nrat,nrow=nrat)
  
  for (i in 1:nrat){
    for (j in 1:nrat){
      
      data_temp<-as.data.frame(data[,c(i,j)])
      
      result[i,j]<-qow(data_temp,"linear")[1]
    }
  }
  
  colnames(result)<-colnames(data)
  rownames(result)<-colnames(data)
  
  corrplot(result, method="shade",addCoef.col = rgb(0,0,0, alpha = 1),cl.lim=c(0,1),number.cex=1.5,type="upper",diag=F)
  
}

correlo4<-function(data){
  
  nrat<-ncol(data)
  
  result<-matrix(NA,ncol=nrat,nrow=nrat)
  
  for (i in 1:nrat){
    for (j in 1:nrat){
      
      data_temp<-as.data.frame(data[,c(i,j)])
      
      result[i,j]<-qow(data_temp,"quadratic")[1]
    }
  }
  
  colnames(result)<-colnames(data)
  rownames(result)<-colnames(data)
  
  corrplot(result, method="shade",addCoef.col = rgb(0,0,0, alpha = 1),cl.lim=c(0,1),number.cex=1.5,type="upper",diag=F)
  
}


correlo5<-function(data){
  
  nrat<-ncol(data)
  
  result<-matrix(NA,ncol=nrat,nrow=nrat)
  
  for (i in 1:nrat){
    for (j in 1:nrat){
      
      data_temp<-as.data.frame(data[,c(i,j)])
      
      result[i,j]<-kappaw(data_temp,"linear",2)[1]
    }
  }
  
  colnames(result)<-colnames(data)
  rownames(result)<-colnames(data)
  
  corrplot(result, method="shade",addCoef.col = rgb(0,0,0, alpha = 1),cl.lim=c(0,1),number.cex=1.5,type="upper",diag=F)
  
}


correlo6<-function(data){
  
  nrat<-ncol(data)
  
  result<-matrix(NA,ncol=nrat,nrow=nrat)
  
  for (i in 1:nrat){
    for (j in 1:nrat){
      
      data_temp<-as.data.frame(data[,c(i,j)])
      
      result[i,j]<-kappaw(data_temp,"linear",3)[1]
    }
  }
  
  colnames(result)<-colnames(data)
  rownames(result)<-colnames(data)
  
  corrplot(result, method="shade",addCoef.col = rgb(0,0,0, alpha = 1),cl.lim=c(0,1),number.cex=1.5,type="upper",diag=F)
  
}

correlo7<-function(data){
  
  nrat<-ncol(data)
  
  result<-matrix(NA,ncol=nrat,nrow=nrat)
  
  for (i in 1:nrat){
    for (j in 1:nrat){
      
      data_temp<-as.data.frame(data[,c(i,j)])
      
      result[i,j]<-kappaw(data_temp,"quadratic",2)[1]
    }
  }
  
  colnames(result)<-colnames(data)
  rownames(result)<-colnames(data)
  
  corrplot(result, method="shade",addCoef.col = rgb(0,0,0, alpha = 1),cl.lim=c(0,1),number.cex=1.5,type="upper",diag=F)
  
}


correlo8<-function(data){
  
  nrat<-ncol(data)
  
  result<-matrix(NA,ncol=nrat,nrow=nrat)
  
  for (i in 1:nrat){
    for (j in 1:nrat){
      
      data_temp<-as.data.frame(data[,c(i,j)])
      
      result[i,j]<-kappaw(data_temp,"quadratic",3)[1]
    }
  }
  
  colnames(result)<-colnames(data)
  rownames(result)<-colnames(data)
  
  corrplot(result, method="shade",addCoef.col = rgb(0,0,0, alpha = 1),cl.lim=c(0,1),number.cex=1.5,type="upper",diag=F)
  
}

