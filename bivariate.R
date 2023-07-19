---
  title: "Bivariate analysis"
author: "Mehwish Zaman"
date: "`r Sys.Date()`"
output: html_document
---

rm(list = ls())
library(tidydr)
library(dplyr)
library(extRemes)
library(evd)
library(ismev)
library(IDF)




partrue<- c(3.2, 3.5,0.47,0.55, 0.45)## True parameter value
d<- c(1,2)
s1<- partrue[2]/(d[1]^(partrue[4])) ### s(d)= s_o * d^(-eta)
s2<- partrue[2]/(d[2]^(partrue[4]))
u1 <- partrue[1]*s1                    ### u(d)= tilda(u)*s(d)
u2 <- partrue[1]*s2
k1 <- partrue[3]
k2<- partrue[3]
mr1 <- c(u1,s1,k1)
mr2<- c(u2,s2,k2)
data_new<-rbvevd(500, dep = partrue[5], model = "log", mar1 = mr1, mar2 = mr2)
colnames(data_new)<- c("d1","d2")
#data_new
fbvevd(data_new, model="log")
sum(dbvevd(data_new, dep=partrue[5], model="log",mar1 = mr1, mar2 = mr2, log = TRUE))
plot(data_new)



######likelihood for bivariate density
x1<- data_new[,1]
x2<- data_new[,2]
data1 <- data.frame(x1,x2)

dep.gev.fit <- function(par, data1, debug=FALSE)
{
  epsi<- 10^(-6)
  mu<-par[1]
  sigma_o<- par[2]
  k<- par[3]
  eta<- par[4]
  a <- par[5]
  if (debug) cat(mu,sigma_o,k,eta,a," ")
  d <- c(1,2)
  s1 <- (sigma_o)/(d[1])^eta
  s2 <- (sigma_o)/(d[2])^eta
  m1 <- mu * s1
  m2 <- mu * s2
  y1<- (pmax(0, (1 + (k * (data1$x1 - m1)/s1))))^(1/k)
  y2<- ( pmax(0, (1 + (k * (data1$x2 - m2)/s2))))^(1/k)
  y <- (1/(y1^(1/a))+1/(y2^(1/a))) 
  if(any(y1<= 0)||any(y2 <= 0 ) || any(s1 <= 0) || any(s2 <= 0)|| any(y <= 0)) return(10^6)
  if(eta <= epsi || eta >= 1-epsi)return(10^6)
  if(any(a <= epsi ) || any(a > 1-epsi)) return(10^6)
  g<- -sum(-log(s1)-log(s2)-(y)^a + (a - 2)*log(y) + log(y^a - 1+ 1/a)- (1/a)*log(y1)-k*log(y1)- (1/a)*log(y2)-k*log(y2), na.rm = TRUE)
  if (debug) cat(g,"\n")
  return(g)
}
dep.gev.fit(par=partrue, data1=data1, debug=TRUE)
optim2 <- optim(par=partrue, dep.gev.fit, data1=data1, debug=TRUE, method = "Nelder-Mead", hessian = T, control=list(maxit=5000))
optim2


#######independent d-GEV model
x1<-data_new[,1]
x2<-data_new[,2]
annualmax <- data.frame(x1,x2)
likehood<-function(par,d,annualmax) {
  ### IDF-parameters :IDF model for d[1,72]:
  mu<-par[1]
  sigma<-par[2]
  gamma<-par[3]
  eta<-par[4]
  d<-c(1,2)
  epsi<-10^(-16)
  lh.mat<-array(0,dim=dim(annualmax))
  Y<-array(0,dim=dim(annualmax))
  ### Loop over rainfall durations #####
  for (i in 1:length(d)) {
    sigma_d<-sigma*(d[i])^(-eta) 
    mu_d<- mu*sigma_d
    Y[,i]<-1+gamma*(annualmax[,i]-mu_d)/sigma_d
    if(any(eta <= epsi || eta > 1-epsi)) return(10^6)
    if (min(Y[,i],na.rm=TRUE)>0) {
      lh.mat[,i]<- -log(sigma_d)-Y[,i]^(-1/gamma)-(1+1/gamma)*log(Y[,i])   ### GEV-loglikelihood for duration d[i]
    } 
    else {
      lh.mat[,i]<--10^10
    }
  }
  ##################################
  ### Make sum of the log-likelihoods over the years and rainfall durations d, and return
  
  return(-sum(lh.mat,na.rm=TRUE))
}
par0<-c(3.20, 3.50 , 0.47 , 0.55 ) ## Inital Value for Optimisation 
likehood(par = par0, d=d, annualmax=annualmax)
result_r<-optim(par0,likehood,annualmax=annualmax, method = "Nelder-Mead", hessian = T, control=list(maxit=5000)) 
result_r


