##### title: duration dependent dgev
####    author: "Mehwish Zaman"

#INDEX Station_id from_date to_date station height geolatitude   geolongitude   station name     federal state
# 1.    01262     19950901    20230605       446     48.3477     11.8134        München-Flughafen     Bayern 
# 2.     03023       19950901 20161231         22     52.5181    7.3081           Lingen                  Niedersachsen  
# 3.     00662 19971022 20230815             81     52.2915   10.4464 Braunschweig                             Niedersachsen      
# 4.      00003 19950901 20110401            202     50.7827    6.0941 Aachen                                   Nordrhein-Westfalen  
#5.      00164 19950901 20230815             54     53.0316   13.9908 Angermünde                               Brandenburg  
# 6.    00183 19950901 20230815             42     54.6791   13.4344 Arkona                                   Mecklenburg-Vorpommern           
# 7.    00198 19950901 20230815            164     51.3744   11.2920 Artern                                   Thüringen        
# 8.    00222 19950901 20230815            393     50.5908   12.7139 Aue                                      Sachsen
# 9.    00282 19950901 20230815            240     49.8743   10.9206 Bamberg                                  Bayern          
# 10. 00298 19950901 20230815              3     54.3405   12.7108 Barth                                    Mecklenburg-Vorpommern  
# 11. 00303 19950901 20230815             55     52.0613   13.4997 Baruth                                   Brandenburg                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
# 12. 00427 19950901 20230815             46     52.3807   13.5306 Berlin Brandenburg                       Brandenburg   
# 13. 01262 19950901 20230815            446     48.3477   11.8134 München-Flughafen                        Bayern   
# 14. 03023 19950901 20161231             22     52.5181    7.3081 Lingen                                   Niedersachsen      
#15.  03552 19950901 20190408             38     52.9037   12.8071 Neuruppin                                Brandenburg           

load(file = "allstations.Rdata")
load(file = "allsthourly.Rdata")
load(file = "AM_all_Stations.Rdata")
library(sf)
library(tmap)
library(dplyr)
library(OpenStreetMap)

########## conditional duration dependence likelihood function for gev model
cdd.gev<- function(x,pr,d){
  epsi<- 10^(-6)
  mu<-pr[1]
  sigma<- pr[2]
  xi<- pr[3]
  eta<- pr[4]
  dep <- pr[5]
  epsi<- 10^(-6)
  if(eta <= epsi || eta >= 1-epsi)return(10^6)
  if(dep <= epsi || dep >= 1-epsi)return(10^6)
  ##### calculates the negative log-likelihood of the first 
  ###column of data(x, matrix of data with more than two columns) using the generalized extreme value distribution with parameters mu, sigma, and k, here d=1.
  nllik<- -sum(evd::dgev(x[,1], loc = mu, scale = sigma, shape = xi, log = TRUE)) 
  for (i in 2:(ncol(x))) {###loop over the durations form d=2 to M=ncol(x)
    s1 <-sigma*(d[i])^(-eta)
    s2 <-sigma*(d[i-1])^(-eta) 
    m1 <- mu*s1
    m2 <- mu*s2
    nllik<- nllik-sum(evd::dbvevd(x[,c(i-1,i)],dep=dep, mar2 = c(m1, s1,xi), mar1 = c(m2,s2,xi), log = TRUE )-evd::dgev(x[,i-1], loc = m2, scale = s2,shape = xi, log=TRUE))
    
  }
  return(nllik)
  
}
d<- 1:72
ff<- st1_maxh[,-1] ### data with multiple columns 
######## for station 1
partrue<- c(1.50, 0.61, 0.82, 0.01, 0.045)
cdd.gev(x=ff,d=d, partrue)
stop1 <- optim( x=ff,fn = cdd.gev,partrue, d=d, method = "Nelder-Mead", hessian = T, control=list(maxit=5000))

############## initial values for station 2
par2<- c(1.75, 0.36, 1.06, 0.04, 0.45)
stop2 <- optim( x=st2max[,-1],fn = cdd.gev,par2, d=d, method = "Nelder-Mead", hessian = T, control=list(maxit=5000))

############## initial values for station 3
par3<- c(1.67, 0.58, 0.68, 0.02, 0.03)
############## initial values for station 4
par4<- c(2.6,0.46, 0.8,0.06, 0.05)
############## initial values for station 5
par5<- c(1.8,0.6, 0.8, 0.015, 0.03)

