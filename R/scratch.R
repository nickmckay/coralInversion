
#set parameters
start.year = 1701
end.year = 2000
start.inst = 1950


#specify coefficients
SSS.coef = 0.60
SST.coef = -0.23


#calc year data 
years = seq(start.year,end.year)
nYears = length(years)
inst.ind = which(years >= start.inst)
n.inst = length(inst.ind)

#generate pseudoData
realSST = rnorm(nYears,sd = 2.5)
realSSS = (rnorm(nYears,sd = 2)-realSST/15)/3 #make it covary


cov(realSST,realSSS)
cor(realSST,realSSS)


#generate calibration data from pseudo data
instSST = realSST[inst.ind]+rnorm(n.inst,sd=.5)
instYears = years[inst.ind]

instSSS = realSSS[inst.ind]+rnorm(n.inst,sd=.3)
instd18Osw = instSSS * SSS.coef

instCov = cov(instSST,instSSS)
#check it. 
plot(years, realSST,type="l")
lines(instYears,instSST,col="red")

#generate pseudo d18O
reald18O = realSST*SST.coef + realSSS*SSS.coef + rnorm(nYears,sd=0.1)

#calculate instrumental Fsw
Fsw = function(d18O.coral, SSS, SST, SST.coef = -0.23, SSS.coef = 0.60){
  d18O.sw = SSS*SSS.coef
  return(  (var(d18O.sw) + 2*SST.coef*cov(SST,d18O.sw)) /
             var(d18O.coral) )
}


inst.Fsw = Fsw(reald18O[inst.ind],instSSS,instSST )


#set up likelihood functions
Fsw.likelihood = function(simFsw,obsFsw,uncFsw = 0.1){
  #likelihood function for Fsw
  #assume normal distribution of Fsw for now
  like = dnorm(simFsw,mean = obsFsw, sd = uncFsw)
  return(like)
}

sst.likelihood = function(simSST,obsSST,uncSST = 2){
scaledDiffs = (simSST-mean(simSST)) - (obsSST - mean(obsSST))
RMSE = sqrt(t(scaledDiffs)%*%scaledDiffs/length(simSST))
#assume normal disribution for SSE for now
like = dnorm(RMSE,sd = uncSST)
  return(like)
}

#covariance likelihood
cov.likelihood = function(simCov,obsCov,unc = 1){
   #assume normal disribution for SSE for now
  like = dnorm(simCov,mean = obsCov, sd = unc)
  return(like)
}

#add sss?

#initialize (based on Fsw)
SSTvar = 2.5
SSSvar = abs(SSTvar*SST.coef/SSS.coef*inst.Fsw) #this may need to be improved.

#start Gibbs sampler.
nIt = 200
sampleStep = 1
thresh = 1e-8 #.Machine$double.eps
allLike = matrix(NA,nrow = nIt*length(years))
nEns = 5
ensSST = matrix(NA,ncol = nEns, nrow = length(years))
ensSSS = ensSST
t = 0
for(e in 1:nEns){#loop through ensemble members
  # initSST = rnorm(nYears,sd = SSTvar)
  # initSSS = (reald18O - SST.coef*initSST) / SSS.coef
  initSSS = rnorm(nYears,sd = SSSvar)
 initSST = (reald18O - SSS.coef*initSSS) / SST.coef
  
  
  SST = initSST
  SSS = initSSS
  
    initLikelihood = Fsw.likelihood(Fsw(d18O.coral = reald18O,SSS = initSSS,SST = initSST),obsFsw = inst.Fsw) *
    sst.likelihood(initSST[inst.ind],instSST) *
    cov.likelihood(cov(SST,SSS),instCov) 
  

  like = initLikelihood
  
  
  print(paste("e = ",as.character(e)))
  
  
  keepRunning = TRUE
  i = 0
  for(i in 1:nIt){ #run through nIt times times.
  #while(keepRunning){
    if(i%%(nIt/5)==0){print(i)}
    
    #loop through all years
    for(y in 1:length(SST)){#update sample for this year...
      # i = i+1 
      
      #propose an update to a year
      proposedSST = SST
      proposedSSS = SSS
      
      proposedSST[y] = SST[y]+rnorm(1,sd = sampleStep)
      proposedSSS[y] = (reald18O[y] - SST.coef*proposedSST[y]) / SSS.coef
      
      #test the likelihood 
      newLike = Fsw.likelihood(Fsw(d18O.coral = reald18O,SSS = proposedSSS,SST = proposedSST),obsFsw = inst.Fsw) *
        sst.likelihood(proposedSST[inst.ind],instSST)*
        cov.likelihood(cov(proposedSST,proposedSSS),instCov) 
      
      
      #accept the proposed change? #modified metropolis scheme here.

      if(newLike/like > 1){
        SST = proposedSST
        SSS = proposedSSS
        like = newLike
      }
      allLike[(i-1)*length(SST)+y] = like
      
      # if(i > 50000){
      #   if(mean(abs(diff(allLike[i-20000:i]))) < thresh){
      #     keepRunning = FALSE
      #   }
      # }
      # 
      # if(i > 1e6){
      #   keepRunning = FALSE
      # }
      
    }
  }
  
  ensSST[,e] = SST
  ensSSS[,e] = SSS
  
}


#some ensemble ploting code
for(e in 1:nEns){

  if(e==1){
    plot(years,ensSST[,e],type="l",col="gray")
  }
  lines(years,ensSST[,e],type="l",col="gray")
  
}

medSST= apply(ensSST,1,median)
lines(years,medSST)
lines(years,realSST,col = "red")


plot(medSST,realSST)
cor(medSST,realSST)

#again for SSS
for(e in 1:nEns){
  
  if(e==1){
    plot(years,ensSSS[,e],type="l",col="gray")
  }
  lines(years,ensSSS[,e],type="l",col="gray")
  
}

medSSS= apply(ensSSS,1,median)
lines(years,medSSS)
lines(years,realSSS,col = "red")



plot(medSSS,realSSS)
cor(medSSS,realSSS)



#test
d18Otest = SST*SST.coef + SSS*SSS.coef
plot(d18Otest,reald18O)


