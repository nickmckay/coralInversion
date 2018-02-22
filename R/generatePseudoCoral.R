pseudoCoral <- function(start.year = 1701, end.year = 2000,start.inst = 1950,SSS.coef = 0.60,SST.coef = -0.23,SSTvar = 2.5,SSSvar = .5,SST.SSS.corr = -.3,instSNR = 3,proxSNR = .5){

#calc year data
years = seq(start.year,end.year)
nYears = length(years)
inst.ind = which(years >= start.inst)
n.inst = length(inst.ind)

#generate "real" instrumental data
realSST = rnorm(nYears,sd = SSTvar)
realSSS = realSST*SST.SSS.corr + (1-SST.SSS.corr)*rnorm(nYears,sd = SSTvar) #make it covary
#scale to appropriate variability
realSSS = scale(realSSS)*SSSvar


#generate calibration data from pseudo data
instSST = realSST[inst.ind]+rnorm(n.inst,sd=SSTvar/instSNR)
instYears = years[inst.ind]

instSSS = realSSS[inst.ind]+rnorm(n.inst,sd=SSSvar/instSNR)
instd18Osw = instSSS * SSS.coef

instCov = cov(instSST,instSSS)

#generate pseudo d18O
proxy.d18O = realSST*SST.coef + realSSS*SSS.coef
#add noise
proxy.d18O = proxy.d18O + rnorm(nYears,sd=sd(proxy.d18O)/proxSNR)

inst.Fsw = Fsw(proxy.d18O[inst.ind],instSSS,instSST )

pseudo = list(years = years, realSST = realSST, realSSS = realSSS, instSSS = instSSS,
              instSST = instSST, instYears = instYears, proxy.d18O = proxy.d18O, inst.Fsw=inst.Fsw)

return(pseudo)

}



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


