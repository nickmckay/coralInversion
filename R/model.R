runModel <- function(d18O,inst.Fsw, SST.coef = -0.23,
                     SSS.coef = 0.6, sampleStep = 1,nIt = 200,
                     nEns = 10, years, instCov,instYears, instSST, instSSS){

  #remove the mean from all datasets
  d18O = as.numeric(scale(d18O, center = T,scale = F))
  instSST = as.numeric(scale(instSST, center = T,scale = F))
  instSSS = as.numeric(scale(instSSS, center = T,scale = F))


  nYears <- length(years)

  #initialize (based on Fsw)
  SSSvar = (inst.Fsw*var(d18O) - 2*SST.coef*instCov)/SSS.coef^2 #from Russon 2013 Fsw equation



  #instrumental variances
  # inst.SSTvar = SSTvar
  # inst.SSSvar = var(instSSS)
  inst.ind <- which(years %in% instYears)

  #start Gibbs sampler.
  #nIt = 200
  #sampleStep = 1
  #thresh = 1e-8 #.Machine$double.eps
  allLike = matrix(NA,nrow = nIt*nYears,ncol = nEns)
  allFsw = allLike
  allCov = allLike
  #nEns = 5
  ensSST = matrix(NA,ncol = nEns, nrow = nYears)
  ensSSS = ensSST
  t = 0
  for(e in 1:nEns){#loop through ensemble members
    # initSST = rnorm(nYears,sd = SSTvar)
    # initSSS = (d18O - SST.coef*initSST) / SSS.coef
    initSSS = rnorm(nYears,sd = sqrt(abs(SSSvar)))
    initSST = (d18O - SSS.coef*initSSS) / SST.coef



    SST = initSST
    SSS = initSSS
    initLikelihood = Fsw.likelihood(Fsw(d18O.coral = d18O,SSS = initSSS,SST = initSST, SST.coef = SST.coef, SSS.coef = SSS.coef ),obsFsw = inst.Fsw) *
      #sst.likelihood(initSST[inst.ind],instSST) *
      cov.likelihood(cov(SST,SSS),instCov)
    #  var.likelihood(var(initSST[inst.ind]),inst.SSTvar) *
    #  var.likelihood(var(initSSS[inst.ind]),inst.SSSvar)

    like = initLikelihood


    print(paste("e = ",as.character(e)))


    ii = 0
    for(i in 1:nIt){ #run through nIt times times.
      #while(keepRunning){
      if(i%%(nIt/5)==0){print(i)}

      #loop through all years
      for(y in 1:length(SST)){#update sample for this year...
        ii = ii+1

        #propose an update to a year
        proposedSST = SST
        proposedSSS = SSS

        proposedSST[y] = SST[y]+rnorm(1,sd = sampleStep)
        proposedSSS[y] = (d18O[y] - SST.coef*proposedSST[y]) / SSS.coef

        allFsw[i,e] <- Fsw(d18O.coral = d18O,SSS = proposedSSS,SST = proposedSST, SST.coef = SST.coef, SSS.coef = SSS.coef)
        allCov[i,e] <- cov(proposedSST,proposedSSS)



        #test the likelihood
        newLike = Fsw.likelihood(allFsw[i,e],obsFsw = inst.Fsw) *
          #sst.likelihood(proposedSST[inst.ind],instSST)*
          cov.likelihood(allCov[i,e],instCov)
        #  var.likelihood(var(proposedSST[inst.ind]),inst.SSTvar) *
        # var.likelihood(var(proposedSSS[inst.ind]),inst.SSSvar)


        #accept the proposed change? #modified metropolis scheme here.
        if(TRUE){
          accept <- newLike/like > 1
        }else{
          accept <- newLike/like > runif(1) #metrop
        }


        if(accept){
          SST = proposedSST
          SSS = proposedSSS
          like = newLike
        }




        allLike[(i-1)*length(SST)+y,e] = like

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
  return(list(ensSST = ensSST, ensSSS = ensSSS,allLike = allLike, allFsw = allFsw, allCov = allCov, initSST = initSST, initSSS = initSSS))
}

