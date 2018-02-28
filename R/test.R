if(FALSE){
test <- pseudoCoral(SSTvar = 2, instSNR = 100,proxSNR = 100,SST.SSS.corr = -.1)

ens.test = runModel(inst.Fsw = Fsw(d18O.coral = test$proxy.d18O,
                                   SSS = test$instSSS,
                                   SST = test$instSST),
                    d18O = test$proxy.d18O,
                    instCov = cov(test$instSSS,test$instSST),
                    nIt = 10000,
                    nEns = 1,
                    instSST = test$instSST,instYears = test$instYears,
                    years = test$years,instSSS = test$instSSS)

}

