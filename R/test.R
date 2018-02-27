if(FALSE){
test = pseudoCoral()

ens.test = runModel(inst.Fsw = Fsw(d18O.coral = test$proxy.d18O,
                                   SSS = test$instSSS,
                                   SST = test$instSST),
                    d18O = test$proxy.d18O,
                    instCov = cov(test$instSSS,test$instSST),
                    nIt = 200,
                    instSST = test$instSST,instYears = test$instYears,
                    years = test$years,instSSS = test$instSSS)

}

