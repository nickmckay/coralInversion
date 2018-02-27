pseudoCoral <- function(start.year = 1701, end.year = 2000,start.inst = 1950,SSS.coef = 0.60,SST.coef = -0.23,SSTvar = 1.5,SSSvar = .5,SST.SSS.corr = -.3,instSNR = 3,proxSNR = .5){

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

