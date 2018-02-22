#
#
# #some ensemble ploting code
# for(e in 1:nEns){
#
#   if(e==1){
#     plot(years,ensSST[,e],type="l",col="gray")
#   }
#   lines(years,ensSST[,e],type="l",col="gray")
#
# }
#
# medSST= apply(ensSST,1,median)
# lines(years,medSST)
# lines(years,realSST,col = "red")
#
#
# plot(medSST,realSST)
# cor(medSST,realSST)
#
# #again for SSS
# for(e in 1:nEns){
#
#   if(e==1){
#     plot(years,ensSSS[,e],type="l",col="gray")
#   }
#   lines(years,ensSSS[,e],type="l",col="gray")
#
# }
#
# medSSS= apply(ensSSS,1,median)
# lines(years,medSSS)
# lines(years,realSSS,col = "red")
#
#
#
# plot(medSSS,realSSS)
# cor(medSSS,realSSS)
#
#
#
# #test
# d18Otest = SST*SST.coef + SSS*SSS.coef
# plot(d18Otest,reald18O)
