
if(FALSE){
  nEns = ncol(ens.test$ensSST)

#some ensemble ploting code
for(e in 1:nEns){

  if(e==1){
    plot(test$years,ens.test$ensSST[,e],type="l",col="gray")
  }
  lines(test$years,ens.test$ensSST[,e],type="l",col="gray")

}

medSST= apply(ens.test$ensSST,1,median)
lines(test$years,medSST)
lines(test$years,test$realSST,col = "red")


#two

plot(medSST,test$realSST)
cor(medSST,test$realSST)

medSSS= apply(ens.test$ensSSS,1,median)
plot(medSSS,test$realSSS)
cor(medSSS,test$realSSS)


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
}
