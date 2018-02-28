if(FALSE){

library(tidyverse)

#read in data
direc = "/Users/npm4/Google Drive/coral salinity project"

coral.data = read_csv(file.path(direc,"coralData.csv"),col_names = c("Year","d18O"))
coral.meta <- read_csv(file.path(direc,"coralMetadata.csv"),col_names = TRUE)
inst.data <- read_csv(file.path(direc,"instData.csv"),col_names = c("Year","SST","SSS"))

#remove NAs
coral.data <- na.omit(coral.data)
inst.data <- na.omit(inst.data)

#
inst.ind2 =which(inst.data$Year>=1950)
cov(inst.data$SSS,inst.data$SST)
#run model
recon = runModel(inst.Fsw = coral.meta$Fsw,
                 d18O = coral.data$d18O,
                 instCov = cov(inst.data$SSS[inst.ind2],inst.data$SST[inst.ind2]),
                 nIt = 1000,
                 nEns = 1,
                 instSST = inst.data$SST,
                 instYears = inst.data$Year,
                 years = coral.data$Year,
                 instSSS = inst.data$SSS,
                 SST.coef = -.22,SSS.coef = coral.meta$a2)



itSeq = seq(100000,nrow(recon$allLike),by=1000)
plot(itSeq,log(recon$allLike[itSeq,1]),type = "l")

#plot some results...
###SST
plot(coral.data$Year,scale(recon$ensSST[,1],center = T,scale = F),type = "l")
lines(inst.data$Year,scale(inst.data$SST,center = T, scale = F),col = "red")

plot(recon$ensSST[inst.ind,1],inst.data$SST[inst.ind2] )
cor(recon$ensSST[inst.ind,1],inst.data$SST[inst.ind2] )


###SSS
plot(coral.data$Year,scale(recon$ensSSS[,1],center = T,scale = F),type = "l")
lines(inst.data$Year,scale(inst.data$SSS,center = T, scale = F),col = "red")


plot(recon$ensSSS[inst.ind,1],inst.data$SSS[inst.ind2] )
cor(recon$ensSSS[inst.ind,1],inst.data$SSS[inst.ind2] )



#plot(recon$ensSST[,1]*-.22+recon$ensSSS* 0.97002*.27, coral.data$d18O)

plot(coral.data$Year,scale(-coral.data$d18O,center = T,scale = T),type = "l")
lines(inst.data$Year,scale(inst.data$SST,center = T, scale = T),col = "red")


inst.ind = which(coral.data$Year %in% inst.data$Year[inst.ind2])


plot(recon$ensSSS[inst.ind,1],inst.data$SSS[inst.ind2] )
cor(recon$ensSSS[inst.ind,1],inst.data$SSS[inst.ind2] )


plot(coral.data$d18O[inst.ind],inst.data$SST[inst.ind2] )
cor(coral.data$d18O[inst.ind],inst.data$SST[inst.ind2]  )

plot(coral.data$d18O[inst.ind],inst.data$SSS[inst.ind2]  )
cor(coral.data$d18O[inst.ind],inst.data$SSS[inst.ind2]  )



#with reconstructed data
plot(coral.data$d18O,recon$ensSST)
cor(coral.data$d18O,recon$ensSST)

plot(coral.data$d18O,recon$ensSSS )
cor(coral.data$d18O,recon$ensSSS)

Fsw(d18O.coral = coral.data$d18O,SSS = recon$ensSSS[,1], SST = recon$ensSST[,1],SST.coef = -.22, SSS.coef = coral.meta$a2)
coral.meta$Fsw

}
