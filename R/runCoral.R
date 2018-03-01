if(FALSE){

library(tidyverse)

#read in data
direc = "/Users/npm4/Google Drive/coral salinity project"

toread = list.dirs(direc)[!(list.dirs(direc) %in% direc)]

for(d in 1:length(toread)){

coral.data = read_csv(file.path(toread[d],"coralData.csv"),col_names = c("Year","d18O"))
coral.meta <- read_csv(file.path(toread[d],"coralMetadata.csv"),col_names = TRUE)
inst.data <- read_csv(file.path(toread[d],"instData.csv"),col_names = c("Year","SST","SSS"))

#remove NAs
coral.data <- na.omit(coral.data)
inst.data <- na.omit(inst.data)

#
inst.ind2 =which(inst.data$Year>=1950)
cov(inst.data$SSS,inst.data$SST)
#run model
recon = runModelSingleCoral(inst.Fsw = coral.meta$Fsw,
                 d18O = coral.data$d18O,
                 instCov = cov(inst.data$SSS[inst.ind2],inst.data$SST[inst.ind2]),
                 nIt = 1000,
                 nEns = 100,
                 instSST = inst.data$SST,
                 instYears = inst.data$Year,
                 years = coral.data$Year,
                 instSSS = inst.data$SSS,
                 SST.coef = -.22,SSS.coef = coral.meta$a2)


summPlot = plotReconSummary(coral.data,coral.meta, inst.data, recon,pdf.opt = F)
ggsave(plot = summPlot,filename = file.path(direc,paste0(coral.meta$iso2kid,".Reconstruction.SummaryPlot.pdf")),width = 16, height = 8)
}

}
