plotReconSummary <- function(coral.data,coral.meta,inst.data,recon,pdf.opt = F){

#see if instrumental data extends beyond coral data
  if(max(coral.data$Year) < max(inst.data$Year)){
    inst.data <- dplyr::filter(inst.data,Year<=max(coral.data$Year))
  }


  inst.ind =which(coral.data$Year %in% inst.data$Year)
coral.data$Year[inst.ind]

#make a plot that shows the ensemble reconstruction, the instrumental data, and the correlations...
timeseries.SST <- plotTimeseriesEnsRibbons(coral.data$Year, recon$ensSST) +
  geom_line(data = inst.data,aes(x = Year,y = scale(SST,scale=F)),color = "red")+
  xlab("Year (AD)")+ ylab("SST anomaly (deg C)")+ggtitle(paste0(coral.meta$iso2kid," - Sea surface temperature"))


cor.sst = corEns(coral.data$Year,recon$ensSST,time2 = inst.data$Year, values2 = inst.data$SST,binstep = 1)

corr.dist.SST <- plotCorrEns(cor.sst$cor.df,corStats = cor.sst$corStats)

ensMed.sst <- apply(recon$ensSST,1,median)
scatter.med.SST <- ggplot()+geom_point(aes(x = ensMed.sst[inst.ind], y = inst.data$SST)) +
  labs(x = "Median reconstructed SST (deg C)",y = "Observed SST (deg C)",title = paste0("Ensemble Median: r = ",as.character(signif(cor(ensMed.sst[inst.ind], inst.data$SST),3)))) + theme_bw()



#repeat for SSS

timeseries.SSS <- plotTimeseriesEnsRibbons(coral.data$Year, recon$ensSSS) +
  geom_line(data = inst.data,aes(x = Year,y = scale(SSS,scale=F)),color = "red")+
  xlab("Year (AD)")+ ylab("SSS anomaly (deg C)")+ggtitle("Sea surface salinity")


cor.SSS = corEns(coral.data$Year,recon$ensSSS,time2 = inst.data$Year, values2 = inst.data$SSS,binstep = 1)

corr.dist.SSS <- plotCorrEns(cor.SSS$cor.df,corStats = cor.SSS$corStats)

ensMed.SSS <- apply(recon$ensSSS,1,median)
scatter.med.SSS <- ggplot()+geom_point(aes(x = ensMed.SSS[inst.ind], y = inst.data$SSS)) +
  labs(x = "Median reconstructed SSS (deg C)",y = "Observed SSS (deg C)",title = paste0("Ensemble Median: r = ",as.character(signif(cor(ensMed.SSS[inst.ind], inst.data$SSS),3)))) + theme_bw()

lay <- rbind(c(1,1,2,3),
             c(4,4,5,6))
summPlot <- gridExtra::grid.arrange(grobs = list(timeseries.SST, scatter.med.SST, corr.dist.SST,timeseries.SSS, scatter.med.SSS, corr.dist.SSS),layout_matrix = lay)

if(pdf.opt){
  ggsave(plot = summPlot,filename = paste0(coral.meta$iso2kid,".Reconstruction.SummaryPlot.pdf"),width = 16, height = 8)
}

return(summPlot)

}
