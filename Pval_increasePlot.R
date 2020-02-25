options(scipen=999)
library(data.table)

# MATRIX SCAN RSAT CEBPB+ATF4
ms<-read.table(pipe('grep "site" matrix-scan_CEBPB+ATF4.ft|cut -f3,4,5,6,9'),sep="\t",stringsAsFactors=FALSE)
pos <- sort(unique(ms[,3]))

# pval
cebpb<-ms[ms[,1]=="CEBPB",c(3,5)]
cebpb[,2]<-(-log10(cebpb[,2]))
atf4<-ms[ms[,1]=="MM1_HSA_K562_CEBPB_+_MM1_HSA_K562_ATF4",c(3,5)]
atf4[,2]<-(-log10(atf4[,2]))
jund<-ms[ms[,1]=="MM1_HSA_K562_CEBPB_+_MM1_HSA_K562_CEBPD",c(3,5)]
jund[,2]<-(-log10(jund[,2]))

cebpb.t<-aggregate(cebpb[,2], list(cebpb$V3), sum)
atf4.t<-aggregate(atf4[,2], list(atf4$V3), sum)
jund.t<-aggregate(jund[,2], list(jund$V3), sum)

pdf("CEBPB_and_ATF4_matrixScan.pdf")
plot(as.numeric(names(pos.t)),pos.t,col="white",main="Peaks overlapping between CEBPB and ATF4 ",
      ylab = "# Predicted binding sites",xlab="Peak center",
     yaxt="n",xaxt="n",xlim=c(1,201), ylim=c(0,700) )
axis(1, labels=seq(-100, 100, by = 50) ,at = seq(1, 201, by = 50), las=1)
axis(2, labels=seq(0, 200, by = 25) ,at = seq(0, 200, by = 25), las=1)
lines(smooth.spline(as.numeric(cebpb.t[,1]),(cebpb.t[,2]),spar=.3),lwd=3)
lines(smooth.spline(as.numeric(atf4.t[,1]),(atf4.t[,2]),spar=.3),lwd=3,col="blue")
lines(smooth.spline(as.numeric(jund.t[,1]),(jund.t[,2]),spar=.3),lwd=3,col="red")
legend("topright",legend=c("CEBPB matrix","CEBPB+ATF4 matrix","CEBPB+CEBPD matrix"),fill=c("black","blue","red"),bty = "n")
dev.off()

#############################################################################################################################
ms<-read.table(pipe('grep "site" MS_cebpb_cebpd.ft'),sep="\t",stringsAsFactors=FALSE)
ms<-ms[ms[,9]<0.01,]
table(ms[,3])
ms<-ms[ms[,5]>80 & ms[,6]<120,]
table(ms[,3])

cebpb<-read.table(pipe('grep "site" MS_cebpb_cebpd.ft|cut -f3,4,5,6,9'),sep="\t",stringsAsFactors=FALSE)
cebpb<-cebpb[cebpb[,1]=="CEBPB",]
cebpb<-cebpb[cebpb[,5]<0.0001,]

CEBPD<-read.table(pipe('grep "site" MS_cebpb_cebpd.ft|cut -f3,4,5,6,9'),sep="\t",stringsAsFactors=FALSE)
CEBPD<-CEBPD[CEBPD[,1]=="MM1_HSA_K562_CEBPB_+_MM1_HSA_K562_CEBPD",]
CEBPD<-CEBPD[CEBPD[,5]<0.0001,]



plot(density(rowMeans(cbind(cebpb[,3],cebpb[,4]))))
lines(density(rowMeans(cbind(CEBPD[,3],CEBPD[,4]))),col="red")


pdf("CEBPB+CEBPD_density_plot_pval<0.0001.pdf")
plot(density(rowMeans(cbind(cebpb[,3],cebpb[,4]))),col="black",xaxt="n",xlab="CEBPB+CEBPD peakCenter",main="pval<0.0001",ylab="Predicted BindingSites")
axis(1, labels=seq(-100, 100, by = 50) ,at = seq(1, 201, by = 50), las=1)
lines(density(rowMeans(cbind(CEBPD[,3],CEBPD[,4]))),col="red")
legend("topright",legend=c("CEBPB matrix","CEBPB+CEBPD matrix"),fill=c("black","red"),bty = "n")
dev.off()

pdf("CEBPB+CEBPD_density_plot_pval<0.1.pdf")
plot(density(rowMeans(cbind(CEBPD[,3],CEBPD[,4]))),col="red",xaxt="n",xlab="CEBPB+CEBPD peakCenter",main="pval<0.1",ylab="Predicted BindingSites")
axis(1, labels=seq(-100, 100, by = 50) ,at = seq(1, 201, by = 50), las=1)
lines(density(rowMeans(cbind(cebpb[,3],cebpb[,4]))),col="black")
legend("topright",legend=c("CEBPB matrix","CEBPB+CEBPD matrix"),fill=c("black","red"),bty = "n")
dev.off()
#############################################################################################################################
# trial without motif only

#getting the fastq and ran the matrixScan
cebpb_cebpd <- commonPeaks(target_peak_id="MM1_HSA_K562_CEBPB",
                                       motif_only_for_target_peak=FALSE,
                                       compared_peak_id="MM1_HSA_K562_CEBPD",
                                       motif_only_for_compared_peak=FALSE,
                                       methylation_profile_in_narrow_region=FALSE)

cebpb_cebpd_peaks <-cebpb_cebpd[,1][[1]]@common_peak # 517
cebpb_cebpd_peaks[,2] <- cebpb_cebpd_peaks[,2]-100
cebpb_cebpd_peaks[,3] <- cebpb_cebpd_peaks[,3]+100
write.table(cebpb_cebpd_peaks,"cebpb_cebpd_peaks_withoutMotif.bed",sep="\t",quote=F,row.names=F,col.names=F)




ms<-read.table(pipe('grep "site" matrix-scan_2020-01-21.191945_9qszqY.ft.1'),sep="\t",stringsAsFactors=FALSE)
ms<-ms[ms[,9]<0.01,]
table(ms[,3])
ms<-ms[ms[,5]>80 & ms[,6]<120,]
table(ms[,3])

cebpb<-read.table(pipe('grep "site" matrix-scan_2020-01-21.191945_9qszqY.ft.1|cut -f3,4,5,6,9'),sep="\t",stringsAsFactors=FALSE)
cebpb<-cebpb[cebpb[,1]=="CEBPB",]
cebpb<-cebpb[cebpb[,5]<0.1,]

CEBPD<-read.table(pipe('grep "site" matrix-scan_2020-01-21.191945_9qszqY.ft.1|cut -f3,4,5,6,9'),sep="\t",stringsAsFactors=FALSE)
CEBPD<-CEBPD[CEBPD[,1]=="MM1_HSA_K562_CEBPB_+_MM1_HSA_K562_CEBPD",]
CEBPD<-CEBPD[CEBPD[,5]<0.1,]

plot(density(rowMeans(cbind(cebpb[,3],cebpb[,4]))))
lines(density(rowMeans(cbind(CEBPD[,3],CEBPD[,4]))),col="red")

plot(density(rowMeans(cbind(CEBPD[,3],CEBPD[,4]))),col="red")
lines(density(rowMeans(cbind(cebpb[,3],cebpb[,4]))),col="black")
