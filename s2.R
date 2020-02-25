library('TFregulomeR')

     
     cebpb_atf4 <- commonPeaks(target_peak_id="MM1_HSA_K562_CEBPB",
                                       motif_only_for_target_peak=TRUE,
                                       compared_peak_id="MM1_HSA_K562_ATF4",
                                       motif_only_for_compared_peak=TRUE,
                                       methylation_profile_in_narrow_region=FALSE)
                                       
     cebpb_cebpd <- commonPeaks(target_peak_id="MM1_HSA_K562_CEBPB",
                                       motif_only_for_target_peak=TRUE,
                                       compared_peak_id="MM1_HSA_K562_CEBPD",
                                       motif_only_for_compared_peak=TRUE,
                                       methylation_profile_in_narrow_region=FALSE)
                                       
     common <- commonPeaks(target_peak_id="MM1_HSA_K562_CEBPB",
                                       motif_only_for_target_peak=TRUE,
                                       compared_peak_id=c("MM1_HSA_K562_ATF4","MM1_HSA_K562_CEBPD"),
                                       motif_only_for_compared_peak=TRUE,
                                       methylation_profile_in_narrow_region=FALSE)
                                       
    common1 <- commonPeaks(target_peak_id="MM1_HSA_K562_ATF4",
                                       motif_only_for_target_peak=TRUE,
                                       compared_peak_id=c("MM1_HSA_K562_CEBPD"),
                                       motif_only_for_compared_peak=TRUE,
                                       methylation_profile_in_narrow_region=FALSE)
                                       
# only a 133 peaks overlap between the 3
# 138 ATF4 & CEBPD

cebpb_atf4_peaks <-cebpb_atf4[,1][[1]]@common_peak # 2560
cebpb_atf4_peaks[,2] <- cebpb_atf4_peaks[,2]-100
cebpb_atf4_peaks[,3] <- cebpb_atf4_peaks[,3]+100

cebpb_cebpd_peaks <-cebpb_cebpd[,1][[1]]@common_peak # 517
cebpb_cebpd_peaks[,2] <- cebpb_cebpd_peaks[,2]-100
cebpb_cebpd_peaks[,3] <- cebpb_cebpd_peaks[,3]+100

write.table(cebpb_atf4_peaks,"cebpb_atf4_peaks.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(cebpb_cebpd_peaks,"cebpb_cebpd_peaks.bed",sep="\t",quote=F,row.names=F,col.names=F)

# FETCH SEQ RSAT
system("wget http://rsat.sb-roscoff.fr//tmp/apache/2019/11/28/cebpb_atf4_peaks_28m5_20191128_180643.fasta")
system("wget http://rsat.sb-roscoff.fr//tmp/apache/2019/11/28/cebpb_cebpd_peaks_m448_20191128_180734.fasta")

# MATRIX SCAN RSAT CEBPB+ATF4
ms<-read.table(pipe('grep "site" matrix-scan_CEBPB+ATF4.ft|cut -f3,4,5,6'),sep="\t",stringsAsFactors=FALSE)
pos <- rowMeans(cbind(ms[,4],ms[,3]))
pos.t <- table(pos)

cebpb<-pos[ms[,1]=="CEBPB"]
cebpb.t<-table(cebpb)

atf4<-pos[ms[,1]=="MM1_HSA_K562_CEBPB_+_MM1_HSA_K562_ATF4"]
atf4.t<-table(atf4)

jund<-pos[ms[,1]=="MM1_HSA_K562_CEBPB_+_MM1_HSA_K562_CEBPD"]
jund.t<-table(jund)

pdf("CEBPB_and_ATF4_matrixScan.pdf")
plot(as.numeric(names(pos.t)),pos.t,col="white",main="Peaks overlapping between CEBPB and ATF4 ",
      ylab = "# Predicted binding sites",xlab="Peak center",
     yaxt="n",xaxt="n",xlim=c(1,201), ylim=c(0,max(cebpb.t,atf4.t,jund.t)) )
axis(1, labels=seq(-100, 100, by = 50) ,at = seq(1, 201, by = 50), las=1)
axis(2, labels=seq(0, 200, by = 25) ,at = seq(0, 200, by = 25), las=1)
lines(smooth.spline(as.numeric(names(cebpb.t)),cebpb.t,spar=.3),lwd=3)
lines(smooth.spline(as.numeric(names(atf4.t)),atf4.t,spar=.3),lwd=3,col="blue")
lines(smooth.spline(as.numeric(names(jund.t)),jund.t,spar=.3),lwd=3,col="red")
legend("topright",legend=c("CEBPB matrix","CEBPB+ATF4 matrix","CEBPB+CEBPD matrix"),fill=c("black","blue","red"),bty = "n")
dev.off()

# MATRIX SCAN RSAT CEBPB+CEBPD
ms<-read.table(pipe('grep "site" matrix-scan_CEBPB+CEBPD.tf|cut -f3,4,5,6'),sep="\t",stringsAsFactors=FALSE)
pos <- rowMeans(cbind(ms[,4],ms[,3]))
pos.t <- table(pos)

cebpb<-pos[ms[,1]=="CEBPB"]
cebpb.t<-table(cebpb)

atf4<-pos[ms[,1]=="MM1_HSA_K562_CEBPB_+_MM1_HSA_K562_ATF4"]
atf4.t<-table(atf4)

jund<-pos[ms[,1]=="MM1_HSA_K562_CEBPB_+_MM1_HSA_K562_CEBPD"]
jund.t<-table(jund)

pdf("CEBPB_and_CEBPD_matrixScan.pdf")
plot(as.numeric(names(pos.t)),pos.t,col="white",main="Peaks overlapping between CEBPB and CEBPD ",
      ylab = "# Predicted binding sites",xlab="Peak center",
     yaxt="n",xaxt="n",xlim=c(1,201), ylim=c(0,max(cebpb.t,atf4.t,jund.t)) )
axis(1, labels=seq(-100, 100, by = 50) ,at = seq(1, 201, by = 50), las=1)
axis(2, labels=seq(0, 100, by = 15) ,at = seq(0, 100, by = 15), las=1)
lines(smooth.spline(as.numeric(names(cebpb.t)),cebpb.t,spar=.3),lwd=3)
lines(smooth.spline(as.numeric(names(atf4.t)),atf4.t,spar=.3),lwd=3,col="blue")
lines(smooth.spline(as.numeric(names(jund.t)),jund.t,spar=.3),lwd=3,col="red")
legend("topright",legend=c("CEBPB matrix","CEBPB+ATF4 matrix","CEBPB+CEBPD matrix"),fill=c("black","blue","red"),bty = "n")
dev.off()

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# keso

# MATRIX SCAN RSAT CEBPB+ATF4
ms<-read.table(pipe('grep "site" matrix-scan_CEBPB+ATF4.ft|cut -f3,4,5,6'),sep="\t",stringsAsFactors=FALSE)
pos <- rowMeans(cbind(ms[,4],ms[,3]))
pos.t <- table(pos)

cebpb<-pos[ms[,1]=="CEBPB"]
cebpb.t<-table(cebpb)

atf4<-pos[ms[,1]=="MM1_HSA_K562_CEBPB_+_MM1_HSA_K562_ATF4"]
atf4.t<-table(atf4)

jund<-pos[ms[,1]=="ATF4"]
jund.t<-table(jund)

meco<-pos[ms[,1]=="CEBPD"]
meco.t<-table(meco)

pdf("CEBPB_and_ATF4_individualMATRIX_matrixScan.pdf")
plot(as.numeric(names(pos.t)),pos.t,col="white",main="Peaks overlapping between CEBPB and ATF4 ",
      ylab = "# Predicted binding sites",xlab="Peak center",
     yaxt="n",xaxt="n",xlim=c(1,201), ylim=c(0,max(cebpb.t,atf4.t,jund.t)) )
axis(1, labels=seq(-100, 100, by = 50) ,at = seq(1, 201, by = 50), las=1)
axis(2, labels=seq(0, 200, by = 25) ,at = seq(0, 200, by = 25), las=1)
lines(smooth.spline(as.numeric(names(cebpb.t)),cebpb.t,spar=.3),lwd=3)
lines(smooth.spline(as.numeric(names(atf4.t)),atf4.t,spar=.3),lwd=3,col="blue")
lines(smooth.spline(as.numeric(names(jund.t)),jund.t,spar=.3),lwd=3,col="darkgoldenrod")
lines(smooth.spline(as.numeric(names(meco.t)),meco.t,spar=.3),lwd=3,col="darkgreen")
legend("topright",legend=c("CEBPB matrix","CEBPB+ATF4 matrix","ATF4 matrix","CEBPD matrix"),fill=c("black","blue","darkgoldenrod","darkgreen"),bty = "n")
dev.off()

# MATRIX SCAN RSAT CEBPB+CEBPD
ms<-read.table(pipe('grep "site" matrix-scan_CEBPB+CEBPD.tf|cut -f3,4,5,6'),sep="\t",stringsAsFactors=FALSE)
pos <- rowMeans(cbind(ms[,4],ms[,3]))
pos.t <- table(pos)

cebpb<-pos[ms[,1]=="CEBPB"]
cebpb.t<-table(cebpb)

atf4<-pos[ms[,1]=="CEBPD"]
atf4.t<-table(atf4)

jund<-pos[ms[,1]=="MM1_HSA_K562_CEBPB_+_MM1_HSA_K562_CEBPD"]
jund.t<-table(jund)

meco<-pos[ms[,1]=="ATF4"]
meco.t<-table(meco)

pdf("CEBPB_and_CEBPD_individualMATRIX_matrixScan.pdf")
plot(as.numeric(names(pos.t)),pos.t,col="white",main="Peaks overlapping between CEBPB and CEBPD ",
      ylab = "# Predicted binding sites",xlab="Peak center",
     yaxt="n",xaxt="n",xlim=c(1,201), ylim=c(0,max(cebpb.t,atf4.t,jund.t)) )
axis(1, labels=seq(-100, 100, by = 50) ,at = seq(1, 201, by = 50), las=1)
axis(2, labels=seq(0, 100, by = 15) ,at = seq(0, 100, by = 15), las=1)
lines(smooth.spline(as.numeric(names(cebpb.t)),cebpb.t,spar=.3),lwd=3)
lines(smooth.spline(as.numeric(names(atf4.t)),atf4.t,spar=.3),lwd=3,col="darkgreen")
lines(smooth.spline(as.numeric(names(jund.t)),jund.t,spar=.3),lwd=3,col="red")
lines(smooth.spline(as.numeric(names(meco.t)),meco.t,spar=.3),lwd=3,col="darkgoldenrod")
legend("topright",legend=c("CEBPB matrix","CEBPB+CEBPD matrix","ATF4 matrix","CEBPD matrix"),fill=c("black","red","darkgoldenrod","darkgreen"),bty = "n")
dev.off()
