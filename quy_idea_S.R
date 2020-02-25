library('TFregulomeR')

     cebpb_hepg2 <- loadPeaks(id = "MM1_HSA_HepG2_CEBPB", includeMotifOnly = TRUE)
     write.table(cebpb_hepg2,"HepG2_cebpb.bed",sep="\t",quote=F,row.names=F,col.names=F)
     
     cebpb_atf4_hepg2 <- commonPeaks(target_peak_id="MM1_HSA_HepG2_CEBPB",
                                       motif_only_for_target_peak=TRUE,
                                       compared_peak_id="MM1_HSA_HepG2_ATF4",
                                       motif_only_for_compared_peak=TRUE,
                                       methylation_profile_in_narrow_region=FALSE)
                                       
     cebpb_cebpd_hepg2 <- commonPeaks(target_peak_id="MM1_HSA_HepG2_CEBPB",
                                       motif_only_for_target_peak=TRUE,
                                       compared_peak_id="MM1_HSA_HepG2_CEBPD",
                                       motif_only_for_compared_peak=TRUE,
                                       methylation_profile_in_narrow_region=FALSE)
                                       
cebpb_atf4_peaks <-cebpb_atf4_hepg2[,1][[1]]@common_peak # 1027
cebpb_atf4_peaks <-cebpb_cebpd_hepg2[,1][[1]]@common_peak # 4560



#########
atf4.o<-read.table(pipe('grep -e "site\tMM1_HSA_K562_CEBPB_+_MM1_HSA_K562_ATF4" matrix-scan_2020-01-17_cebpb_atf4_peaks.ft|cut -f4,5,6,9'),
               sep="\t",stringsAsFactors=FALSE)

atf4<-atf4.o[atf4.o[,4]<0.00001,]

atf4[,4]<-(-log10(atf4[,4]))

center_pos <- rowMeans(atf4[,2:3])

atf4_box <- data.frame(pos = center_pos, pval = atf4[,4])

binMap <- cut(atf4_box$pos, 
              breaks = seq(0, 200, by = 10), 
              labels = seq(1, 200, by = 10)
             )

par(mfrow=c(1,1))

boxplot(atf4_box$pval~binMap,ylab="-log10 Pval",xaxt="n",main="CEBPB+ATF4",ylim=c(3,8),col="blue")
axis(1, labels=seq(-100, 100, by = 50) ,at = c(1,5,10,15,20), las=1)

#atf4.t<-aggregate(atf4_box$pval, list(atf4_box$pos), mean)
###############################

cebpb.o<-read.table(pipe('grep -e "site\tCEBPB\t" matrix-scan_2020-01-17_cebpb_atf4_peaks.ft|cut -f4,5,6,9'),
               sep="\t",stringsAsFactors=FALSE)

cebpb<-cebpb.o[cebpb.o[,4]<0.00001,]

cebpb[,4]<-(-log10(cebpb[,4]))

center_pos <- rowMeans(cebpb[,2:3])

cebpb_box <- data.frame(pos = center_pos, pval = cebpb[,4])

binMap <- cut(cebpb_box$pos, 
              breaks = seq(0, 200, by = 10), 
              labels = seq(1, 200, by = 10)
             )

boxplot(cebpb_box$pval~binMap,ylab="-log10 Pval",xaxt="n",main="CEBPB",ylim=c(3,8),add=TRUE,col="red")
axis(1, labels=seq(-100, 100, by = 50) ,at = c(1,5,10,15,20), las=1)

#cebpb.t<-aggregate(cebpb_box$pval, list(cebpb_box$pos), mean)




################

################
atf4<-read.table(pipe('grep -e "site\tMM1_HSA_K562_CEBPB_+_MM1_HSA_K562_ATF4" matrix-scan_2020-01-17_cebpb_atf4_peaks.ft|cut -f4,5,6,9'),
               sep="\t",stringsAsFactors=FALSE)
atf4[,4]<-(-log10(atf4[,4]))
pairs<-rep(1:2, dim(atf4)[1]/2)
center_pos <- rowMeans(atf4[,2:3])
atf4_box <- data.frame(pos = center_pos[pairs==1], pval = center_pos[pairs==1])

p1<-atf4[pairs==1,4]
p2<-atf4[pairs==2,4]

atf4_box$pval<-pmax(p1, p2)
atf4.t<-aggregate(atf4_box$pval, list(atf4_box$pos), mean)

################
atf4<-read.table(pipe('grep -e "site\tMM1_HSA_K562_CEBPB_+_MM1_HSA_K562_CEBPD" matrix-scan_2020-01-17_cebpb_atf4_peaks.ft|cut -f4,5,6,9'),
               sep="\t",stringsAsFactors=FALSE)
atf4[,4]<-(-log10(atf4[,4]))
pairs<-rep(1:2, dim(atf4)[1]/2)
center_pos <- rowMeans(atf4[,2:3])
atf4_box <- data.frame(pos = center_pos[pairs==1], pval = center_pos[pairs==1])

p1<-atf4[pairs==1,4]
p2<-atf4[pairs==2,4]

atf4_box$pval<-pmax(p1, p2)
cebpd.t<-aggregate(atf4_box$pval, list(atf4_box$pos), mean)
################
pdf("mean_pvals.pdf")
plot(as.numeric(cebpd.t[,1]),cebpd.t[,2],col="white",main="Peaks overlapping between CEBPB and ATF4 ",
      ylab = "Mean Pval",xlab="Peak center",
     yaxt="n",xaxt="n",xlim=c(1,201), ylim=c(0,5) )
axis(1, labels=seq(-100, 100, by = 50) ,at = seq(1, 201, by = 50), las=1)
axis(2, labels=seq(0, 200, by = 25) ,at = seq(0, 200, by = 25), las=1)
lines(smooth.spline(as.numeric(cebpd.t[,1]),(cebpd.t[,2]),spar=.3),lwd=3)
lines(smooth.spline(as.numeric(atf4.t[,1]),(atf4.t[,2]),spar=.3),lwd=3,col="blue")
legend("topright",legend=c("CEBPB+ATF4 matrix","CEBPB+CEBPD matrix"),fill=c("blue","black"),bty = "n")
dev.off()

pdf("sum_pvals.pdf")
plot(as.numeric(cebpd.t[,1]),cebpd.t[,2],col="white",main="Peaks overlapping between CEBPB and ATF4 ",
      ylab = "Sum Pvals",xlab="Peak center",
     yaxt="n",xaxt="n",xlim=c(1,201), ylim=c(1000,2500) )
axis(1, labels=seq(-100, 100, by = 50) ,at = seq(1, 201, by = 50), las=1)
lines(smooth.spline(as.numeric(cebpd.t[,1]),(cebpd.t[,2]),spar=.3),lwd=3)
lines(smooth.spline(as.numeric(atf4.t[,1]),(atf4.t[,2]),spar=.3),lwd=3,col="blue")
legend("topright",legend=c("CEBPB+ATF4 matrix","CEBPB+CEBPD matrix"),fill=c("blue","black"),bty = "n")
dev.off()




################################


################
atf4<-read.table(pipe('grep -e "site\tCEBPB\t" matrix-scan_CEBPB+CEBPD.tf|cut -f4,5,6,9'),
               sep="\t",stringsAsFactors=FALSE)
atf4[,4]<-(-log10(atf4[,4]))
pairs<-rep(1:2, dim(atf4)[1]/2)
center_pos <- rowMeans(atf4[,2:3])
atf4_box <- data.frame(pos = center_pos[pairs==1], pval = center_pos[pairs==1])

p1<-atf4[pairs==1,4]
p2<-atf4[pairs==2,4]

atf4_box$pval<-pmax(p1, p2)
cebpb.t<-aggregate(atf4_box$pval, list(atf4_box$pos), mean)

################
atf4<-read.table(pipe('grep -e "site\tMM1_HSA_K562_CEBPB_+_MM1_HSA_K562_CEBPD" matrix-scan_CEBPB+CEBPD.tf|cut -f4,5,6,9'),
               sep="\t",stringsAsFactors=FALSE)
atf4[,4]<-(-log10(atf4[,4]))
pairs<-rep(1:2, dim(atf4)[1]/2)
center_pos <- rowMeans(atf4[,2:3])
atf4_box <- data.frame(pos = center_pos[pairs==1], pval = center_pos[pairs==1])

p1<-atf4[pairs==1,4]
p2<-atf4[pairs==2,4]

atf4_box$pval<-pmax(p1, p2)
cebpd.t<-aggregate(atf4_box$pval, list(atf4_box$pos), mean)
################
