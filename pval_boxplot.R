###############################

atf4.o<-read.table(pipe('grep -e "site\tMM1_HSA_K562_CEBPB_+_MM1_HSA_K562_ATF4" matrix-scan_2020-01-17_cebpb_atf4_peaks.ft|cut -f4,5,6,9'),
               sep="\t",stringsAsFactors=FALSE)

atf4<-atf4.o[atf4.o[,4]<0.001,]

atf4[,4]<-(-log10(atf4[,4]))

center_pos <- rowMeans(atf4[,2:3])

atf4_box <- data.frame(pos = center_pos, pval = atf4[,4])

#binMap <- cut(atf4_box$pos, 
#              breaks = seq(0, 200, by = 10), 
#              labels = seq(1, 200, by = 10)
#             )

atf4_binMap <- cut(atf4_box$pos, 
              breaks = seq(50, 150, by = 10), 
              labels = seq(51, 150, by = 10)
             )
###############################

ttff.o<-read.table(pipe('grep -e "site\tATF4\t" matrix-scan_2020-01-17_cebpb_atf4_peaks.ft|cut -f4,5,6,9'),
               sep="\t",stringsAsFactors=FALSE)

ttff<-ttff.o[ttff.o[,4]<0.001,]

ttff[,4]<-(-log10(ttff[,4]))

center_pos <- rowMeans(ttff[,2:3])

ttff_box <- data.frame(pos = center_pos, pval = ttff[,4])

ttff_binMap <- cut(ttff_box$pos, 
              breaks = seq(50, 150, by = 10), 
              labels = seq(51, 150, by = 10)
             )
###############################

cebpb.o<-read.table(pipe('grep -e "site\tCEBPB\t" matrix-scan_2020-01-17_cebpb_atf4_peaks.ft|cut -f4,5,6,9'),
               sep="\t",stringsAsFactors=FALSE)

cebpb<-cebpb.o[cebpb.o[,4]<0.001,]

cebpb[,4]<-(-log10(cebpb[,4]))

center_pos <- rowMeans(cebpb[,2:3])

cebpb_box <- data.frame(pos = center_pos, pval = cebpb[,4])

cebpb_binMap <- cut(cebpb_box$pos, 
              breaks = seq(50, 150, by = 10), 
              labels = seq(51, 150, by = 10)
             )

###############################
atf4_box_dim<-data.frame(atf4_box,Matrix=rep("CEBPB+ATF4",dim(atf4_box)[1]) )
cebpb_box_dim<-data.frame(cebpb_box,Matrix=rep("CEBPB",dim(cebpb_box)[1]) )
ttff_box_dim<-data.frame(ttff_box,Matrix=rep("ATF4",dim(ttff_box)[1]) )


atf4_box_dim$pos <- atf4_binMap
cebpb_box_dim$pos <- cebpb_binMap
ttff_box_dim$pos <- ttff_binMap

databox <- rbind(atf4_box_dim,cebpb_box_dim,ttff_box_dim)
databox <- databox[!is.na(databox[,1]),]

library(ggplot2)

pdf("matrixScanpvalue_ATF4_ttff.pdf",height=3,width=8)
ggplot(databox, aes(x = pos, y = pval, fill = Matrix)) + geom_boxplot(outlier.shape = NA) + xlab("Distance from CEBPB and ATF4 overlapping peaks") + ylab("-Log10 Pvalues")+ 
scale_x_discrete(labels=c("51" = "-50 -40", 
                          "61" = "-40 -30",
                          "71" = "-30 -20",
                          "81" = "-20 -10",
                          "91" = "-10 0",
                          "101" = "0 +10",
                          "111" = "+10 +20",
                          "121" = "+20 +30",
                          "131" = "+30 +40",
                          "141" = "+40 +50"))
dev.off()

############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
cebpd.o<-read.table(pipe('grep -e "site\tMM1_HSA_K562_CEBPB_+_MM1_HSA_K562_CEBPD" ~/CSI/FPWM/tests/matrix-scan_CEBPB+CEBPD.tf|cut -f4,5,6,9'),
               sep="\t",stringsAsFactors=FALSE)

cebpd<-cebpd.o[cebpd.o[,4]<0.0001,]

cebpd[,4]<-(-log10(cebpd[,4]))

center_pos <- rowMeans(cebpd[,2:3])

cebpd_box <- data.frame(pos = center_pos, pval = cebpd[,4])

cebpd_binMap <- cut(cebpd_box$pos, 
              breaks = seq(50, 150, by = 10), 
              labels = seq(51, 150, by = 10)
             )
#################################################
cebpb.o<-read.table(pipe('grep -e "site\tCEBPB\t" ~/CSI/FPWM/tests/matrix-scan_CEBPB+CEBPD.tf|cut -f4,5,6,9'),sep="\t",stringsAsFactors=FALSE)

cebpb<-cebpb.o[cebpb.o[,4]<0.0001,]

cebpb[,4]<-(-log10(cebpb[,4]))

center_pos <- rowMeans(cebpb[,2:3])

cebpb_box <- data.frame(pos = center_pos, pval = cebpb[,4])

cebpb_binMap <- cut(cebpb_box$pos, 
              breaks = seq(50, 150, by = 10), 
              labels = seq(51, 150, by = 10)
             )

###############################
cebpd_box_dim<-data.frame(cebpd_box,Matrix=rep("CEBPB+CEBPD",dim(cebpd_box)[1]) )
cebpb_box_dim<-data.frame(cebpb_box,Matrix=rep("CEBPB",dim(cebpb_box)[1]) )

cebpd_box_dim$pos <- cebpd_binMap
cebpb_box_dim$pos <- cebpb_binMap

databox <- rbind(cebpd_box_dim,cebpb_box_dim)
databox <- databox[!is.na(databox[,1]),]

library(ggplot2)

pdf("matrixScanpvalue_CEBPD.pdf",height=3,width=8)
ggplot(databox, aes(x = pos, y = pval, fill = Matrix)) + geom_boxplot(outlier.shape = NA) + xlab("Distance from CEBPB and CEBPD overlapping peaks") + ylab("-Log10 Pvalues")+ 
scale_x_discrete(labels=c("51" = "-50 -40", 
                          "61" = "-40 -30",
                          "71" = "-30 -20",
                          "81" = "-20 -10",
                          "91" = "-10 0",
                          "101" = "0 +10",
                          "111" = "+10 +20",
                          "121" = "+20 +30",
                          "131" = "+30 +40",
                          "141" = "+40 +50"))
dev.off()
