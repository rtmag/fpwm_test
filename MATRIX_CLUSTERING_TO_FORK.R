library('TFregulomeR')
cebpa <- loadPeaks(id = "MM1_HSA_HepG2_CEBPA", includeMotifOnly = TRUE)
cebpa[,2] <- cebpa[,2]-50
cebpa[,3] <- cebpa[,3]+50
write.table(cebpa[,1:3],"HepG2_CEBPA.bed",sep="\t",quote=F,row.names=F,col.names=F)
################################################################################
################################################################################
################################################################################
# Example 1 MAFF K562
FPWM_MAFF <- createFPWM(mainTF ="MAFF",
                            partners = c("MAFK","NFE2"),
                                cell = "K562", 
                        forkPosition = 7)

write.FPWM( FPWM = FPWM_MAFF, 
            format = "transfac", 
            fileName = "FPWM_Counts_MAFF+MAFK_NFE2+K562.transfac")

FPWM_MAFF_p <- createFPWM(mainTF ="MAFF",
                            partners = c("MAFK","NFE2"),
                                cell = "K562", 
                        forkPosition = 7,
                       probabilityMatrix = TRUE)

write.FPWM( FPWM = FPWM_MAFF_p, 
            format = "transfac", 
            fileName = "FPWM_Probability_MAFF+MAFK_NFE2+K562.transfac")

MAFF_MAFK_K562 <- commonPeaks(target_peak_id="MM1_HSA_K562_MAFF",
                                       motif_only_for_target_peak=TRUE,
                                       compared_peak_id="MM1_HSA_K562_MAFK",
                                       motif_only_for_compared_peak=TRUE,
                                       methylation_profile_in_narrow_region=FALSE)
MAFF_MAFK_K562_peaks <-MAFF_MAFK_K562[,1][[1]]@common_peak # 9790
write.table(MAFF_MAFK_K562_peaks[,1:3],"MAFF_MAFK_K562_9790peaks.bed",sep="\t",quote=F,row.names=F,col.names=F)


MAFF_NFE2_K562 <- commonPeaks(target_peak_id="MM1_HSA_K562_MAFF",
                                       motif_only_for_target_peak=TRUE,
                                       compared_peak_id="MM1_HSA_K562_NFE2",
                                       motif_only_for_compared_peak=TRUE,
                                       methylation_profile_in_narrow_region=FALSE)
MAFF_NFE2_K562_peaks <-MAFF_NFE2_K562[,1][[1]]@common_peak # 4821
write.table(MAFF_NFE2_K562_peaks[,1:3],"MAFF_NEF2_K562_4821peaks.bed",sep="\t",quote=F,row.names=F,col.names=F)





# Example 2 JUND HepG2
FPWM_JUND <- createFPWM(mainTF ="JUND",
                            partners = c("ATF2","FOSL2"),
                                cell = "HepG2", 
                        forkPosition = 6)

write.FPWM( FPWM = FPWM_JUND, 
            format = "transfac", 
            fileName = "FPWM_Counts_JUND+ATF2_FOSL2+HepG2.transfac")

FPWM_JUND_p <- createFPWM(mainTF ="JUND",
                            partners = c("ATF2","FOSL2"),
                                cell = "HepG2", 
                        forkPosition = 6,
                       probabilityMatrix = TRUE)

write.FPWM( FPWM = FPWM_JUND_p, 
            format = "transfac", 
            fileName = "FPWM_Probability_JUND+ATF2_FOSL2+HepG2.transfac")


JUND_ATF2_HepG2 <- commonPeaks(target_peak_id="MM1_HSA_HepG2_JUND",
                                       motif_only_for_target_peak=TRUE,
                                       compared_peak_id="MM1_HSA_HepG2_ATF2",
                                       motif_only_for_compared_peak=TRUE,
                                       methylation_profile_in_narrow_region=FALSE)
JUND_ATF2_HepG2_peaks <-JUND_ATF2_HepG2[,1][[1]]@common_peak # 6325
write.table(JUND_ATF2_HepG2_peaks[,1:3],"JUND_ATF2_HepG2_6325peaks.bed",sep="\t",quote=F,row.names=F,col.names=F)


JUND_FOSL2_HepG2 <- commonPeaks(target_peak_id="MM1_HSA_HepG2_JUND",
                                       motif_only_for_target_peak=TRUE,
                                       compared_peak_id="MM1_HSA_HepG2_FOSL2",
                                       motif_only_for_compared_peak=TRUE,
                                       methylation_profile_in_narrow_region=FALSE)
JUND_FOSL2_HepG2_peaks <-JUND_FOSL2_HepG2[,1][[1]]@common_peak # 3160
write.table(JUND_FOSL2_HepG2_peaks[,1:3],"JUND_FOSL2_HepG2_3160peaks.bed",sep="\t",quote=F,row.names=F,col.names=F)





K562_CEBPB_transfac <- searchMotif(id = "MM1_HSA_K562_CEBPB", motif_format = "TRANSFAC")

K562_CEBPB_transfac@MMmotif@motif_matrix <- apply(K562_CEBPB_transfac@MMmotif@motif_matrix,2,function(x) x/sum(x))

exportMMPFM(fun_output = K562_CEBPB_transfac, fun = "searchMotif",
                 save_motif_PFM = TRUE, save_betaScore_matrix = FALSE)
                                                  
                                                  
jund_walt <- intersectPeakMatrix(peak_id_x = "MM1_HSA_HepG2_JUND",
                                              peak_id_y = c("MM1_HSA_HepG2_ATF2","MM1_HSA_HepG2_FOSL2"),
                                              motif_only_for_id_y = TRUE, 
                                              methylation_profile_in_narrow_region = TRUE)
                                                  
                                                  

                                                  
exportMMPFM(fun_output = jund_walt, 
            fun = "intersectPeakMatrix", 
            save_motif_PFM = TRUE, 
            save_betaScore_matrix = FALSE)
                                                  
                                                  
jund_walt <- intersectPeakMatrix(peak_id_x = "MM1_HSA_HepG2_JUND",
                                              peak_id_y = c("MM1_HSA_HepG2_ATF2","MM1_HSA_HepG2_FOSL2"),
                                              motif_only_for_id_y = TRUE, 
                                              methylation_profile_in_narrow_region = TRUE)
                                                  
fileConn <- file("jund_walt.transfac")

		transfac_vector <- c()
		transfac_vector <- c( transfac_vector, paste0("AC ",jund_walt[1,1][[1]]@id) )
		transfac_vector <- c( transfac_vector, "XX" )
		transfac_vector <- c( transfac_vector, paste0("ID ",jund_walt[1,1][[1]]@id), "XX" )
		transfac_vector <- c( transfac_vector, paste0("DE ")  )
		transfac_vector <- c( transfac_vector, paste("PO","A","C","G","T",sep="\t")  )
    
    # probability to counts:
    matrix <- jund_walt[1,1][[1]]@MethMotif_x@MMmotif@motif_matrix * jund_walt[1,1][[1]]@MethMotif_x@MMmotif@nsites   
                                                  
    for ( jx in 1:dim(matrix)[1] ){
				transfac_vector <- c( transfac_vector, paste(matrix[jx,],collapse="\t") )
			}
                                                  
		transfac_vector <- c( transfac_vector, c("XX","CC program: FPWM")  )
		transfac_vector <- c( transfac_vector, paste0("CC numberOfSites: ",jund_walt[1,1][[1]]@MethMotif_x@MMmotif@nsites)  )
		transfac_vector <- c( transfac_vector, c("XX","//")  )

    writeLines(transfac_vector, fileConn)
    close(fileConn)

						  
						  
######################################################################################################################################################				  
##### CEBPB + ATF4; CEBPD in K562			  
# CEBPB + ATF4 peaks
CEBPB_ATF4_K562 <- commonPeaks(target_peak_id="MM1_HSA_K562_CEBPB",
                                       motif_only_for_target_peak=TRUE,
                                       compared_peak_id="MM1_HSA_K562_ATF4",
                                       motif_only_for_compared_peak=TRUE,
                                       methylation_profile_in_narrow_region=FALSE)
CEBPB_ATF4_K562_peaks <-CEBPB_ATF4_K562[,1][[1]]@common_peak # 2560
write.table(CEBPB_ATF4_K562_peaks[,1:3],"CEBPB_ATF4_K562_2560peaks.bed",sep="\t",quote=F,row.names=F,col.names=F)

# CEBPB + CEBPD peaks
CEBPB_CEBPD_K562 <- commonPeaks(target_peak_id="MM1_HSA_K562_CEBPB",
                                       motif_only_for_target_peak=TRUE,
                                       compared_peak_id="MM1_HSA_K562_CEBPD",
                                       motif_only_for_compared_peak=TRUE,
                                       methylation_profile_in_narrow_region=FALSE)
CEBPB_CEBPD_K562_peaks <-CEBPB_CEBPD_K562[,1][[1]]@common_peak # 517
write.table(CEBPB_CEBPD_K562_peaks[,1:3],"CEBPB_CEBPD_K562_517peaks.bed",sep="\t",quote=F,row.names=F,col.names=F)
						  
###########################################################################
cebpb_inter <- intersectPeakMatrix(peak_id_x = "MM1_HSA_K562_CEBPB",
                                              peak_id_y = c("MM1_HSA_K562_CEBPB","MM1_HSA_K562_ATF4","MM1_HSA_K562_CEBPD"),
                                              motif_only_for_id_y = TRUE,
                                              methylation_profile_in_narrow_region = TRUE)
						  
fileConn <- file("K562_CEBPB_AND_ATF4_CEBPD.transfac")
		transfac_vector <- c()
		transfac_vector <- c( transfac_vector, paste0("AC ",cebpb_inter[1,1][[1]]@id) )
		transfac_vector <- c( transfac_vector, "XX" )
		transfac_vector <- c( transfac_vector, paste0("ID ",cebpb_inter[1,1][[1]]@id), "XX" )
		transfac_vector <- c( transfac_vector, paste0("DE ")  )
		transfac_vector <- c( transfac_vector, paste("PO","A","C","G","T",sep="\t")  )
                                         
    matrix <- cebpb_inter[1,1][[1]]@MethMotif_x@MMmotif@motif_matrix * cebpb_inter[1,1][[1]]@MethMotif_x@MMmotif@nsites   

    for ( jx in 1:dim(matrix)[1] ){
				transfac_vector <- c( transfac_vector, paste(matrix[jx,],collapse="\t") )
			}
						  
		transfac_vector <- c( transfac_vector, c("XX","CC program: forkedTF")  )
		transfac_vector <- c( transfac_vector, paste0("CC numberOfSites: ",cebpb_inter[1,1][[1]]@MethMotif_x@MMmotif@nsites)  )
		transfac_vector <- c( transfac_vector, paste0("CC numberOfPeaks: ",cebpb_inter[1,1][[1]]@MethMotif_x@MMmotif@nPeaks)  )
		transfac_vector <- c( transfac_vector, c("XX","//")  )
################						  
		transfac_vector <- c( transfac_vector, paste0("AC ",cebpb_inter[1,2][[1]]@id) )
		transfac_vector <- c( transfac_vector, "XX" )
		transfac_vector <- c( transfac_vector, paste0("ID ",cebpb_inter[1,2][[1]]@id), "XX" )
		transfac_vector <- c( transfac_vector, paste0("DE ")  )
		transfac_vector <- c( transfac_vector, paste("PO","A","C","G","T",sep="\t")  )
                                         
    matrix <- cebpb_inter[1,2][[1]]@MethMotif_x@MMmotif@motif_matrix * cebpb_inter[1,2][[1]]@MethMotif_x@MMmotif@nsites   

    for ( jx in 1:dim(matrix)[1] ){
				transfac_vector <- c( transfac_vector, paste(matrix[jx,],collapse="\t") )
			}
						  
		transfac_vector <- c( transfac_vector, c("XX","CC program: forkedTF")  )
		transfac_vector <- c( transfac_vector, paste0("CC numberOfSites: ",cebpb_inter[1,2][[1]]@MethMotif_x@MMmotif@nsites)  )
		transfac_vector <- c( transfac_vector, paste0("CC numberOfPeaks: ",cebpb_inter[1,2][[1]]@MethMotif_x@MMmotif@nPeaks)  )
		transfac_vector <- c( transfac_vector, c("XX","//")  )
################
		transfac_vector <- c( transfac_vector, paste0("AC ",cebpb_inter[1,3][[1]]@id) )
		transfac_vector <- c( transfac_vector, "XX" )
		transfac_vector <- c( transfac_vector, paste0("ID ",cebpb_inter[1,3][[1]]@id), "XX" )
		transfac_vector <- c( transfac_vector, paste0("DE ")  )
		transfac_vector <- c( transfac_vector, paste("PO","A","C","G","T",sep="\t")  )
                                         
    matrix <- cebpb_inter[1,3][[1]]@MethMotif_x@MMmotif@motif_matrix * cebpb_inter[1,3][[1]]@MethMotif_x@MMmotif@nsites   

    for ( jx in 1:dim(matrix)[1] ){
				transfac_vector <- c( transfac_vector, paste(matrix[jx,],collapse="\t") )
			}
						  
		transfac_vector <- c( transfac_vector, c("XX","CC program: forkedTF")  )
		transfac_vector <- c( transfac_vector, paste0("CC numberOfSites: ",cebpb_inter[1,3][[1]]@MethMotif_x@MMmotif@nsites)  )
		transfac_vector <- c( transfac_vector, paste0("CC numberOfPeaks: ",cebpb_inter[1,3][[1]]@MethMotif_x@MMmotif@nPeaks)  )
		transfac_vector <- c( transfac_vector, c("XX","//")  )						  
						  
						  
    writeLines(transfac_vector, fileConn)
    close(fileConn)
						  
						  
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################

cebpb_inter <- intersectPeakMatrix(peak_id_x = "MM1_HSA_HepG2_JUND",
                                              peak_id_y = c("MM1_HSA_HepG2_JUND","MM1_HSA_HepG2_ATF2","MM1_HSA_HepG2_FOSL2"),
                                              motif_only_for_id_y = TRUE,
                                              methylation_profile_in_narrow_region = TRUE)
						  
fileConn <- file("HepG2_JUND_AND_ATF2_FOSL2.transfac")
		transfac_vector <- c()
		transfac_vector <- c( transfac_vector, paste0("AC ",cebpb_inter[1,1][[1]]@id) )
		transfac_vector <- c( transfac_vector, "XX" )
		transfac_vector <- c( transfac_vector, paste0("ID ",cebpb_inter[1,1][[1]]@id), "XX" )
		transfac_vector <- c( transfac_vector, paste0("DE ")  )
		transfac_vector <- c( transfac_vector, paste("PO","A","C","G","T",sep="\t")  )
                                         
    matrix <- cebpb_inter[1,1][[1]]@MethMotif_x@MMmotif@motif_matrix * cebpb_inter[1,1][[1]]@MethMotif_x@MMmotif@nsites   

    for ( jx in 1:dim(matrix)[1] ){
				transfac_vector <- c( transfac_vector, paste(matrix[jx,],collapse="\t") )
			}
						  
		transfac_vector <- c( transfac_vector, c("XX","CC program: forkedTF")  )
		transfac_vector <- c( transfac_vector, paste0("CC numberOfSites: ",cebpb_inter[1,1][[1]]@MethMotif_x@MMmotif@nsites)  )
		transfac_vector <- c( transfac_vector, paste0("CC numberOfPeaks: ",cebpb_inter[1,1][[1]]@MethMotif_x@MMmotif@nPeaks)  )
		transfac_vector <- c( transfac_vector, c("XX","//")  )
################						  
		transfac_vector <- c( transfac_vector, paste0("AC ",cebpb_inter[1,2][[1]]@id) )
		transfac_vector <- c( transfac_vector, "XX" )
		transfac_vector <- c( transfac_vector, paste0("ID ",cebpb_inter[1,2][[1]]@id), "XX" )
		transfac_vector <- c( transfac_vector, paste0("DE ")  )
		transfac_vector <- c( transfac_vector, paste("PO","A","C","G","T",sep="\t")  )
                                         
    matrix <- cebpb_inter[1,2][[1]]@MethMotif_x@MMmotif@motif_matrix * cebpb_inter[1,2][[1]]@MethMotif_x@MMmotif@nsites   

    for ( jx in 1:dim(matrix)[1] ){
				transfac_vector <- c( transfac_vector, paste(matrix[jx,],collapse="\t") )
			}
						  
		transfac_vector <- c( transfac_vector, c("XX","CC program: forkedTF")  )
		transfac_vector <- c( transfac_vector, paste0("CC numberOfSites: ",cebpb_inter[1,2][[1]]@MethMotif_x@MMmotif@nsites)  )
		transfac_vector <- c( transfac_vector, paste0("CC numberOfPeaks: ",cebpb_inter[1,2][[1]]@MethMotif_x@MMmotif@nPeaks)  )
		transfac_vector <- c( transfac_vector, c("XX","//")  )
################
		transfac_vector <- c( transfac_vector, paste0("AC ",cebpb_inter[1,3][[1]]@id) )
		transfac_vector <- c( transfac_vector, "XX" )
		transfac_vector <- c( transfac_vector, paste0("ID ",cebpb_inter[1,3][[1]]@id), "XX" )
		transfac_vector <- c( transfac_vector, paste0("DE ")  )
		transfac_vector <- c( transfac_vector, paste("PO","A","C","G","T",sep="\t")  )
                                         
    matrix <- cebpb_inter[1,3][[1]]@MethMotif_x@MMmotif@motif_matrix * cebpb_inter[1,3][[1]]@MethMotif_x@MMmotif@nsites   

    for ( jx in 1:dim(matrix)[1] ){
				transfac_vector <- c( transfac_vector, paste(matrix[jx,],collapse="\t") )
			}
						  
		transfac_vector <- c( transfac_vector, c("XX","CC program: forkedTF")  )
		transfac_vector <- c( transfac_vector, paste0("CC numberOfSites: ",cebpb_inter[1,3][[1]]@MethMotif_x@MMmotif@nsites)  )
		transfac_vector <- c( transfac_vector, paste0("CC numberOfPeaks: ",cebpb_inter[1,3][[1]]@MethMotif_x@MMmotif@nPeaks)  )
		transfac_vector <- c( transfac_vector, c("XX","//")  )						  
						  
						  
    writeLines(transfac_vector, fileConn)
    close(fileConn)
