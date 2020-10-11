library(TFregulomeR)
library(forkedTF)

all_record <- dataBrowser()

mm_table <- all_record[all_record$source=="MethMotif",]

mm_table$cell_tissue_name<- gsub(" ","",mm_table$cell_tissue_name)

for( i in 556:dim(mm_table)[1] ){
  miniCofactorReport( TF = mm_table$TF[i], cell = mm_table$cell_tissue_name[i],cobinding_threshold=.1 )
}

####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
matrixScan_data_gen<-function(main="",cofact=""){
  
    cebpb_inter <- intersectPeakMatrix(peak_id_x = main,
                                              peak_id_y = c(main,cofact),
                                              motif_only_for_id_y = TRUE,
                                              methylation_profile_in_narrow_region = TRUE)
  
    fileConn <- file(paste0(main,"_AND_",cofact,".transfac"))
  
		transfac_vector <- c()
    tmp_id <- gsub("\\[|\\]","",cebpb_inter[1,1][[1]]@id)
		transfac_vector <- c( transfac_vector, paste0("AC ",tmp_id) )
		transfac_vector <- c( transfac_vector, "XX" )
		transfac_vector <- c( transfac_vector, paste0("ID ",tmp_id), "XX" )
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

    tmp_id <- gsub("\\[|\\]","",cebpb_inter[1,2][[1]]@id)
		transfac_vector <- c( transfac_vector, paste0("AC ",tmp_id) )
		transfac_vector <- c( transfac_vector, "XX" )
		transfac_vector <- c( transfac_vector, paste0("ID ",tmp_id), "XX" )
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
  
		writeLines(transfac_vector, fileConn)
    close(fileConn)
  
  
cp_res <- commonPeaks(target_peak_id=main,
                                       motif_only_for_target_peak=TRUE,
                                       compared_peak_id=cofact,
                                       motif_only_for_compared_peak=TRUE,
                                       methylation_profile_in_narrow_region=FALSE)
cp_res_peaks <-cp_res[,1][[1]]@common_peak # 6325
tmp_title <- paste0(main,"_AND_",cofact,"_",dim(cp_res_peaks)[1],"_common_peaks.bed")
write.table(cp_res_peaks[,1:3],tmp_title,sep="\t",quote=F,row.names=F,col.names=F)
}

###  TEST CASES
# HepG2 USF2 cofactor: TFE3
matrixScan_data_gen(main="MM1_HSA_HepG2_USF2",cofact="MM1_HSA_HepG2_TFE3")
# HepG2 USF1 cofactor: TFE3
matrixScan_data_gen(main="MM1_HSA_HepG2_USF1",cofact="MM1_HSA_HepG2_TFE3")
# HCT116 ELF1 cofactor: FOSL1
matrixScan_data_gen(main="MM1_HSA_HCT116_ELF1",cofact="MM1_HSA_HCT116_FOSL1")
# A549 SREBF1 cofactor: JUNB
matrixScan_data_gen(main="MM1_HSA_A549_SREBF1",cofact="MM1_HSA_A549_JUNB")
