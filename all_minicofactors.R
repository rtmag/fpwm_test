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
cp_res_peaks[,2] <- cp_res_peaks[,2]-50
cp_res_peaks[,3] <- cp_res_peaks[,3]+50	
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

######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
# Bin TFBS p-values
pval_in_bins <- function(ft_file, pssm_name, pval_thr = 0.001){
  # Create grep cmd
  if(str_detect(ft_file, ".gz")){
    grep_cmd <- str_c(
      "gzcat ", ft_file, 
      "| grep -e 'site\t",pssm_name,"\t' ",  # Need to come back
      "| cut -f4,5,6,9"
    )
  } else{
    grep_cmd <- str_c(
      "grep -e 'site\\t",pssm_name,"\\t' ",ft_file,  # Need to come back
      "| cut -f4,5,6,9"
    )   
  }
  message(grep_cmd)
  # Load feature file
  ft_out <- grep_cmd %>% 
    pipe() %>% 
    read.table(
      sep = "\t", 
      stringsAsFactors = FALSE
    )
  
  # Filter  p-values
  ft_filt <- ft_out %>% 
    dplyr::filter(V4 <= pval_thr)
  
  # Convert to -log10(p-val)
  ft_filt <- ft_filt %>% 
    dplyr::mutate( V4 = (-log10(V4)) ) %>%
    rename(
      strand = V1 , 
      start = V2 ,
      end = V3 , 
      pval = V4 
    )
  
  # Create center pos
  ft_filt <- ft_filt %>%
    dplyr::mutate( pos = (start + end)/2 )
  
  # Creat bins and add name
  ft_filt <- ft_filt %>%
    dplyr::mutate( 
      binMap = pos %>%
        cut(
          breaks = seq(50,150, by = 10),
          labels = seq(51,150, by = 10)
          ),
      matrix_name = pssm_name,
    )
  
  # Filter out NA values
  ft_filt <- ft_filt %>% 
    dplyr::filter( !(is.na( binMap )  ) )

}

##

# Load feature file 1
pssm1_df <- pval_in_bins(ft_file = "/Users/wone/CSI/FPWM_walter/matrixScan_minicof/ft/matrix-scan_HepG2_USF1_AND_TFE3.ft", 
                         pssm_name = "MM1_HSA_HepG2_USF1_AND_MM1_HSA_HepG2_USF1")
# Load feature file 2
pssm2_df <- pval_in_bins(ft_file = "/Users/wone/CSI/FPWM_walter/matrixScan_minicof/ft/matrix-scan_HepG2_USF1_AND_TFE3.ft", 
                         pssm_name = "MM1_HSA_HepG2_USF1_AND_MM1_HSA_HepG2_TFE3")

pssm1_df <- pssm1_df %>% mutate( matrix_name = "USF1"  )
pssm2_df <- pssm2_df %>% mutate( matrix_name = "USF1+TFE3" )


databox2 <- rbind(pssm1_df,pssm2_df) %>% 
  mutate( matrix_name =  matrix_name %>% as.factor() )

p2 <- databox2 %>% 
  ggplot(aes(x = binMap, y = pval, fill = matrix_name)) + 
  geom_boxplot(outlier.shape = NA) + 
  xlab("JUND and FOSL2 overlapping peaks bins") + 
  ylab("-Log10 Pvalues") + 
  labs(fill='Motifs') +
  scale_x_discrete(
    labels=c(
      "51" = "-50 -40",
      "61" = "-40 -30",
      "71" = "-30 -20",
      "81" = "-20 -10",
      "91" = "-10 0",
      "101" = "0 +10",
      "111" = "+10 +20",
      "121" = "+20 +30",
      "131" = "+30 +40",
      "141" = "+40 +50")
    ) +
  theme_bw()

# Plot
p2
# Save as PDF
ggsave(filename = "../img/matrix_scan_pval_JUND_and_FOSL2_peaks.pdf")

# Get number of TFBSs
databox2 %>% 
  mutate(matrix_name = matrix_name %>% 
           factor(levels = c("JUND+ATF2","JUND", "JUND+FOSL2")
         )) %>% 
  ggplot( aes(x=matrix_name, fill = matrix_name) ) +
  geom_bar(stat="count", width=0.5)+
  labs(fill='Motifs') +
  theme_bw()
