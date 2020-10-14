library(TFregulomeR)
library(forkedTF)
library(stringr)
library(dplyr)
library(ggplot2)

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

##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################
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
	
	ft_filt$pos <- ft_filt$pos+100
  
  # Creat bins and add name
  ft_filt <- ft_filt %>%
    dplyr::mutate( 
      binMap = pos %>%
        cut(
          breaks = seq(1,101, by = 10),
          labels = seq(51,150, by = 10)
          ),
      matrix_name = pssm_name,
    )
  
  # Filter out NA values
	return(ft_filt)

}

##############################################################################################################################################################################
plotMS <- function(df1,df2,mat1,mat2,colors){
	pssm1_df <- df1 %>% mutate( matrix_name = mat1  )
	pssm2_df <- df2 %>% mutate( matrix_name = mat2 )

	databox2 <- rbind(pssm1_df,pssm2_df) %>% 
	  mutate( matrix_name =  matrix_name %>% as.factor() )

	p2 <- databox2 %>% 
	  ggplot(aes(x = binMap, y = pval, fill = matrix_name)) + 
	  geom_boxplot(outlier.shape = NA) + 
	  xlab(paste0(mat2," peaks center")) + 
	  ylab("-Log10 Pvalues") + 
	  labs(fill='Motifs') +
	  scale_x_discrete(
	    labels=c(
	      "51" = "-50",
	      "61" = "-40",
	      "71" = "-30",
	      "81" = "-20",
	      "91" = "-10",
	      "101" = "+10",
	      "111" = "+20",
	      "121" = "+30",
 	     "131" = "+40",
	      "141" = "+50")
	    ) +
 	 theme_bw() + scale_fill_manual(values=colors)
	
	return(p2)
	}
##############################################################################################################################################################################
# Load feature file 1
pssm1_df <- pval_in_bins(ft_file = "/Users/wone/CSI/FPWM_walter/matrixScan_minicof/ft/matrix-scan_HepG2_USF1_AND_TFE3.ft", 
                         pssm_name = "MM1_HSA_HepG2_USF1_AND_MM1_HSA_HepG2_USF1",pval_thr=.01)
# Load feature file 2
pssm2_df <- pval_in_bins(ft_file = "/Users/wone/CSI/FPWM_walter/matrixScan_minicof/ft/matrix-scan_HepG2_USF1_AND_TFE3.ft", 
                         pssm_name = "MM1_HSA_HepG2_USF1_AND_MM1_HSA_HepG2_TFE3",pval_thr=.01)

pp <- plotMS(pssm1_df,pssm2_df,mat1="USF1",mat2="USF1+TFE3")
pp
ggsave(filename = "USF1+TFE3_HepG2_matrixScan_pval.pdf")

table(databox2$matrix_name)
##############################################################################################################################################################################
pssm1_df <- pval_in_bins(ft_file = "/Users/wone/CSI/FPWM_walter/matrixScan_minicof/ft/matrix-scan_HCT116_ELF1_AND_FOSL1.ft", 
                         pssm_name = "MM1_HSA_HCT116_ELF1_AND_MM1_HSA_HCT116_ELF1",pval_thr=.01)
# Load feature file 2
pssm2_df <- pval_in_bins(ft_file = "/Users/wone/CSI/FPWM_walter/matrixScan_minicof/ft/matrix-scan_HCT116_ELF1_AND_FOSL1.ft", 
                         pssm_name = "MM1_HSA_HCT116_ELF1_AND_MM1_HSA_HCT116_FOSL1",pval_thr=.01)

pp <- plotMS(pssm1_df,pssm2_df,mat1="ELF1",mat2="ELF1+FOSL1")
pp
ggsave(filename = "ELF1+FOSL1_HCT116_matrixScan_pval.pdf")
table(databox2$matrix_name)
##############################################################################################################################################################################
pssm1_df <- pval_in_bins(ft_file = "/Users/wone/CSI/FPWM_walter/matrixScan_minicof/ft/matrix-scan_HepG2_USF2_AND_TFE3.ft", 
                         pssm_name = "MM1_HSA_HepG2_USF2_AND_MM1_HSA_HepG2_USF2",pval_thr=.01)
# Load feature file 2
pssm2_df <- pval_in_bins(ft_file = "/Users/wone/CSI/FPWM_walter/matrixScan_minicof/ft/matrix-scan_HepG2_USF2_AND_TFE3.ft", 
                         pssm_name = "MM1_HSA_HepG2_USF2_AND_MM1_HSA_HepG2_TFE3",pval_thr=.01)

pp <- plotMS(pssm1_df,pssm2_df,mat1="USF2",mat2="USF2+TFE3")
pp
ggsave(filename = "USF2+TFE3_HepG2_matrixScan_pval.pdf")
table(databox2$matrix_name)
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
# #F8766D (red) #619CFF (blue)  #00BA38 (green) 
###################
xx <- read.table(pipe('grep -w "site" matrix-scan_atf2peaks.ft |cut -f 3,9'))
bartable <- table(xx[xx[,2]<.0001,1])
names(bartable) <- c("JUND","JUND+ATF2","JUND+FOSL2")

pdf("matrix-scan_atf2peaks_barplot.pdf")
barplot(bartable,main="JUND+ATF2 peaks",ylab="Number of predicted binding sites", col=c("#F8766D","#619CFF","#00BA38"),cex.axis=1.5,cex.lab=1.5,cex.names=1.5,cex.main=1.5)
dev.off()

pssm1_df <- pval_in_bins(ft_file = "matrix-scan_atf2peaks.ft", 
                         pssm_name = "JUND",pval_thr=.01)
pssm2_df <- pval_in_bins(ft_file = "matrix-scan_atf2peaks.ft", 
                         pssm_name = "MM1_HSA_HepG2_JUND_overlapped_with_MM1_HSA_HepG2_ATF2",pval_thr=.01)
pp <- plotMS(pssm1_df,pssm2_df,mat1="JUND",mat2="JUND+ATF2",colors=c("#F8766D","#619CFF"))
ggsave(filename = "JUND+ATF2_HepG2_matrixScan_pval.pdf",height=3,width=4)
###################
xx <- read.table(pipe('grep -w "site" matrix-scan_fosl2peaks.ft |cut -f 3,9'))
bartable <- table(xx[xx[,2]<.0001,1])
names(bartable) <- c("JUND","JUND+ATF2","JUND+FOSL2")

pdf("matrix-scan_fosl2peaks_barplot.pdf")
barplot(bartable,main="JUND+FOSL2 peaks",ylab="Number of predicted binding sites", col=c("#F8766D","#619CFF","#00BA38"),cex.axis=1.5,cex.lab=1.5,cex.names=1.5,cex.main=1.5)
dev.off()

pssm1_df <- pval_in_bins(ft_file = "matrix-scan_fosl2peaks.ft", 
                         pssm_name = "JUND",pval_thr=.01)
pssm2_df <- pval_in_bins(ft_file = "matrix-scan_fosl2peaks.ft", 
                         pssm_name = "MM1_HSA_HepG2_JUND_overlapped_with_MM1_HSA_HepG2_FOSL2",pval_thr=.01)
pp <- plotMS(pssm1_df,pssm2_df,mat1="JUND",mat2="JUND+FOSL2",colors=c("#F8766D","#00BA38"))
ggsave(filename = "JUND+FOSL2_HepG2_matrixScan_pval.pdf",height=3,width=4)
